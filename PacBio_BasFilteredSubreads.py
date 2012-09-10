#!/usr/bin/env python

# Copyright (C) 2012 Genome Research Limited -- See full notice at end
# of module.

# Python script to extract filtered subreads from a bas.h5 file, and
# write them to a fasta (or, optionally, fastq) file. By default, it
# will reproduce the sec/data/filtered_subreads.fasta file for the set
# described by the input .bas.h5 file. Different filtering parameters
# can be specified as command line parameters.

# The main reason for the existence of this script is to convince
# myself that I understand how filtering works. But it also provides a
# quick way to extract a set of reads from a .bas.h5 file.

# The option to produce a fastq file is not easily available elsewhere
# (as of SMRTanalysis 1.2.3). It uses the QualityValue field from the
# bas.h5 data, and ascii-encodes it as Q+33.

import sys
import optparse
import H5BasFile
from tt_log import logger

# Defaults for command line options:

DEF_SCORE_THRESHOLD  = 750
DEF_HQ_LENGTH        = 50
DEF_INSERT_THRESHOLD = 1          # default = no threshold, as done by primary analysis
DEF_FASTA_LEN        = 50

PHRED_SCALER = 33

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    basFilename = args[0]
    logger.debug("bas file: %s" % basFilename)
    bf = H5BasFile.BasFile (basFilename)

    numReads = 0
    numBases = 0

    for hole in bf.holeNumbers():

        if bf.isSequencingZMW(hole):
            
            if bf.productivity(hole) == 1:
                
                HQStart, HQEnd, HQScore = bf.HQregion(hole)[2:5]
                HQLen = HQEnd - HQStart

                if HQLen >= opt.length:
                    if HQScore >= opt.score:

                        for region in bf.holeRegions(hole):

                            regionHole, regionType, regionStart, regionEnd, regionScore = region

                            if regionType == 1:                   # if insert

                                start = max(HQStart, regionStart)
                                end   = min(HQEnd,   regionEnd)

                                if end - start >= opt.insert:

                                    if opt.fastq:
                                        printFastq (bf, hole, start, end)
                                    else:
                                        printFasta (bf, hole, start, end, opt.flen)

                                    numReads += 1
                                    numBases += end-start

    logger.debug("wrote %d reads, %d bases" % (numReads, numBases))
    logger.debug("complete")

def getCleanedSequence (bf, hole, start, end):
    '''Get the sequence between start and end and return a fiddled version of it. Not currently used.'''

    sequence = bf.getBasecallField("Basecall",      hole, start, end)
    insertQ  = bf.getBasecallField("InsertionQV",   hole, start, end)
    duration = bf.getBasecallField("WidthInFrames", hole, start, end)

    intBases = []

    for ix in xrange(end-start):
        if duration[ix] >= 4 or insertQ[ix] >= 5:
            intBases.append(chr(sequence[ix]))

    return ''.join(intBases)

def printFasta (bf, hole, start, end, flen):

    movie = bf.movieName()
    print ">%s/%d/%d_%d" % (movie, hole, start, end)

####    sequence = getCleanedSequence(bf, hole, start, end)
    sequence = bf.getSequence(hole, start, end)
    for ix in xrange(0,len(sequence),flen):
        print sequence[ix:ix+flen]

def printFastq (bf, hole, start, end):
    '''Print sequence and qualities in fastq format.'''

    movie = bf.movieName()
    print "@%s/%d/%d_%d" % (movie, hole, start, end)

    sequence = bf.getSequence(hole, start, end)
    print sequence
    print '+'

    qualities = bf.getBasecallField("QualityValue", hole, start, end)
    phredQuals = ''.join(chr(Q+PHRED_SCALER) for Q in qualities)
    print phredQuals

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file>',
                                   description='Output (to stdout) filtered subreads in fasta or fastq format.')

    parser.add_option ('--score',     type='int', help='Minimum HQ region score (def: %default)')
    parser.add_option ('--length',    type='int', help='Minimum HQ region length (def: %default)')
    parser.add_option ('--insert',    type='int', help='Minimum insert length (def: %default)')
    parser.add_option ('--fastq',     action='store_true', help='Output data as fastq file (def: fasta)')
    parser.add_option ('--fasta-len', type='int', help='Number of bases in output fasta line (def: %default)',
                                      dest='flen')

    parser.set_defaults (score=DEF_SCORE_THRESHOLD,
                         length=DEF_HQ_LENGTH,
                         insert=DEF_INSERT_THRESHOLD,
                         flen=DEF_FASTA_LEN)

    opt, args = parser.parse_args()

    return opt, args

if __name__ == "__main__":
    main()

# Copyright (C) 2012 Genome Research Limited
#
# This library is free software. You can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
