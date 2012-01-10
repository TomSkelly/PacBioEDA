#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Python script to print summary information about the contents of a
# *.bas.h5 file, as produced by the PacBio instrument. Provides counts
# of reads, subreads and bases for a hierarchy of increasingly
# restrictive subsets of subreads:

#       Total contents of bas file
#         Sequencing ZMWs
#           Productivity-0 Sequencing ZMWs
#           Productivity-1 Sequencing ZMWs
#             Reads with HQ region length > threshold
#               Reads with HQ score > threshold
#                 Reads with average insert size > threshold (i.e., not adapter dimers)
#                   Reads which aligned (if a .cmp.h5 file was supplied)
#           Productivity-2 Sequencing ZMWs

# All of the above except average insert size are criteria applied by
# the filtering step of secondary analysis.

# Inputs are a .bas.h5 file and (optionally) a .cmp.h5 file, specified
# as command line parameters. Filtering parameters can be overridden
# on the command line. See the --help option for details. Output is to
# stdout.

import sys
import optparse
import H5BasFile
import H5CmpFile
from tt_log import logger

DEF_SCORE_THRESHOLD  = 750
DEF_HQ_LENGTH        = 50
DEF_INSERT_THRESHOLD = 100

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    basFilename = args[0]
    logger.debug("bas file: %s" % basFilename)
    bf = H5BasFile.BasFile (basFilename)

    cf = None                      # no cmp file?

    if len(args) > 1:
        
        cmpFilename = args[1]
        logger.debug("cmp file: %s" % cmpFilename)
        cf = H5CmpFile.CmpFile (cmpFilename,
                                set=bf.setNumber(), 
                                strobe=bf.strobeNumber(),
                                maxHole=bf.numZMWs())

    totalC   = Counter('Total')
    seqC     = Counter('--Sequencing')
    prod0C   = Counter('----Productivity-0')
    prod1C   = Counter('----Productivity-1')
    HQLenC   = Counter('------HQ Len >= %s' % opt.length)
    HQScoreC = Counter('--------HQ Score >= %s' % opt.score)
    adaptC   = Counter('----------Insert >= %s' % opt.insert)
    alignC   = Counter('------------Aligned')
    prod2C   = Counter('----Productivity-2')

    longest    = 0
    longestZMW = None

    for hole in xrange(bf.numZMWs()):

        numBases = bf.readLen(hole)
        zProd    = bf.productivity(hole)

        if numBases > longest:
            longest    = numBases
            longestZMW = hole

        HQStart, HQEnd, HQScore = bf.HQregion(hole)[2:5]
        HQLen = HQEnd - HQStart

        numSubreads     = 0
        maxSubreadLen   = 0
        cumSubreadLen   = 0
        alignedSubreads = 0
        alignedBases    = 0

        for region in bf.holeRegions(hole):

            regionHole, regionType, start, end, score = region
            if regionType == 1:                   # if insert
                numSubreads += 1
                maxSubreadLen  = max (end-start, maxSubreadLen)
                cumSubreadLen += max (end-start, 0)   # clip negative lengths to zero

                if cf is not None:
                    align = cf.getAlignmentAsDict (hole, start, end) # alignment record for this region

                    if align is not None:                            # if the region aligned
                        alignedSubreads += 1
                        alignedBases    += align['rEnd'] - align ['rStart']

        # What follows is a series of increasingly restrictive
        # criteria for a useful subread. Keep track of the number of
        # ZMWs, subreads, and bases which pass the successive tests,
        # and the length and ZMW provenance of the longest accepted
        # subread.

        totalC.incr (1, numSubreads, numBases)
        totalC.longest (hole, maxSubreadLen)

        if bf.isSequencingZMW(hole):              # sequencing ZMW?
            seqC.incr (1, numSubreads, numBases)
            seqC.longest (hole, maxSubreadLen)

            if zProd == 0:
                prod0C.incr (1, numSubreads, numBases)
                prod0C.longest (hole, maxSubreadLen)

            elif zProd == 2:
                prod2C.incr (1, numSubreads, numBases)
                prod2C.longest (hole, maxSubreadLen)

            elif zProd == 1:                      # productivity 1 gets broken down further
                prod1C.incr (1, numSubreads, numBases)
                prod1C.longest (hole, maxSubreadLen)

                if HQLen >= opt.length:
                    HQLenC.incr (1, numSubreads, numBases)
                    HQLenC.longest (hole, maxSubreadLen)

                    if HQScore >= opt.score:
                        HQScoreC.incr (1, numSubreads, numBases)
                        HQScoreC.longest (hole, maxSubreadLen)

                        # A very short average insert size probably indicates an adapter dimer.

                        if cumSubreadLen >= numSubreads * opt.insert:
                            adaptC.incr (1, numSubreads, numBases)
                            adaptC.longest (hole, maxSubreadLen)

                            if alignedSubreads > 0:
                                alignC.incr (1, alignedSubreads, alignedBases)
                                alignC.longest (hole, alignedBases)

    print
    print "file: ", basFilename
    print
    print "longest read was ZMW %d at %d bases" % (longestZMW, longest)
    print
    print "statistics for subreads:"
    print

    Counter.title();

    if cf is not None:                 # if we processed a .cmp.h5 file
        for cntr in (totalC, seqC, prod0C, prod1C, HQLenC, HQScoreC, adaptC, alignC, prod2C):
            cntr.longPrint()
    else:
        for cntr in (totalC, seqC, prod0C, prod1C, HQLenC, HQScoreC, adaptC, prod2C):
            cntr.longPrint()
    print

    logger.debug("complete")

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file> [<cmp_file>]')

    parser.add_option ('--score',  type='int', help='Minimum HQ region score (def: %default)')
    parser.add_option ('--length', type='int', help='Minimum HQ region length (def: %default)')
    parser.add_option ('--insert', type='int', help='Minimum average insert length (def: %default)')

    parser.set_defaults (score=DEF_SCORE_THRESHOLD,
                         length=DEF_HQ_LENGTH,
                         insert=DEF_INSERT_THRESHOLD)

    opt, args = parser.parse_args()

    return opt, args

class Counter (object):

    def __init__ (self, name):

        self._name        = name
        self._numReads    = 0
        self._numSubreads = 0
        self._numBases    = 0
        self._whichZMW    = None
        self._maxLen      = 0          # WARNING: assumes value >= 0

    def incr (self, reads, regions, bases):

        self._numReads    += reads
        self._numSubreads += regions
        self._numBases    += bases

    def longest (self, which, value):

        if value > self._maxLen:
            self._maxLen   = value
            self._whichZMW = which

    @staticmethod
    def title ():

        print "reads subreads      bases     max    ZMW  criterion"
        print

    def longPrint (self):

        print "%5d  %7d  %9d  %6d  %5d  %s" \
            % (self._numReads, self._numSubreads, self._numBases, self._maxLen, self._whichZMW, self._name)

if __name__ == "__main__":
    main()

# Copyright (C) 2011 Genome Research Limited
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
