#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Print details of sub-read alignments for one movie from a cmp.h5
# file. Alignments are printed in the order in which they appear in
# the cmp file, which may have been sorted in several ways. Note that
# one of the corresponding bas.h5 files must be specified as the first
# command line parameter, to determine the desired movie.

# See the --help option for details. Output is to stdout.

import sys
import optparse
import H5BasFile
import H5CmpFile
from tt_log import logger

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    if len(args) != 2:
        logger.error ("please specify bas.h5 and cmp.h5 files as parameters. See --help")
        sys.exit()

    basFilename = args[0]
    logger.debug("bas file: %s" % basFilename)
    bf = H5BasFile.BasFile (basFilename)

    cmpFilename = args[1]
    logger.debug("cmp file: %s" % cmpFilename)
    cf  = H5CmpFile.CmpFile (fileName=cmpFilename)
    cmp = H5CmpFile.CmpMovie (cmpObject=cf,
                              movieName=bf.movieName(),
                              maxHole=bf.maxZMW())

    cf.printDetails()

    print " AlnID  RG   Hole Set Stb  SubRd Seq Ref St      Start        End RefStrt  RefEnd    OffStrt     OffEnd"
    print

    if opt.sort == 'hole':

        for align in cmp.getAlignmentsByHole():
            printAlign(align)

    else :                                                    # else, must be 'none' 

        for align in cmp.getAllAlignments():                  # generator function, returns a dict
            printAlign(align)

    logger.debug("complete")

def printAlign (align):

    print "%6d  %2d  %5d  %2d  %2d  %5d  %2d  %2d  %1d  %9d  %9d  %6d  %6d  %9d  %9d   %5d  %5d" \
        % (align['AlignmentId'],
           align['ReadGroupId'],
           align['HoleNumber'],
           align['SetNumber'],
           align['StrobeNumber'],
           align['SubreadId'],
           align['RefSeqId'],
           align['contig'],                # chicanery here, see H5CmpFile.py
           align['RCRefStrand'],
           align['tStart'],
           align['tEnd'],
           align['rStart'],
           align['rEnd'],
           align['offset_begin'],
           align['offset_end'],
           align['rEnd'] - align['rStart'],
           align['tEnd'] - align['tStart'])

    return

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file> <cmp_file>',
                                   description='Print (to stdout) summary information about the contents of a cmp.h5 file.')

    parser.add_option ('--sort', type='choice', choices=['hole', 'none'],
                       help='desired sort order: none/hole (def: %default)')

    parser.set_defaults (sort='none')

    opt, args = parser.parse_args()

    return opt, args

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
