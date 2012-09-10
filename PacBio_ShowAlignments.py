#!/usr/bin/env python

# Copyright (C) 2012 Genome Research Limited -- See full notice at end
# of module.

# Print sub-read alignment strings for one movie from a cmp.h5
# file. Note that one of the corresponding bas.h5 files must be
# specified as the first command line parameter, to determine the
# desired movie (but see the comment below).

# Printed offsets are zero-based.

# See the --help option for details. Output is to stdout.

import sys
import optparse
import H5BasFile
import H5CmpFile
from tt_log import logger

DEF_LINE_LEN = 60

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    if len(args) != 2:
        logger.error ("please specify bas.h5 and cmp.h5 files as parameters. See --help")
        sys.exit()

    # TODO: Actually, all we need from the bas file is the movie name
    # (maxHole will default to something clever). We don't need to
    # open the bas file to determine the movie name: It's part of the
    # filename. The only real reason we specify a bas file as the
    # first parameter is to match the command line interface of other
    # scripts.

    basFilename = args[0]
    logger.debug("bas file: %s" % basFilename)
    bf = H5BasFile.BasFile (basFilename)
    movie = bf.movieName()

    cmpFilename = args[1]
    logger.debug("cmp file: %s" % cmpFilename)
    cf  = H5CmpFile.CmpFile (fileName=cmpFilename)
    cmp = H5CmpFile.CmpMovie (cmpObject=cf,
                              movieName=movie,
                              maxHole=bf.maxZMW())

    if opt.ZMW is not None:                      # did we ask for a specific ZMW?

        for align in cmp.getAlignmentsForHole(opt.ZMW):

            printAlignment (align, cmp, opt.flen)

    else:                                        # else, print all ZMWs

        for align in cmp.getAlignmentsByHole():

            printAlignment (align, cmp, opt.flen)

    logger.debug("complete")

def printAlignment (align, cmp, flen):

    movie  = cmp.movieName()
    hole   = align['HoleNumber']
    rStart = align['rStart']
    rEnd   = align['rEnd']

    print "%s/%d/%d_%d" % (movie, hole, rStart, rEnd);    # read name
    print

    refString, readString = cmp.getAlignmentStrings (align)

    refDashes  = 0
    readDashes = 0

    for ix in xrange(0,len(refString),flen):

        print "   %5d  " % (ix - refDashes),
        print refString[ix:ix+flen]
        print "   %5d  " % (ix - readDashes),
        print readString[ix:ix+flen]
        print

        refDashes  += refString.count('-', ix, ix+flen)
        readDashes += readString.count('-', ix, ix+flen)

    return

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file> <cmp_file>',
                                   description='Print (to stdout) alignment strings from a cmp.h5 file.')

    parser.add_option ('--ZMW',      type='int', help='ZMW to print (def: all ZMWs)')
    parser.add_option ('--line-len', type='int', help='Number of bases in output line (def: %default)',
                       dest='flen')

    parser.set_defaults (flen=DEF_LINE_LEN)

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
