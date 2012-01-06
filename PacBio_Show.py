#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Read a PacBio bas.h5 file, print a range of sequence from a
# specified ZMW in fasta format. Command line parameters are bas.h5
# file name and ZMW#. Options include 0-based start and stop offsets.

import sys
import optparse
import H5BasFile
from tt_log import logger

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    basFilename = args[0]
    logger.debug("bas file: %s" % basFilename)
    basfile = H5BasFile.BasFile (basFilename)

    try:
        hole = int(args[1])
    except ValueError:
        logger.error('ERROR: second parameter must be an integer ZMW number')
        sys.exit()

    if not opt.reverse:
        sequence = basfile.getSequence(hole, opt.start, opt.end)    # end==None gets the whole read
    else:
        sequence = basfile.getRevCompSequence(hole, opt.start, opt.end)

    movie = basfile.movieName()
    length = len(sequence)
    print ">%s/%d/%d_%d" % (movie, hole, opt.start, opt.start+length)

    for ix in xrange(0,length,opt.flen):
        print sequence[ix:ix+opt.flen]

    logger.debug("complete")

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file> <ZMW#>')

    parser.add_option ('--start',     type='int', help='0-based start offset of output into read (def: %default)')
    parser.add_option ('--end',       type='int', help='0-based end offset of output into read (def: all)')
    parser.add_option ('--fasta-len', type='int', help='Number of bases in output fasta line (def: %default)', dest='flen')
    parser.add_option ('--reverse', action='store_true', help='output reverse-complement of region')

    parser.set_defaults (start=0,
                         end=None,
                         reverse=False,
                         flen=60)

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
