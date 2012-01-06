#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Read a PacBio bas.h5 file, extract the sequence for a ZMW specified
# on the command line, align it to the adapter sequence, and output a
# plot of Smith-Waterman alignment scores at each position in the
# read. By default, output is to file ZMW-zzzzz.png, where zzzzz is
# the zero-padded ZMW number -- this can be overridden on the command
# line. Optional command line parameters specify a 0-based inclusive
# range of bases within the read.

import sys
import optparse
import H5BasFile
import SWAligner
from tt_log import logger

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    basFilename = args[0]
    basfile = H5BasFile.BasFile (basFilename)

    aln = SWAligner.Aligner()
    aln.setRead (H5BasFile.ADAPTER)

    try:
        hole = int(args[1])
    except ValueError:
        logger.error('ERROR: second parameter must be an integer ZMW number')
        sys.exit()

    start = opt.start
    end   = basfile.readLen(hole) if opt.end is None else opt.end

    sequence = basfile.getSequence(hole, start, end)

    aln.setRef (sequence)

    score = aln.fillMatrix()

    refString, readString = aln.alignmentStrings()

    print "best alignment score: %d" % score
    print refString
    print readString
    print

    allScores = aln.allScores()

    if opt.output is not None:
        outfile = opt.output
    else:
        outfile = "ZMW-%05d.png" % hole

    graphScores (xrange(start,end), 
                 allScores, 
                 outfile,
                 "ZMW %d (%d to %d)" % (hole, start, end),
                 len(H5BasFile.ADAPTER))

    logger.debug("complete")

def graphScores (range, scores, outfile, title=None, ymax=None):

    fig = plt.figure()

    if title != None:
        fig.suptitle(title, fontsize=14, fontweight='bold')

    ax = fig.add_subplot(1,1,1)
    ax.plot (range, scores)
    if ymax != None:
        plt.ylim (0, ymax)

    plt.savefig (outfile)

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file> <ZMW#>')

    parser.add_option ('--start', type='int', help='0-based start offset of output into read (def: %default)')
    parser.add_option ('--end',   type='int', help='0-based end offset of output into read (def: all)')
    parser.add_option ('--output',            help='Output file name (def: %default)')

    parser.set_defaults (start=0,
                         end=None,
                         output=None)

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
