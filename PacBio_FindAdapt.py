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

# By default, different colours are used to indicate bases in the HQ
# region, and in regions called adapters by the analysis
# software. Colouring can be turned off by command line option, since
# finding the HQ region is a time-consuming process.

import sys
import optparse
import H5BasFile
import SWAligner
from tt_log import logger

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COL_NOT_HQ = '#00e000'     # colour for non-HQ region
COL_HQ     = '#408000'     # HQ region
COL_ADAPT  = '#e000e0'     # adapter

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    basFilename = args[0]
    basfile = H5BasFile.BasFile (basFilename)

    try:
        hole = int(args[1])
    except ValueError:
        logger.error('ERROR: second parameter must be an integer ZMW number')
        sys.exit()

    start = opt.start
    end   = basfile.readLen(hole) if opt.end is None else opt.end

    sequence = basfile.getSequence(hole, start, end)

    aln = SWAligner.Aligner()
    aln.setRef (sequence)
    aln.setRead (H5BasFile.ADAPTER)
    aln.fillMatrix()
    allScores = aln.allScores()

    range = xrange(start,end)
    title = "ZMW %d (%d to %d)" % (hole, start, end)

    plt.suptitle(title, fontsize=14, fontweight='bold')

    plt.plot (range, allScores, COL_NOT_HQ, zorder=1, label='non-HQ')

    if not opt.nocol:           # finding HQ region takes a long time, so optionally turn it off

        # There doesn't seem to be a way to separately specify a
        # colour for each point in a plot. So we'll plot in one
        # colour, then overlay subregions of that with another
        # colour. Plot commands are rendered in increasing zorder.

        HQStart, HQEnd = basfile.HQregion(hole)[2:4]
        HQRange = xrange(HQStart, HQEnd)
        HQScores = allScores[HQStart:HQEnd]

        plt.plot (HQRange, HQScores, COL_HQ, zorder=2, label='HQ')

        label = 'adapter';             # I will only say this once...

        for region in basfile.holeRegions(hole):     # loop through the regions looking for adapters

            regionType, regionStart, regionEnd = region[1:4]

            if regionType == 0:        # an adapter?
                regionRange = xrange(regionStart, regionEnd)
                regionScores = allScores[regionStart:regionEnd]
                plt.plot (regionRange, regionScores, COL_ADAPT, zorder=3, label=label)

                label = '_nolegend_'   # don't generate multiple legend entries

        plt.legend(loc='best', prop={'size':10})     # add a legend box to figure

    plt.ylim (0, len(H5BasFile.ADAPTER))

    if opt.output is not None:
        outfile = opt.output
    else:
        outfile = "ZMW-%05d.png" % hole
    plt.savefig (outfile)

    logger.debug("complete")

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file> <ZMW#>')

    parser.add_option ('--start', type='int', help='0-based start offset of output into read (def: %default)')
    parser.add_option ('--end',   type='int', help='0-based end offset of output into read (def: all)')
    parser.add_option ('--nocol', action='store_true', help='do not colour plot by region (faster!)')
    parser.add_option ('--output',            help='output file name (def: ZMW-nnnnn.png)')

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
