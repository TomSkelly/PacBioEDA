#!/usr/bin/env python

# Create a time-stepped series of plots showing the spatial
# distribution on the SMRTcell of ZMWs in their HQ regions.

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

import sys
import optparse
import H5BasFile
from tt_log import logger

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEF_START = 0
DEF_STEP  = 300        # 300 * 6 * 4 = 7200 = a 2-hour movie
DEF_NROWS = 6
DEF_NCOLS = 4
DEF_OUTPUT='HQ-history.png'

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    basFilename = args[0]
    logger.debug("bas file: %s" % basFilename)
    bf = H5BasFile.BasFile (basFilename)

    HQTimes = []

    for hole in bf.holeNumbers():       # find start and end times of HQ region for each ZMW

        if bf.isSequencingZMW(hole) and bf.productivity(hole) == 1:

            HQStart, HQEnd = bf.HQregion(hole)[2:4]

            if HQEnd > 0:
                start = float(bf.elapsedFrames(hole, 0, HQStart)) / H5BasFile.frameRate
                end   = float(bf.elapsedFrames(hole, 0, HQEnd))   / H5BasFile.frameRate
                HQTimes.append((hole,start,end))

    logger.debug("found %d sequencing prod-1 ZMWs" % len(HQTimes))

    coords = bf.cellCoords()

    plt.figure(figsize=(opt.ncols*3,opt.nrows*3))
    plotno = 0
    plotTime = float(opt.start)

    for row in xrange(opt.nrows):
        for col in xrange(opt.ncols):

            logger.debug("time: %6.1f" % plotTime)

            scatterX   = []
            scatterY   = []
            scatterCol = []
            numPreHQ   = 0
            numInHQ    = 0
            numPostHQ  = 0

            for hole,start,end in HQTimes:

                X,Y = bf.holeXY(hole)
                scatterX.append(X)
                scatterY.append(Y)

                if start > plotTime:         # if HQ has not yet started
                    scatterCol.append('yellow')
                    numPreHQ += 1
                elif end > plotTime:         # if HQ has not yet ended
                    scatterCol.append('green')
                    numInHQ += 1
                else:                        # else HQ has ended
                    scatterCol.append('red')
                    numPostHQ += 1

            plotno += 1
            ax = plt.subplot(opt.nrows, opt.ncols, plotno)

            ax.scatter (scatterX, scatterY, c=scatterCol, s=1, edgecolor='face')

            plt.title("Time: %5.0f  %d/%d/%d" % (plotTime, numPreHQ, numInHQ, numPostHQ),
                      fontsize='x-small')
            ax.axis(coords)
            ax.set_xticklabels([])
            ax.set_yticklabels([])

            plotTime += float(opt.step)


    if opt.title is not None:
        plt.suptitle(opt.title)

    plt.savefig (opt.output)

    logger.debug("complete")

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file>')

    parser.add_option ('--start', type='int', help='Time of first plot in seconds (def: %default)')
    parser.add_option ('--step',  type='int', help='Time interval between plots in seconds (def: %default)')
    parser.add_option ('--nrows', type='int', help='Number of plots per row (def: %default)')
    parser.add_option ('--ncols', type='int', help='Number of plots per column (def: %default)')
    parser.add_option ('--title',             help='Title for top of figure (def: %default)')
    parser.add_option ('--output',            help='Output file name (def: %default)')

    parser.set_defaults (start=DEF_START,
                         step=DEF_STEP,
                         nrows=DEF_NROWS,
                         ncols=DEF_NCOLS,
                         title=None,
                         output=DEF_OUTPUT)

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
