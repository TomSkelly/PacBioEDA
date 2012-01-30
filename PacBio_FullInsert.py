#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Read a PacBio bas.h5 file, find 'full' insert regions, and produce a
# histogram of their lengths. Full inserts are insert regions which
# cover the entire sequenced milecule. I.e., inserts which run from
# one adapter to the other.

import sys
import optparse
import numpy as np                # for argmax
import H5BasFile
from tt_log import logger

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEF_BINS   = 100
DEF_MIN    = 20
DEF_MAX    = 999999999            # impossibly large
DEF_OUTPUT = 'full_inserts.png'

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    basFilename = args[0]
    bf = H5BasFile.BasFile (basFilename)

    sizes          = []
    totInserts     = 0
    shortInserts   = 0
    longInserts    = 0
    partialInserts = 0

    for hole in xrange(bf.numZMWs()):

        if bf.isSequencingZMW(hole):                    # sequencing ZMW?
            if bf.productivity(hole) == 1:              # prod-1?

                # We need to find inserts which lie between
                # adapters. since inserts come before adapters in the
                # bas.h5 file, we'll have to scan the regions list
                # twice. (We could save the insert information, but
                # then we'd have to scan that list.)

                firstAdapt = 999999999
                lastAdapt  = 0

                for region in bf.holeRegions(hole):     # loop through the regions looking for adapters

                    regionType, regionStart, regionEnd = region[1:4]

                    if regionType == 0:                 # an adapter?

                        firstAdapt = min(regionEnd,  firstAdapt)
                        lastAdapt  = max(regionStart, lastAdapt)
                
                for region in bf.holeRegions(hole):     # next loop through the regions looking for inserts

                    regionType, regionStart, regionEnd = region[1:4]

                    if regionType == 1:                 # an insert?

                        totInserts += 1

                        if regionStart >= firstAdapt and regionEnd <= lastAdapt:

                            regionSize = regionEnd - regionStart
                            if regionSize < opt.min:
                                shortInserts += 1
                            elif regionSize > opt.max:
                                longInserts += 1
                            else:
                                sizes.append(regionSize)

                        else:
                            partialInserts += 1
                
    n, bins, patches = plt.hist(sizes, bins=opt.bins)

    plt.savefig (opt.output)

    maxBin = np.argmax(n)              # find the fullest bin

    logger.debug("found %d productivity-1 inserts" % totInserts)
    logger.debug("ignored %d as partial"           % partialInserts)
    logger.debug("ignored %d as too short"         % shortInserts)
    logger.debug("ignored %d as too long"          % longInserts)
    logger.debug("plotted %d"                      % (totInserts-partialInserts-shortInserts-longInserts))
    logger.debug("bin size is %7.2f"               % (bins[1]-bins[0]))
    logger.debug("max bin was %7.1f to %7.1f"      % (bins[maxBin], bins[maxBin+1]))

    if opt.details:
        logger.debug("details:")
        for ix in xrange(len(n)):      # beware: len(bins) == len(n)+1
            print "%7.1f %5d" % (bins[ix], n[ix])

    logger.debug("complete")

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file>')

    parser.add_option ('--bins',   type='int', help='number of bins in histogram (def: %default)')
    parser.add_option ('--min',    type='int', help='minimum insert size to accept (def: %default)')
    parser.add_option ('--max',    type='int', help='maximum insert size to accept (def: any)')
    parser.add_option ('--output',             help='output file name (def: %default)')
    parser.add_option ('--details', action='store_true', help='show me the numbers')

    parser.set_defaults (bins=DEF_BINS,
                         min=DEF_MIN,
                         max=DEF_MAX,
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
