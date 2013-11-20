#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Create a time plot showing the number of ZMWs not yet in
# High-Quality region, in HQ, and post-HQ.

# The outer loop here will loop over ZMWs. For each ZMW, we'll first
# retrieve its HQ region start and stop base offsets into the
# read. We'll use the elapsedFrames method of H5BasFile to compute the
# number of frames in the HQ region (counting pulse durations and
# inter-pulse times), and use that time information to update our
# histogram.

# See also PacBio_HQHistory.py, which presents the same story as a
# time series of plots showing the spatial distribution of HQ ZMWs on
# the SMRTcell.

import sys
import optparse
import H5BasFile
from tt_log import logger

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEF_BINSIZE = 30
DEF_OUTPUT='HQ.png'

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    basFilename = args[0]
    bf = H5BasFile.BasFile (basFilename)

    numZ          = bf.numZMWs()

    framesPerBin  = opt.bin * H5BasFile.frameRate

    preHQ  = []
    inHQ   = []
    postHQ = []

    numEmpty = 0

    for hole in bf.holeNumbers():

        HQStart, HQEnd = bf.HQregion(hole)[2:4]

        if HQEnd == 0:
            numEmpty += 1
        else:

            preHQTime  = bf.elapsedFrames(hole, 0,   HQStart)
            inHQTime   = bf.elapsedFrames(hole, HQStart, HQEnd)
            postHQTime = bf.elapsedFrames(hole, HQEnd, bf.readLen(hole))

            HQStartBin = int(preHQTime  / framesPerBin)
            HQEndBin   = int(inHQTime   / framesPerBin) + HQStartBin
            postBin    = int(postHQTime / framesPerBin) + HQEndBin

            incrementRange (preHQ,  0, HQStartBin)
            incrementRange (inHQ,   HQStartBin, HQEndBin)
            incrementRange (postHQ, HQEndBin, postBin)

    logger.debug("%d of %d ZMWs had no HQ region" % (numEmpty, numZ))

    # Make sure the bin arrays are all the same size.

    maxLen = max(len(preHQ), len(inHQ), len(postHQ))
    if len(preHQ) < maxLen:
        preHQ.extend( [0] * (maxLen-len(preHQ)))
    if len(inHQ) < maxLen:
        inHQ.extend(  [0] * (maxLen-len(inHQ)))
    if len(postHQ) < maxLen:
        postHQ.extend([0] * (maxLen-len(postHQ)))

    # Stack the plots. Although it destroys the chronological
    # progression of states, we'll put the count of inHQ ZMWs on the
    # bottom, so it's zero-based, because it's the most interesting of
    # the three numbers.

    for ix in xrange(maxLen):
        preHQ[ix]  += inHQ[ix]
        postHQ[ix] += preHQ[ix]

    xCoords = xrange(0, len(postHQ)*opt.bin, opt.bin)

    plt.plot(xCoords, preHQ,  label='pre HQ',  color='yellow')
    plt.plot(xCoords, inHQ,   label='in HQ',   color='green')
    plt.plot(xCoords, postHQ, label='post HQ', color='red')

    plt.fill_between(xCoords, preHQ,  inHQ,  facecolor='yellow')
    plt.fill_between(xCoords, inHQ,   0,     facecolor='green')
    plt.fill_between(xCoords, postHQ, preHQ, facecolor='red')

    plt.title('ZMWs in High-Quality region vs time')
    plt.legend(loc='best', prop={'size':10})
    plt.savefig (opt.output)

    logger.debug("complete")

def incrementRange (list, start, end):

    listLen = len(list)
    if listLen < end:
        list.extend([0] * (end-listLen))

    for ix in xrange(start, end):
        list[ix] += 1

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file>')

    parser.add_option ('--bin',   type='int', help='plot bin size in seconds (def: %default)')
    parser.add_option ('--output',            help='Output file name (def: %default)')

    parser.set_defaults (bin=DEF_BINSIZE,
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
