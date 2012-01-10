#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Create a plot of ZMW productivity by x/y position on the
# SMRTcell. First parameter is input .bas.h5 file. Output png file is
# optional command line parameter, defaulting to productivity.png.

import sys
import optparse
import numpy as np
import h5py
from tt_log import logger

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEF_OUTPUT = 'productivity.png'

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    infile_name = args[0]
    infile = h5py.File (infile_name, 'r')

    colours = ('grey', 'red', 'green')
    legends = ('non-seq', 'prod-0', 'prod-1')

    top = h5py.Group (infile, '/')

    ZMW        = top["PulseData/BaseCalls/ZMW"]
    ZMWMetrics = top["PulseData/BaseCalls/ZMWMetrics"]
    
    holeStatus = ZMW["HoleStatus"]
    holeXY     = ZMW["HoleXY"]
    holeProd   = ZMWMetrics["Productivity"]

    nonseqHoles = holeStatus[:]!=0      # ZMWs other than sequencing
    prod0Holes  = np.logical_and(holeProd[:]==0, np.logical_not(nonseqHoles))
    prod1Holes  = np.logical_and(holeProd[:]==1, np.logical_not(nonseqHoles))

    holesByType = (nonseqHoles, prod0Holes, prod1Holes)

    for which in xrange(len(holesByType)):

        whichHoles = holesByType[which]
        howMany = sum(whichHoles)
        logger.debug("%5d  %s" % (howMany, legends[which]));
        if howMany > 0:
            plt.scatter (holeXY[whichHoles,0], holeXY[whichHoles,1], \
                             s=1, c=colours[which], edgecolor='face', \
                             label="%5d  %s" % (howMany, legends[which]))

    plt.axis    ('equal')
    plt.legend  (scatterpoints=3, prop={'size':8})
    plt.savefig (opt.output)

    infile.close()

    logger.debug("complete")

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file>')

    parser.add_option ('--output', help='Output file name (def: %default)')

    parser.set_defaults (output=DEF_OUTPUT)

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
