#!/usr/bin/env python

# Support for common operations reading a PacBio bas.h5 file.

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

import sys
import re
import h5py
from tt_log import logger

# Useful constants
holeStatusTable = ('SEQ', 'ANTIH', 'FIDUC', 'SUSP', 'ANTIM', 'FDZMW', 'FBZMW', 'ANTIB', 'OUT')
regionTypeTable   = ('AD', 'IN', 'HQ')
frameRate = 75.0001831055      # frames per second

ADAPTER = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"       # adapter sequence

COMPLEMENTS = [None] * 256
COMPLEMENTS[ord('A')] = 'T'
COMPLEMENTS[ord('C')] = 'G'
COMPLEMENTS[ord('G')] = 'C'
COMPLEMENTS[ord('T')] = 'A'
COMPLEMENTS[ord('N')] = 'N'

class BasFile (object):

    def __init__ (self, filename):

        logger.debug("creating BasFile object")

        self._filename   = filename
        self._infile     = h5py.File (filename, 'r')
        self._top        = h5py.Group (self._infile, '/')

        self._basecalls  = self._top["PulseData/BaseCalls"]
        self._ZMW        = self._top["PulseData/BaseCalls/ZMW"]
        self._ZMWMetrics = self._top["PulseData/BaseCalls/ZMWMetrics"]
        self._regions    = self._top["PulseData/Regions"]
        self._numRegions = self._regions.shape[0]

        self._holeStatus   = self._ZMW["HoleStatus"]
        self._productivity = self._ZMWMetrics["Productivity"]


####        self._Basecall        = self._basecalls["Basecall"]           
####        self._DeletionQV      = self._basecalls["DeletionQV"]         
####        self._DeletionTag     = self._basecalls["DeletionTag"]        
####        self._InsertionQV     = self._basecalls["InsertionQV"]        
        self._PreBaseFrames   = self._basecalls["PreBaseFrames"]      
####        self._SubstitutionQV  = self._basecalls["SubstitutionQV"]     
####        self._SubstitutionTag = self._basecalls["SubstitutionTag"]    
####        self._QualityValue    = self._basecalls["QualityValue"]       
        self._WidthInFrames   = self._basecalls["WidthInFrames"]      
                                      
        self._basecallIndex = None
        self._regionIndex   = None
        self._HQIndex       = None
        self._coords        = None

        match = re.search('_s(\d+)_p(\d+)$', self.movieName())
        self._setNumber    = int(match.group(1))        # match.group(0) is entire match
        self._strobeNumber = int(match.group(2))
        
    def __del__ (self):
        self._infile.close()

    def movieName (self):
        return self._top["ScanData/RunInfo"].attrs["MovieName"]

    def setNumber (self):
        return self._setNumber

    def strobeNumber (self):
        return self._strobeNumber

    def basecalls (self):
        return self._basecalls

    def ZMW (self):
        return self._ZMW

    def ZMWMetrics (self):
        return self._ZMWMetrics

    def regions (self):
        return self._regions

    def numRegions (self):
        return self._numRegions

    def numZMWs (self):
        return len(self._ZMW["HoleNumber"])

    def holeStatusStr (self, hole):               # returns a string
        return holeStatusTable[self._holeStatus[hole]]

    def isSequencingZMW (self, hole):
        return self._holeStatus[hole] == 0

    def productivity (self, hole):
        return self._productivity[hole]

    def readLen (self, hole):
        return self._ZMW["NumEvent"][hole]

    def holeXY (self, hole):
        holeXY = self._ZMW["HoleXY"]
        return holeXY[hole,:]

    def basecallIndex (self):
        '''Create and cache a table by hole number of starting indexes into PulseData/BaseCalls.'''

        if self._basecallIndex == None:

            logger.debug("creating basecall index")
            
            numEvent   = self._ZMW["NumEvent"]
            holeNumber = self._ZMW["HoleNumber"]
            numZ       = self.numZMWs()
            index      = 0
            self._basecallIndex = [0] * numZ

            # The loop below includes a check that the HoleNumbers run
            # from 0 to N-1. I.e., HoleNumber[ix] == ix. Otherwise,
            # there is no way to get from the hole number of a region
            # back to the ZMW entry (other than searching for it).

            for ix in xrange(numZ):
                
                if holeNumber[ix] != ix:
                    raise RuntimeError("Hole number != index at %d" % ix)

                self._basecallIndex[ix] = index
                index += numEvent[ix]

            logger.debug("processed %d basecalls" % index)

        return self._basecallIndex

    def cellCoords (self):
        '''Find and cache the minimum and maximum X/Y coordinates on the SMRTcell'''

        if self._coords is None:

            logger.debug("finding SMRTcell coordinates")

            numZ       = self.numZMWs()
            holeXY = self._ZMW["HoleXY"]

            minX = maxX = holeXY[0,0]
            minY = maxY = holeXY[0,1]

            for ix in xrange(numZ):

                x,y = holeXY[ix,:]

                if x < minX:
                    minX = x
                elif x > maxX:
                    maxX = x

                if y < minY:
                    minY = y
                elif y > maxY:
                    maxY = y

            self._coords = (minX, maxX, minY, maxY)

            logger.debug("SMRTcell is (%d,%d,%d,%d)" % (minX, maxX, minY, maxY))

        return self._coords

    def _fillRegionTables (self):
        '''Create and cache tables by hole number of starting and HQ indexes into PulseData/Regions.'''

        if self._regionIndex == None:

            logger.debug("creating region index")
            
            self._regionIndex = [0] * self.numZMWs()
            self._HQIndex     = [0] * self.numZMWs()

            regions  = self._regions
            index    = 0
            lastHole = -1        # init to non-matching value

            for line in regions:

                hole, regionType = line[0:2]

                if hole != lastHole:                    # start of new hole?
                    self._regionIndex[hole] = index
                    lastHole = hole

                if regionType == 2:                     # HQ region for this hole?
                    self._HQIndex[hole] = index

                index += 1
                    
            logger.debug("processed %d regions" % index)

        return self._regionIndex

    def regionIndex (self):
        '''Return lookup table by hole number of starting index into PulseData/Regions.'''

        if self._regionIndex == None:
            self._fillRegionTables()          # creates self._regionIndex and self._HQIndex

        return self._regionIndex              # this is a list, indexed by hole

    def HQregion (self, hole):
        '''Return the HQ region tuple from the regions dataset for a specified hole.'''

        if self._HQIndex == None:
            self._fillRegionTables()          # creates self._regionIndex and self._HQIndex

        return self._regions[self._HQIndex[hole]]   # this is a region tuple: hole, type, start, end, score

    def holeRegions (self, hole):
        '''Generator which returns regions for specified hole, one by one'''

        regions     = self._regions           # the regions dataset
        regionIndex = self.regionIndex()      # table of first region indexes by hole#
        nextRegionIndex = regionIndex[hole]   # first region for this hole

        while nextRegionIndex < self._numRegions and regions[nextRegionIndex][0] == hole:
            yield regions[nextRegionIndex]
            nextRegionIndex += 1

        return

    def getBasecallField (self, field, hole, start=0, end=None):
        '''Return one field of the basecalls dataset from a region of a ZMW as an array.'''

        if end == None:
            end = self.readLen(hole)

        baseData = self._basecalls[field]
        index = self.basecallIndex()[hole]

        return baseData[index+start:index+end]
        
    def getSequence (self, hole, start=0, end=None):
        '''Return the basecalls from a region of a ZMW as an ascii string.'''

        # The logic of getBasecallField is duplicated here, rather
        # than called, for efficiency, since this routine is expected
        # to be called a lot.

        if end == None:
            end = self.readLen(hole)

        baseData = self._basecalls["Basecall"]
        index = self.basecallIndex()[hole]
        intBases = [chr(baseData[x]) for x in xrange(index+start, index+end)]

        return ''.join(intBases)
        
    def getRevCompSequence (self, hole, start=0, end=None):
        '''Return the REVERSE_COMPLEMENTED basecalls from a region of a ZMW as an ascii string.'''

        if end == None:
            end = self.readLen(hole)

        baseData = self._basecalls["Basecall"]
        index = self.basecallIndex()[hole]

        # Example: xrange(10,15) is [10,11,12,13,14]. xrange(14,9,-1) is [14,13,12,11,10]

        intBases = [COMPLEMENTS[baseData[x]] for x in xrange(index+end-1, index+start-1, -1)]

        return ''.join(intBases)
        
    def elapsedFrames (self, hole, start=0, end=None):
        '''Compute the number of frames in a sub-region of a read.'''

        frames = 0

        if end == None:
            end = self._ZMW["NumEvent"][hole]

        index      = self.basecallIndex()[hole]
        startIndex = index + start
        endIndex   = index + end

        if end > start:
            frames = sum(self._PreBaseFrames[startIndex:endIndex]) \
                   + sum(self._WidthInFrames[startIndex:endIndex])

        return frames

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
