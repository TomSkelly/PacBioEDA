#!/usr/bin/env python

# Support for common operations reading a PacBio bas.h5 file.

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

import sys
import numpy as np
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

        self._consZMW    = self._top["PulseData/ConsensusBaseCalls/ZMW"]
        self._consPasses = self._top["PulseData/ConsensusBaseCalls/Passes"]

        self._consPassIndex = list(np.cumsum(self._consPasses["NumPasses"]))
        self._consPassIndex.insert(0,0)      # ZMW 0 starts at 0

####        self._Basecall        = self._basecalls["Basecall"]           
####        self._DeletionQV      = self._basecalls["DeletionQV"]         
####        self._DeletionTag     = self._basecalls["DeletionTag"]        
####        self._InsertionQV     = self._basecalls["InsertionQV"]        
        self._PreBaseFrames   = self._basecalls["PreBaseFrames"]      
####        self._SubstitutionQV  = self._basecalls["SubstitutionQV"]     
####        self._SubstitutionTag = self._basecalls["SubstitutionTag"]    
####        self._QualityValue    = self._basecalls["QualityValue"]       
        self._WidthInFrames   = self._basecalls["WidthInFrames"]      
                                      
        self._basecallIndex  = None
        self._consensusIndex = None
        self._regionIndex    = None
        self._HQIndex        = None
        self._coords         = None

        self._sanityChecked = False

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

    def cellCoords (self):
        '''Find and cache the minimum and maximum X/Y coordinates on the SMRTcell'''

        if self._coords is None:

            logger.debug("finding SMRTcell coordinates")

            holeXY = self._ZMW["HoleXY"]

            minX = min(holeXY[:,0])
            maxX = max(holeXY[:,0])
            minY = min(holeXY[:,1])
            maxY = max(holeXY[:,1])

            self._coords = (minX, maxX, minY, maxY)

            logger.debug("SMRTcell is (%d,%d,%d,%d)" % (minX, maxX, minY, maxY))

        return self._coords

    def elapsedFrames (self, hole, start=0, end=None):
        '''Compute the number of frames in a sub-region of a read.'''

        if end == None:
            end = self._ZMW["NumEvent"][hole]

        frames = 0
        index      = self._getBasecallIndex()[hole]
        startIndex = index + start
        endIndex   = index + end

        if end > start:
            frames = sum(self._PreBaseFrames[startIndex:endIndex]) \
                   + sum(self._WidthInFrames[startIndex:endIndex])

        return frames

    def _getBasecallIndex (self):
        '''Create and cache a table by hole number of starting indexes into PulseData/BaseCalls.'''

        if self._basecallIndex == None:

            logger.debug("creating basecall index")
            
            numEvent   = self._ZMW["NumEvent"]
            numZ       = self.numZMWs()
            index      = 0                        # index into basecalls array
            self._basecallIndex = [0] * numZ

            for ix in xrange(numZ):               # for all ZMWs
                self._basecallIndex[ix] = index
                index += numEvent[ix]

            logger.debug("processed %d basecalls" % index)

        return self._basecallIndex

    def _fillRegionTables (self):
        '''Create and cache tables by hole number of starting and HQ indexes into PulseData/Regions.'''

        # Both the start index and the HQ index require walking the
        # regions table. Might as well do both in one pass.

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
        index = self._getBasecallIndex()[hole]

        return baseData[index+start:index+end]
        
    def getSequence (self, hole, start=0, end=None):
        '''Return the basecalls from a region of a ZMW as an ascii string.'''

        # The logic of getBasecallField is duplicated here, rather
        # than called, for efficiency, since this routine is expected
        # to be called a lot.

        if end == None:
            end = self.readLen(hole)

        baseData = self._basecalls["Basecall"]
        index = self._getBasecallIndex()[hole]
        intBases = [chr(baseData[x]) for x in xrange(index+start, index+end)]

        return ''.join(intBases)
        
    def getRevCompSequence (self, hole, start=0, end=None):
        '''Return the REVERSE_COMPLEMENTED basecalls from a region of a ZMW as an ascii string.'''

        if end == None:
            end = self.readLen(hole)

        baseData = self._basecalls["Basecall"]
        index = self._getBasecallIndex()[hole]

        # Example: xrange(10,15) is [10,11,12,13,14]. xrange(14,9,-1) is [14,13,12,11,10]

        intBases = [COMPLEMENTS[baseData[x]] for x in xrange(index+end-1, index+start-1, -1)]

        return ''.join(intBases)

        
    # The methods which follow are for accessing the consensus
    # basecalls. Most of these routines have logic which is very
    # similar to that of the routines for accessing the raw reads. But
    # the structure, naming conventions, and semantics of the
    # consensus data are different enough from the raw case that
    # combining the access logic would have lead to a lot of
    # complication. So I've opted to duplicate the logic.

    # The structure of the consensus-related data in a bas.h5 file looks like:

    #   PulseData                (Group)
    #     ConsensusBaseCalls     (Group)
    #       Basecall             [0..K]
    #       DeletionQV           [0..K]
    #       DeletionTag          [0..K]
    #       InsertionQV          [0..K]
    #       QualityValue         [0..K]
    #       SubstitutionQV       [0..K]
    #       SubstitutionTag      [0..K]
    #       Passes               (Group)
    #         NumPasses          [0..N]   <-- Note
    #         AdapterHitAfter    [0..P]
    #         AdapterHitBefore   [0..P]
    #         PassDirection      [0..P]
    #         PassNumBases       [0..P]
    #         PassStartBase      [0..P]
    #       ZMW                  (Group)
    #         HoleNumber         [0..N]
    #         HoleStatus         [0..N]
    #         HoleXY             [0..N,0..1]
    #         NumEvent           [0..N]

    # where N is the highest ZMW number, K is the total number of
    # consensus bases, and P is the total number of passes.  Given the
    # way the datasets are indexed, the nesting is a bit awkward
    # (IMHO).

    def _getConsensusBasecallIndex (self):
        '''Create and cache a table by hole number of starting indexes into PulseData/ConsensusBaseCalls.'''

        # Dual of _getBasecallIndex.

        if self._consensusIndex == None:

            self.ZMWSanityClause()

            logger.debug("creating consensus index")
            
            numEvent   = self._consZMW["NumEvent"]
            holeNumber = self._consZMW["HoleNumber"]
            numZ       = self.numZMWs()
            index      = 0                        # index into basecalls array
            self._consensusIndex = [0] * numZ

            for ix in xrange(numZ):               # for all ZMWs
                self._consensusIndex[ix] = index
                index += numEvent[ix]

            logger.debug("processed %d consensus basecalls" % index)

        return self._consensusIndex

    def ZMWSanityClause (self):
        '''The trouble with a lot of code today is that there ain't no sanity clause.'''

        # Check that the subfields of PulseData/BaseCalls/ZMW and
        # PulseData/ConsensusBaseCalls/ZMW match, except for
        # NumEvent. This is only a sanity check, and (assuming it
        # passes) has no effect on the class's function. I.e., it
        # doesn't need to be called, if you're the trusting sort.

        # There is also a check that HoleNumber[ix] == ix. Otherwise,
        # there is no way to get from the hole number of a region back
        # to the ZMW entry (other than searching for it).

        if not self._sanityChecked:               # we only need to check once

            logger.debug("performing sanity check.")
            
            holeNumber     = self._ZMW["HoleNumber"]
            holeStatus     = self._ZMW["HoleStatus"]
            holeXY         = self._ZMW["HoleXY"]

            holeNumberCons = self._consZMW["HoleNumber"]
            holeStatusCons = self._consZMW["HoleStatus"]
            holeXYCons     = self._consZMW["HoleXY"]

            numZ = self.numZMWs()

            if np.any(holeNumber[:] != xrange(numZ)):
                raise RuntimeError("Hole number != index")

            if np.any(holeNumber[:] != holeNumberCons[:]):
                raise RuntimeError("Hole number != consensus hole number")

            if np.any(holeStatus[:] != holeStatusCons[:]):
                raise RuntimeError("Hole status != consensus hole status")

            if np.any(holeXY[:] != holeXYCons[:]):
                raise RuntimeError("Hole XY != consensus hole XY")

            self._sanityChecked = True
            logger.debug("sanity check passed")

        return                      # nothing returned -- if we return at all, we passed!

    def holeConsensusPasses (self, hole):
        '''Generator which returns consensus passes for specified hole, one by one'''

        # Dual of holeRegions, sort of.

        passIndex = self._consPassIndex
        passes    = self._consPasses
        numPasses = passes["NumPasses"]

        for passNo in xrange(numPasses[hole]):
            ix = passIndex[hole] + passNo
            yield {"AdapterHitAfter":  passes["AdapterHitAfter"][ix], \
                   "AdapterHitBefore": passes["AdapterHitBefore"][ix], \
                   "PassDirection":    passes["PassDirection"][ix], \
                   "PassNumBases":     passes["PassNumBases"][ix], \
                   "PassStartBase":    passes["PassStartBase"][ix]}

        return

    # TODO: Need to implement getConsensusBasecallField, getConsensusSequence, and getRevCompConsensusSequence.


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
