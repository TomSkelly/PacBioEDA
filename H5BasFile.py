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

        self._consBasecalls = self._top["PulseData/ConsensusBaseCalls"]
        self._consZMW       = self._top["PulseData/ConsensusBaseCalls/ZMW"]
        self._consPasses    = self._top["PulseData/ConsensusBaseCalls/Passes"]

        self._PreBaseFrames = self._basecalls["PreBaseFrames"]      
        self._WidthInFrames = self._basecalls["WidthInFrames"]      

        self._maxZMW = max(self._ZMW["HoleNumber"])     # this takes a surprisingly long time to compute
                                      
        self._ZMWIndex      = None
        self._fillZMWIndex()          # why now? See note in _fillZMWIndex

        self._consPassIndex = None
        self._fillPassIndex()         # same here

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

    def holeStatusStr (self, hole):               # returns a string
        index = self._ZMWIndex[hole]
        return holeStatusTable[self._holeStatus[index]]

    def isSequencingZMW (self, hole):
        index = self._ZMWIndex[hole]
        return self._holeStatus[index] == 0

    def productivity (self, hole):
        index = self._ZMWIndex[hole]
        return self._productivity[index]

    def readLen (self, hole):
        index = self._ZMWIndex[hole]
        return self._ZMW["NumEvent"][index]

    def holeXY (self, hole):
        index = self._ZMWIndex[hole]
        holeXY = self._ZMW["HoleXY"]
        return holeXY[index,:]

    def numZMWs (self):        # N = number of ZMWs reported in file, not necessarily numbered 0..N-1
        return len(self._ZMW["HoleNumber"])

    def maxZMW (self):         # maximum ZMW# (not max+1, which would perhaps be more pythonic...)
        return self._maxZMW

    def holeNumbers (self):
        '''Return HDF5 array of hole numbers which actually exist in /PulseData/BaseCalls/ZMW/HoleNumber.'''
        return self._ZMW["HoleNumber"]

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
            end = self.readLen(hole)

        frames = 0
        index      = self._getBasecallIndex()[hole]
        startIndex = index + start
        endIndex   = index + end

        if end > start:
            frames = sum(self._PreBaseFrames[startIndex:endIndex]) \
                   + sum(self._WidthInFrames[startIndex:endIndex])

        return frames

    def _fillZMWIndex (self):
        '''Create and cache a table by hole number of starting indexes into PulseData/BaseCalls/ZMW.'''

        # ZMWIndex is always filled when the BasFile object is
        # created, unlike the case of basecallIndex and regionIndex,
        # which are filled only when they are first accessed. ZMWIndex
        # will be needed for any meaningful bas file access. Creating
        # it initially eliminates the need to check, on every access,
        # whether it exists.

        numHoles = self.numZMWs()
        maxHole  = self.maxZMW()

        logger.debug("creating ZMW index for %d holes, last hole is %d" % (numHoles, maxHole))
            
        self._ZMWIndex = [None] * (maxHole+1)

        holeNumber = self._ZMW["HoleNumber"]
        numEvent   = self._ZMW["NumEvent"]

        for ix in xrange(numHoles):
            hole = holeNumber[ix]
            self._ZMWIndex[hole] = ix

        logger.debug("complete")

    def _getBasecallIndex (self):
        '''Create and cache a table by hole number of starting indexes into PulseData/BaseCalls.'''

        if self._basecallIndex == None:

            numHoles = self.numZMWs()
            maxHole  = self.maxZMW()

            logger.debug("creating basecall index for %d holes, last hole is %d" % (numHoles, maxHole))
            
            self._basecallIndex = [None] * (maxHole+1)

            holeNumber = self._ZMW["HoleNumber"]
            numEvent   = self._ZMW["NumEvent"]
            index      = 0                        # index into basecalls array

            for ix in xrange(numHoles):
                hole = holeNumber[ix]
                self._basecallIndex[hole] = index
                index += numEvent[ix]

            logger.debug("processed %d basecalls" % index)

        return self._basecallIndex

    def _fillRegionTables (self):
        '''Create and cache tables by hole number of starting and HQ indexes into PulseData/Regions.'''

        # Both the start index and the HQ index require walking the
        # regions table. Might as well do both in one pass.

        if self._regionIndex == None:

            numHoles = self.numZMWs()
            maxHole  = self.maxZMW()

            logger.debug("creating region index for %d holes, last hole is %d" % (numHoles, maxHole))
            
            self._regionIndex = [None] * (maxHole+1)
            self._HQIndex     = [None] * (maxHole+1)

            regions  = self._regions
            index    = 0
            lastHole = -1        # init to non-matching value

            for line in regions:

                hole, regionType = line[0:2]

                if hole < 0 or hole > maxHole:          # sanity check hole# from region table
                    raise RuntimeError("hole number %d out of range" % (hole))

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

        # Note that, per python convention, range returned is (start..end-1).

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

        # Note that, per python convention, range returned is (start..end-1).

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

    def consReadLen (self, hole):
        index = self._ZMWIndex[hole]
        return self._consZMW["NumEvent"][index]

    def _fillPassIndex (self):
        '''Create and cache a table by hole number of starting indexes into PulseData/ConsensusBaseCalls/Passes.'''

        # Same story here as for fillZMWIndex: It's quick, call it
        # from the constructor, then there's no need to check whether
        # it exists every time we access it.

        # Note that we'll create a consPassIndex entry for every hole
        # which exists, whether or not it has any consensus
        # passes. That means we rely on the user to check that
        # numPasses is non-zero before using consPassIndex.

        numHoles = self.numZMWs()
        maxHole  = self.maxZMW()

        logger.debug("creating consensus pass index for %d holes, last hole is %d" % (numHoles, maxHole))
            
        self._consPassIndex = [None] * (maxHole+1)

        holeNumber = self._consZMW["HoleNumber"]
        passes     = self._consPasses
        numPasses  = passes["NumPasses"]
        index      = 0            # index into arrays in passes group

        for ix in xrange(numHoles):
            hole = holeNumber[ix]
            self._consPassIndex[hole] = index
            index += numPasses[ix]

        logger.debug("processed %d passes" % index)

    def _getConsensusBasecallIndex (self):
        '''Create and cache a table by hole number of starting indexes into PulseData/ConsensusBaseCalls.'''

        # Dual of _getBasecallIndex.

        if self._consensusIndex == None:

            self.ZMWSanityClause()     # make sure PulseData/BaseCalls matches PulseData/ConsensusBaseCalls

            logger.debug("creating consensus index")
            
            numHoles = self.numZMWs()
            maxHole  = self.maxZMW()

            self._consensusIndex = [None] * (maxHole+1)

            numEvent   = self._consZMW["NumEvent"]
            holeNumber = self._consZMW["HoleNumber"]
            index      = 0                            # index into basecalls array

            for ix in xrange(numHoles):               # for all ZMWs
                hole = holeNumber[ix]
                self._consensusIndex[hole] = index
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

        if not self._sanityChecked:               # we only need to check once

            logger.debug("performing sanity check.")
            
            holeNumber     = self._ZMW["HoleNumber"]
            holeStatus     = self._ZMW["HoleStatus"]
            holeXY         = self._ZMW["HoleXY"]

            holeNumberCons = self._consZMW["HoleNumber"]
            holeStatusCons = self._consZMW["HoleStatus"]
            holeXYCons     = self._consZMW["HoleXY"]

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
        index = self._ZMWIndex[hole]

        for passNo in xrange(numPasses[index]):
            ix = passIndex[hole] + passNo
            yield {"AdapterHitAfter":  passes["AdapterHitAfter"][ix], \
                   "AdapterHitBefore": passes["AdapterHitBefore"][ix], \
                   "PassDirection":    passes["PassDirection"][ix], \
                   "PassNumBases":     passes["PassNumBases"][ix], \
                   "PassStartBase":    passes["PassStartBase"][ix]}

        return

    def getConsensusBasecallField (self, field, hole, start=0, end=None):
        '''Return one field of the basecalls dataset from a region of a ZMW as an array.'''

        if end == None:
            end = self.consReadLen(hole)

        baseData = self._consBasecalls[field]
        index = self._getConsensusBasecallIndex()[hole]

        return baseData[index+start:index+end]
        
    def getConsensusSequence (self, hole, start=0, end=None):
        '''Return the basecalls from a region of a ZMW as an ascii string.'''

        if end == None:
            end = self.consReadLen(hole)

        baseData = self._consBasecalls["Basecall"]
        index = self._getConsensusBasecallIndex()[hole]
        intBases = [chr(baseData[x]) for x in xrange(index+start, index+end)]

        return ''.join(intBases)
        
    def getRevCompConsensusSequence (self, hole, start=0, end=None):
        '''Return the REVERSE_COMPLEMENTED basecalls from a region of a ZMW as an ascii string.'''

        if end == None:
            end = self.consReadLen(hole)

        baseData = self._consBasecalls["Basecall"]
        index = self._getConsensusBasecallIndex()[hole]

        # Example: xrange(10,15) is [10,11,12,13,14]. xrange(14,9,-1) is [14,13,12,11,10]

        intBases = [COMPLEMENTS[baseData[x]] for x in xrange(index+end-1, index+start-1, -1)]

        return ''.join(intBases)


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
