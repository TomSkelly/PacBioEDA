#!/usr/bin/env python

# Support for common operations reading a PacBio bas.h5 file.

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Historical note:

# In version 2.0.0 of the PacBio software, the concept of a bas.h5
# file changed. Pre-2.0.0, a bas.h5 file contained THE STUFF. Now THE
# STUFF is divided among three bax.h5 files, each containing N/3 of
# the ZMWs. The bas.h5 file is just a thin level of indirection on top
# of that, containing only the mapping of ZMW# to bax.h5 file.

# How should we accommodate that change? Hopefully in a manner which
# is both backward compatible and efficient, without being too ugly.

# One possibility would be to turn the existing BasFile class into a
# Baxfile class, and to create a new BasFile class whose job is to
# manage several (3) BasFiles. Each method of the new BasFile would
# just decide which BaxFile object to access, and invoke the
# corresponding method on that object. That would be conceptually
# straightforward, but ugly (due to the duplication of methods) and
# slow (requiring a double method invocation where only one was needed
# before). The duplication could probably be dealt with by using a
# decorator -- but the approach would still be slow.

# So I've gone with Plan B: Define a new BaxFile object to contain the
# data that used to be kept as attributes of the BasFile. BasFile
# methods, rather than accessing "self.whatever", must now access
# bax.whatever, having first decided which of several bax's contains
# the data of interest.

import os
import sys
import numpy as np
import re                       # for regular expressions
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

class BasFile (object):    # see header comments describing this class

    def __init__ (self, filename, CCSDir=None):

        logger.debug('creating BasFile object')

        self._filename  = filename
        self._CCSDir    = CCSDir
        self._infile    = h5py.File (filename, 'r')
        self._top       = self._infile
        self._baxfile   = list()
        self._coords    = None

        self._consensusIndex = None
        self._consPassIndex  = None

        if 'MultiPart' not in self._top:       # if this is an old-style bas file
            bf = BaxFile(filename)             # file will contain its own CCSdata
            self._baxfile.append(bf)           # only one file (this one) in the list

        else:                                  # else it's an index to a set of bax files

            h5Dir = os.path.dirname (os.path.abspath (self._filename))

            for baxfileName in self._top['MultiPart/Parts']:      # for each bax file

                fqBaxfileName = os.path.join(h5Dir, baxfileName)  # fq = fully qualified
                bf = BaxFile(fqBaxfileName, CCSDir=CCSDir)
                self._baxfile.append(bf)       # add file to list

        self.fillCombinedFields()              # need to compute this first, we'll need it later
        self.fillZMWIndexes()
        self.fillMovieName()
        self.fillRegionIndexes()

    def fillCombinedFields (self):
        '''Called from __init__ to compute aggregated fields across all ax files.'''

        self._maxZMW  = 0
        self._numZMWs = 0
        self._hasConsensus = True

        for bf in self._baxfile:
            self._maxZMW = max (self._maxZMW, bf._maxZMW)     # largest ZMW#
            self._numZMWs += bf.numZMWs()                     # count of ZMWs
            if not bf.hasConsensus():
                self._hasConsensus = False

        logger.debug("largest ZMW# is %d" % self._maxZMW)

    def fillZMWIndexes (self):
        '''Called from __init__ to compute _baxByHole, _ZMWindex and _basecallIndex.'''

        # _baxByHole points to the bax object for the bax file which
        # contains a given ZMW.

        # _ZMWIndex is the offset for a given ZMW into the datasets of
        # the PulseData/BaseCalls/ZMW group in the bax file pointed to
        # by _baxByHole. 

        # Likewise, _basecallIndex is the offset to the first entry
        # for a given ZMW in the datasets of the PulseData/BaseCalls
        # group. Each of those datasets contains numEvent entries for
        # a given ZMW.

        # The alert reader will probably notice that we never examine
        # the MultiPart/HoleLookup dataset in the bas.h5
        # file. Instead, we derive the ZMW#-to-baxfile mapping by
        # looking at the ZMW#s contained in each bax file, in the loop
        # below.

        self._baxByHole     = [None] * (self._maxZMW+1)
        self._ZMWIndex      = [None] * (self._maxZMW+1)
        self._basecallIndex = [None] * (self._maxZMW+1)

        for bf in self._baxfile:

            ix = 0
            eventIndex = 0
            numEvent   = bf._ZMW["NumEvent"]

            for hole in bf.holeNumbers():

                self._baxByHole[hole] = bf
                self._ZMWIndex[hole] = ix
                self._basecallIndex[hole] = eventIndex
                eventIndex += numEvent[ix]
                ix += 1

            logger.debug("%s processed %d ZMWs"      % (bf.shortName(), ix))
            logger.debug("%s processed %d basecalls" % (bf.shortName(), eventIndex))

    def fillRegionIndexes (self):
        '''Called from __init__ to compute _regionIndex and _HQIndex.'''

        # _regionIndex is the offset to the first entry for a given
        # ZMW in the PulseData/Regions dataset in the bas file where
        # the hole resides.

        # _HQIndex is the offset to the HQ-region entry in
        # PulseData/Regions for a given ZMW. (There is only one
        # ... for now.)

        self._regionIndex = [None] * (self._maxZMW+1)
        self._HQIndex     = [None] * (self._maxZMW+1)

        for bf in self._baxfile:

            regions = bf._regions
            regionIndex = 0
            lastHole = -1                               # init to non-matching value

####            for line in regions:                             #### doing it this way was *very* slow
####                hole, regionType = line[0:2]

            for hole, regionType in regions[:,0:2]:

                if hole != lastHole:                    # start of new hole?

                    if hole < 0 or hole > self._maxZMW: # sanity check hole# from region table
                        raise RuntimeError("hole number %d out of range" % (hole))

                    self._regionIndex[hole] = regionIndex
                    lastHole = hole

                if regionType == 2:                     # HQ region for this hole?
                    self._HQIndex[hole] = regionIndex

                regionIndex += 1
                    
            logger.debug("%s processed %d regions" % (bf.shortName(), regionIndex))

    def fillMovieName (self):
        '''Called from __init__ to compute _movieName, _setNumber and _strobeNumber.'''

        self._movieName = None

        for bf in self._baxfile:

            if self._movieName is None:
                self._movieName = bf._movieName
            elif self._movieName != bf._movieName:
                raise RuntimeError ('bax file movie names are different')

        match = re.search('_s(\d+)_[pX](\d+)$', self._movieName)
        self._setNumber    = int(match.group(1))        # match.group(0) is entire match
        self._strobeNumber = int(match.group(2))

    # ------------ End of initialization routines --------------------

    def __del__ (self):
        self._infile.close()

    def holeNumbers (self):
        '''Generator function to return the sequence of hole numbers which exist in the bax file collection.'''
        for hole in xrange(len(self._baxByHole)):
            if self._baxByHole[hole] is not None:
                yield hole

        return

    def movieName (self):
        return self._movieName

    def setNumber (self):
        return self._setNumber

    def strobeNumber (self):
        return self._strobeNumber

    def maxZMW (self):         # maximum ZMW# (not max+1, which would perhaps be more pythonic...)
        return self._maxZMW

    def numZMWs (self):        # number of ZMWs reported in file, not necessarily numbered 0..N-1
        return self._numZMWs

    def hasConsensus (self):
        return self._hasConsensus

    def holeStatusStr (self, hole):               # returns a string
        bf    = self._baxByHole[hole]
        index = self._ZMWIndex[hole]
        return holeStatusTable[bf._holeStatus[index]]

    def isSequencingZMW (self, hole):
        bf    = self._baxByHole[hole]
        index = self._ZMWIndex[hole]
        return bf._holeStatus[index] == 0

    def productivity (self, hole):
        bf    = self._baxByHole[hole]
        index = self._ZMWIndex[hole]
        return bf._productivity[index]

    def readLen (self, hole):
        bf    = self._baxByHole[hole]
        index = self._ZMWIndex[hole]
        return bf._ZMW["NumEvent"][index]

    def holeXY (self, hole):
        bf    = self._baxByHole[hole]
        index = self._ZMWIndex[hole]
        holeXY = bf._ZMW["HoleXY"]
        return holeXY[index,:]

    def HQregion (self, hole):
        '''Return the HQ region tuple from the regions dataset for a specified hole.'''
        bf = self._baxByHole[hole]
        return bf._regions[self._HQIndex[hole]]   # this is a region tuple: hole, type, start, end, score

    def holeRegions (self, hole):
        '''Generator which returns regions for specified hole, one by one'''

        bf = self._baxByHole[hole]
        regions = bf._regions           # the regions dataset
        nextRegionIndex = self._regionIndex[hole]   # first region for this hole

        while nextRegionIndex < bf._numRegions and regions[nextRegionIndex][0] == hole:
            yield regions[nextRegionIndex]
            nextRegionIndex += 1

        return

    def cellCoords (self):
        '''Find and cache the minimum and maximum X/Y coordinates on the SMRTcell'''

        if self._coords is None:

            logger.debug("finding SMRTcell coordinates")

            minX = 0
            maxX = 0
            minY = 0
            maxY = 0

            for bf in self._baxfile:

                holeXY = bf._ZMW["HoleXY"]

                minX = min(minX, min(holeXY[:,0]))
                maxX = max(maxX, max(holeXY[:,0]))
                minY = min(minY, min(holeXY[:,1]))
                maxY = max(maxY, max(holeXY[:,1]))

            self._coords = (minX, maxX, minY, maxY)

            logger.debug("SMRTcell is (%d,%d,%d,%d)" % (minX, maxX, minY, maxY))

        return self._coords

    def elapsedFrames (self, hole, start=0, end=None):
        '''Compute the number of frames in a sub-region of a read.'''

        if end == None:
            end = self.readLen(hole)

        frames     = 0
        bf         = self._baxByHole[hole]
        index      = self._basecallIndex[hole]
        startIndex = index + start
        endIndex   = index + end

        if end > start:
            frames = sum(bf._PreBaseFrames[startIndex:endIndex]) \
                   + sum(bf._WidthInFrames[startIndex:endIndex])

        return frames

    def getBasecallField (self, field, hole, start=0, end=None):
        '''Return one field of the basecalls dataset from a region of a ZMW as an array.'''

        # Note that, per python convention, range returned is (start..end-1).

        if end == None:
            end = self.readLen(hole)

        bf       = self._baxByHole[hole]
        index    = self._basecallIndex[hole]
        baseData = bf._basecalls[field]

        return baseData[index+start:index+end]
        
    def getSequence (self, hole, start=0, end=None):
        '''Return the basecalls from a region of a ZMW as an ascii string.'''

        # The logic of getBasecallField is duplicated here, rather
        # than called, for efficiency, since this routine is expected
        # to be called a lot.

        # Note that, per python convention, range returned is (start..end-1).

        if end == None:
            end = self.readLen(hole)

        bf       = self._baxByHole[hole]
        index    = self._basecallIndex[hole]
        baseData = bf._basecalls["Basecall"]
        intBases = [chr(baseData[x]) for x in xrange(index+start, index+end)]

        return ''.join(intBases)
        
    def getRevCompSequence (self, hole, start=0, end=None):
        '''Return the REVERSE_COMPLEMENTED basecalls from a region of a ZMW as an ascii string.'''

        if end == None:
            end = self.readLen(hole)

        bf       = self._baxByHole[hole]
        index    = self._basecallIndex[hole]
        baseData = bf._basecalls["Basecall"]

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
        '''Consensus read length for a given ZMW.'''

        bf    = self._baxByHole[hole]
        index = self._ZMWIndex[hole]
        return bf._consZMW["NumEvent"][index]

    def numConsensusPasses (self, hole):
        '''Number of consensus passes for a given ZMW.'''

        bf    = self._baxByHole[hole]
        index = self._ZMWIndex[hole]
        return bf._consPasses["NumPasses"][index]

    def getConsPassIndex (self):
        '''Return consPassIndex, initializing it if necessary.'''

        if self._consPassIndex is None:
            self.fillConsensusIndexes()

        return self._consPassIndex

    def getConsensusIndex (self):
        '''Return consensusIndex, initializing it if necessary.'''

        if self._consensusIndex is None:
            self.fillConsensusIndexes()

        return self._consensusIndex

    def fillConsensusIndexes (self):
        '''Compute _consPassIndex and _consensusIndex. Called only when those arrays are needed.'''

        self._consPassIndex  = [None] * (self._maxZMW+1)
        self._consensusIndex = [None] * (self._maxZMW+1)

        for bf in self._baxfile:

            bf.ZMWSanityClause()          # sanity check the consensus datasets

            ix = 0
            passIndex = 0
            consIndex = 0
            numPasses = bf._consPasses["NumPasses"]
            numEvent  = bf._consZMW["NumEvent"]

            for hole in bf.holeNumbers():
                self._consPassIndex[hole] = passIndex
                passIndex += numPasses[ix]
                self._consensusIndex[hole] = consIndex
                consIndex += numEvent[ix]
                ix += 1

            logger.debug("%s processed %d consensus passes"    % (bf.shortName(), passIndex))
            logger.debug("%s processed %d consensus basecalls" % (bf.shortName(), consIndex))

    def holeConsensusPasses (self, hole):
        '''Generator which returns consensus passes for specified hole, one by one'''

        # Dual of holeRegions, sort of.

        # If the logic here seems a bit contorted, it's because of the
        # strange way the consensus datasets are nested. See the table
        # above.

        bf        = self._baxByHole[hole]
        passes    = bf._consPasses
        numPasses = passes["NumPasses"]

        passIndex = self.getConsPassIndex()
        index     = self._ZMWIndex[hole]

        for passNo in xrange(numPasses[index]):
            ix = passIndex[hole] + passNo
            yield {"AdapterHitAfter":  passes["AdapterHitAfter"][ix], \
                   "AdapterHitBefore": passes["AdapterHitBefore"][ix], \
                   "PassDirection":    passes["PassDirection"][ix], \
                   "PassNumBases":     passes["PassNumBases"][ix], \
                   "PassStartBase":    passes["PassStartBase"][ix]}

        return

    def getConsensusBasecallField (self, field, hole, start=0, end=None):
        '''Return a subset of a specified dataset in hte ConsensusBaseCalls group as an array.'''

        if end == None:
            end = self.consReadLen(hole)

        bf       = self._baxByHole[hole]
        index    = self.getConsensusIndex()[hole]
        baseData = bf._consBasecalls[field]

        return baseData[index+start:index+end]
        
    def getConsensusSequence (self, hole, start=0, end=None):
        '''Return the consensus basecalls from a region of a ZMW as an ascii string.'''

        if end == None:
            end = self.consReadLen(hole)

        bf       = self._baxByHole[hole]
        index    = self.getConsensusIndex()[hole]
        baseData = bf._consBasecalls["Basecall"]
        intBases = [chr(baseData[x]) for x in xrange(index+start, index+end)]

        return ''.join(intBases)
        
    def getRevCompConsensusSequence (self, hole, start=0, end=None):
        '''Return the REVERSE_COMPLEMENTED basecalls from a region of a ZMW as an ascii string.'''

        if end == None:
            end = self.consReadLen(hole)

        bf       = self._baxByHole[hole]
        index    = self.getConsensusIndex()[hole]
        baseData = bf._consBasecalls["Basecall"]

        # Example: xrange(10,15) is [10,11,12,13,14]. xrange(14,9,-1) is [14,13,12,11,10]

        intBases = [COMPLEMENTS[baseData[x]] for x in xrange(index+end-1, index+start-1, -1)]

        return ''.join(intBases)

    def getConsensusAverageQuality (self, hole, start=0, end=None):
        '''Return the average (in Fred-space) of the quality scores for a consensus read.'''

        if end == None:
            end = self.consReadLen(hole)
        if end <= start:           # avoid zero-divide
            raise RuntimeError ('invalid range (%d-%d) specified' % (start, end))

        bf       = self._baxByHole[hole]
        index    = self.getConsensusIndex()[hole]
        baseData = bf._consBasecalls['QualityValue']

        return sum(baseData[index+start:index+end]) / (end-start)

class BaxFile (object):   # see header comments describing this class

    fileNum = 0           # class variable to count files

    def __init__ (self, filename, CCSDir=None):

        BaxFile.fileNum += 1
        self._shortName = 'bax-%d' % BaxFile.fileNum

        logger.debug("creating BaxFile object for %s = %s" % (self._shortName, filename))

        self._filename      = filename
        self._CCSDir        = CCSDir
        self._infile        = h5py.File (filename, 'r')
        self._top           = self._infile         # h5py 2.0.1 change!

        self._pulsedata     = self._top["PulseData"]
        self._basecalls     = self._top["PulseData/BaseCalls"]
        self._ZMW           = self._top["PulseData/BaseCalls/ZMW"]
        self._regions       = self._top["PulseData/Regions"]
        self._productivity  = self._top["PulseData/BaseCalls/ZMWMetrics/Productivity"]
        self._movieName     = self._top["ScanData/RunInfo"].attrs["MovieName"]
        self._holeStatus    = self._ZMW["HoleStatus"]
        self._numRegions    = self._regions.shape[0]

        self._PreBaseFrames = self._basecalls["PreBaseFrames"]      
        self._WidthInFrames = self._basecalls["WidthInFrames"]      

        self._maxZMW        = max(self._ZMW["HoleNumber"])     # this takes a surprisingly long time to compute
                                      
        self._sanityChecked = False

        self.findCCSFile()

    def __del__ (self):
        self._infile.close()

    def shortName (self):
        return self._shortName

    def ZMW (self):
        return self._ZMW

    def numZMWs (self):        # N = number of ZMWs reported in file, not necessarily numbered 0..N-1
        return len(self._ZMW["HoleNumber"])

    def holeNumbers (self):
        '''Return HDF5 array of hole numbers which actually exist in /PulseData/BaseCalls/ZMW/HoleNumber.'''
        return self._ZMW["HoleNumber"]

    def hasConsensus (self):
        return self._hasConsensus

    def findCCSFile (self):
        '''Given a directory to look in, find the ccs.h5 file that contains consensus reads for this bax file.'''

        self._hasConsensus  = False                         # until proven otherwise

        if "PulseData/ConsensusBaseCalls" in self._top:     # if this is an older bax file, in contains its own CCS data

            self._consBasecalls = self._top["PulseData/ConsensusBaseCalls"]
            self._consZMW       = self._top["PulseData/ConsensusBaseCalls/ZMW"]
            self._consPasses    = self._top["PulseData/ConsensusBaseCalls/Passes"]
            self._hasConsensus  = True

        elif self._CCSDir is not None:

            CCSFilename   = os.path.basename(self._filename).replace('bax', 'ccs')
            fqCCSFilename = os.path.join(self._CCSDir, CCSFilename)

            if os.path.exists(fqCCSFilename):

                self._CCSFile = h5py.File (fqCCSFilename, 'r')
                self._consBasecalls = self._CCSFile["PulseData/ConsensusBaseCalls"]
                self._consZMW       = self._CCSFile["PulseData/ConsensusBaseCalls/ZMW"]
                self._consPasses    = self._CCSFile["PulseData/ConsensusBaseCalls/Passes"]
                self._hasConsensus  = True
                logger.debug('BaxFile %s found CCS file %s' % (self._shortName, fqCCSFilename))

            else:
                logger.warning('%s: no CCS file found corresponding to %s' % (self._shortName, self._filename))

        else:
            logger.info('BaxFile %s does not contain CCS data (rel 2.1.0 and later). Use --ccs' % self._shortName)

    def ZMWSanityClause (self):
        '''The trouble with a lot of code today is that there ain't no sanity clause.'''

        # Check that the subfields of PulseData/BaseCalls/ZMW and
        # PulseData/ConsensusBaseCalls/ZMW match, except for
        # NumEvent. This is only a sanity check, and (assuming it
        # passes) has no effect on the class's function. I.e., it
        # doesn't need to be called, if you're the trusting sort.

        # This routine operates on the data of a single bax file, so
        # it remains a method of the BaxFile class.

        if not self._sanityChecked:               # we only need to check once

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
####            logger.debug('%s sanity check passed' % self._shortName)

        return                      # nothing returned -- if we return at all, we passed!


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
