#!/usr/bin/env python

# Package for accessing data from PacBio cmp.h5 alignment files.

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# From /AlnInfo/AlnIndex annotations (in early versions of cmp.h5 files):

#    ColName00 = AlignmentId
#    ColName01 = ReadGroupId
#    ColName02 = MovieId
#    ColName03 = RefSeqId
#    ColName04 = tStart
#    ColName05 = tEnd
#    ColName06 = RCRefStrand
#    ColName07 = HoleNumber
#    ColName08 = SetNumber
#    ColName09 = StrobeNumber
#    ColName10 = SubreadId
#    ColName11 = rStart
#    ColName12 = rEnd
#    ColName13 = MapQV
#    ColName14 = nM
#    ColName15 = nMM
#    ColName16 = nIns
#    ColName17 = nDel
#    ColName18 = offset_begin
#    ColName19 = offset_end
#    ColName20 = nBackRead
#    ColName21 = nReadOverlap

import sys
import h5py
from tt_log import logger

INDEX_COLS = ('AlignmentId',
              'ReadGroupId',
              'MovieId',
              'RefSeqId',
              'tStart',
              'tEnd',
              'RCRefStrand',
              'HoleNumber',
              'SetNumber',
              'StrobeNumber',
              'SubreadId',
              'rStart',
              'rEnd',
              'MapQV',
              'nM',
              'nMM',
              'nIns',
              'nDel',
              'offset_begin',
              'offset_end',
              'nBackRead',
              'nReadOverlap')

class CmpFile (object):

    def __init__ (self, filename, set=1, strobe=0, maxHole=None):

        logger.debug("creating CmpFile object for set %d strobe %d" % (set, strobe))

        self._filename     = filename
        self._setNumber    = set
        self._strobeNumber = strobe

        self._infile       = h5py.File (filename, 'r')
        self._top          = h5py.Group (self._infile, '/')
        self._index        = self._top['AlnInfo/AlnIndex']
        self._subreadMap   = None

        if maxHole == None:
            self._maxHole = max(self._index[:,7])     # largest *mapped* hole (for any set), may not be max hole!
        else:
            self._maxHole = maxHole

    def __del__ (self):
        self._infile.close()

    def index (self):
        return self._index

    def getSubreadMap (self):
        '''Create and cache a lookup table which maps a ZMW number to its mapped subreads'''

        # Entries in AlnInfo/AlnIndex are ordered by alignment
        # position, so the subread entries for a given ZMW will not be
        # contiguous. Here we construct a list, indexed by ZMW#, of
        # lists containing indexes into AlnInfo/AlnIndex for all
        # subreads for the ZMW.

        if self._subreadMap == None:

            logger.debug("creating subread map")
            
            ishape = self._index.shape
            self._subreadMap = [None] * (self._maxHole+1)

            set    = self._setNumber
            strobe = self._strobeNumber
            kept   = 0

            for ix in xrange(ishape[0]):

                # TODO: create multi-level index for all sets/strobes in one go

                if self._index[ix,8] == set and self._index[ix,9] == strobe:

                    kept += 1

                    hole = self._index[ix,7]
                    if self._subreadMap[hole] == None:
                        self._subreadMap[hole] = [ix]
                    else:
                        self._subreadMap[hole].append(ix)

            logger.debug("processed %d entries, kept %d" % (ishape[0], kept))

        return self._subreadMap

    def getAlignmentIndex (self, hole, start, end):
        '''Return the index of the alignment which falls within start-end for specified hole.'''

        map = self.getSubreadMap()

        if map[hole] != None:

            index = self._index

            for ix in map[hole]:            # find the alignment record for this region
                if index[ix,11] >= start and index[ix,12] <= end:
                    return ix

        return None

    def getAlignmentAsDict (self, hole, start, end):
        '''Return info for the specified alignment as a dict.'''

        ix = self.getAlignmentIndex (hole, start, end)

        if ix is None:
            return None

        index = self._index

        return dict(zip(INDEX_COLS, index[ix,:]))

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
