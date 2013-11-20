#!/usr/bin/env python

# Copyright (C) 2012 Genome Research Limited -- See full notice at end
# of module.

# Package for accessing data from PacBio cmp.h5 alignment files.

# A bit of design philosophy: The front door to a cmp file is the
# /AlnInfo/AlnIndex dataset. It contains one row for each aligned
# sub-read. Each row includes 22 columns as shown below. Since a cmp
# file can be sorted in several ways, no assumptions can be made about
# the order of the entries.

# So how should this module provide user access to AlnIndex data? We
# could implement that as a two-step interface, where the user must
# first determine the indexes of the entries of interest, then request
# the data for each of those entries. However, having determined an
# index, the user's next action will almost certainly be to request
# the data for it. So here we combine the two steps: this module's
# methods for accessing AlnIndex return the data itself, rather than
# an index to it. This eliminates the need for the user to know about
# entry indexes at all.

# There are several such access methods, mostly implemented as
# generator functions. See the doc strings for getAllAlignments,
# getAlignmentsByHole, getAlignmentsForHole and getAlignmentByPosition.

# The cmp file also contains numerous groups with names like
# /refnnnnn/<moviename>. That's where the basecall-level information
# resides. The mapping from an alignment entry's ReadGroupId to the
# basecall group is via /AlnGroup/Path. The getReadGroups method
# manages that mapping.

# NOTE that as of Feb 2012 the interface to the cmp file provided by
# this module has changed. Remember that a cmp file can contain
# alignments for many movies. Prior to the change, a CmpFile object
# could access data from only one movie. In most cases, that was quite
# adequate since, in most cases, the cmp file is being accessed in
# combination with a bas.h5 basecall file, which contains data for one
# movie only. However, it is possible that a user would like to access
# data from multiple movies. So we now provide a double-barrelled
# interface which defines file and movie objects. We've added more
# general functionality, at the cost of increased complexity both here
# and for the user.

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

# ALSO: Nothing in the above list is the actual chromosome/contig
# number of the alignment. Given that the caller will probably want
# that, we sneak it in in getAlignmentAsDict, as key 'contig' in the
# returned dictionary. See getAlignmentAsDict for the details.

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

BASE_MAP = ('-',  'A',  'C', None,       # encoding for alnArray 4-bit subfields
            'G', None, None, None,
            'T', None, None, None,
            None, None, None, 'N',
            )

class CmpFile (object):
    '''File-level interface to a cmp.h5 alignments file.'''


    def __init__ (self, fileName):

        logger.debug("creating CmpFile object for %s" % (fileName))

        self._fileName = fileName
        self._infile   = h5py.File (fileName, 'r')
####        self._top      = h5py.Group (self._infile, '/')
        self._top      = self._infile         # h5py 2.0.1 change!

        return

    def __del__ (self):
        self._infile.close()

    def top (self):
        return self._top

    def movieList (self):
        '''Generator function to return the list of movie names included in the cmp file.'''

        list = self._top['MovieInfo/Name']
        for movie in list:
            yield movie

        return

    def printDetails (self):
        '''Debug routine to print a bunch of stuff from the file to the log. Not production grade.'''

        movName = self._top['MovieInfo/Name']
        movID   = self._top['MovieInfo/ID']
        for ix in xrange(len(movName)):
            logger.debug("movie %d (%d): %s" % (ix, movID[ix], movName[ix]))

        path  = self._top['AlnGroup/Path']
        ID    = self._top['AlnGroup/ID']
        for ix in xrange(len(path)):
            logger.debug("path %d (%d): %s" % (ix, ID[ix], path[ix]))

####        logger.debug("" % ())

        return

    # It's not clear that any other methods of CmpFile would be
    # useful. We could do a ReadGroup-to-MovieName method. But if the
    # user knows the ReadGroup, he presumably already knows the
    # MovieName -- that's how he got there.

class CmpMovie (object):
    '''Access to data from one movie in a cmp.h5 file.'''

    def __init__ (self, cmpObject, movieName, maxHole=None):

        logger.debug("creating CmpMovie object for %s" % (movieName))

        self._movieName    = movieName
        self._top          = cmpObject.top()
        self._index        = self._top['AlnInfo/AlnIndex']
        self._refMap       = self._top['RefGroup/RefInfoID']
        self._subreadMap   = None
        self._readGroups   = None

        if maxHole is None:
            self._maxHole = max(self._index[:,7])     # largest *mapped* hole (for any set), may not be max hole!
        else:
            self._maxHole = maxHole

        # Search the /MovieInfo group looking for the movie of
        # interest, and save its ID. We'll use that later to select
        # alignments for this movie from /AlnInfo/AlnIndex.

        movName = self._top['MovieInfo/Name']
        movID   = self._top['MovieInfo/ID']
        self._movieID = None
        for ix in xrange(len(movName)):
            if movName[ix] == self._movieName:
                self._movieID = movID[ix]
                logger.debug("movie ID is %d" % movID[ix])
                break
        if self._movieID is None:
            logger.debug("movie not found in /MovieInfo/Name")
            raise RuntimeError

        # Should we call getSubreadMap at this point? It's fairly
        # time-consuming. Let's put it off until it's actually needed.

        return

    def movieName (self):
        return self._movieName

    def getAlignmentAsDict (self, ix):      # access to AlnIndex info, given an ix
        '''Return info for the specified alignment as a dict.'''

        if ix is None:
            return None

        index = self._index

        ret = dict(zip(INDEX_COLS, index[ix,:]))

        # In the alignment record, the actual chromosome/contig
        # number, as defined by the reference file (the 'refNNNNNN'
        # designator in the contig name, when reprocessed by
        # SMRTportal), is not actually present. ReadGroupId (col 1) is
        # an index into /AlnGroup/Path (see getReadGroups). Reads from
        # two different movies mapping to the same chr will have
        # different ReadGroupIds. And the order of chrs for two movies
        # is not the same. RefSeqId (col 3) *will* be the same for
        # reads mapping the same chr, and should be used as an index
        # into RefGroup/RefInfoID to get the actual chr number.

        # Rather than confront the user with the need to understand
        # all that, we'll sneak an extra entry into the dict we create
        # here. Key 'contig' will contain the chr number.

#### All the above appears to have been invalidated in release 1.4.0. 

####        ret['contig'] = self._refMap[ret['RefSeqId']-1]     # RefSeqId runs 1..N, RefGroup/RefInfoID is 0..N-1
        ret['contig'] = ret['RefSeqId']

        return ret

    def getAllAlignments (self):
        '''Generator function to return all indexes, in whatever order the file is sorted in.'''

        ishape = self._index.shape

        for ix in xrange(ishape[0]):
            yield self.getAlignmentAsDict (ix)

    def getAlignmentsByHole (self):
        '''Generator function to return all alignments, sorted by hole/read position.'''

        map   = self.getSubreadMap()
        index = self._index

        for hole in xrange(self._maxHole):

            # We could call getAlignmentsForHole here, but quicker to
            # do the business right here, rather than invoke a method.

            if map[hole] != None:

                for ix in sorted(map[hole], key=lambda n: index[n,11]):    # col 11 is rStart
                    yield self.getAlignmentAsDict (ix)

        return

    def getAlignmentsForHole (self, hole):
        '''Generator function to return alignments for a hole, in order as they appear in the read.'''

        map = self.getSubreadMap()

        if hole <= self._maxHole and map[hole] is not None:

            index = self._index

            for ix in sorted(map[hole], key=lambda n: index[n,11]):    # col 11 is rStart
                yield self.getAlignmentAsDict (ix)

        return

    def getAlignmentByPosition (self, hole, start, end):
        '''Return the alignment which falls within start-end for specified hole.'''

        map = self.getSubreadMap()

        if hole <= self._maxHole and map[hole] is not None:

            index = self._index

            for ix in map[hole]:            # find the alignment record for this region
                if index[ix,11] >= start and index[ix,12] <= end:
                    return self.getAlignmentAsDict (ix)     # return: this is NOT a generator function

        return None

    def getAlignmentStrings (self, align):      # 'align' is an alignment dict
        '''Return read and reference aligments strings from specified alignment.'''

        readGroups = self.getReadGroups()       # list of /refnnnnn/<moviename> groups, indexed by ReadGroupId
        RG         = align['ReadGroupId']
        alnArray   = readGroups[RG]['AlnArray']

        refString  = []
        readString = []

        # An alnArray entry is an 8-bit integer encoding read and references bases in 4 bits each.

        for bases in alnArray[align['offset_begin']:align['offset_end']]:
            refString.append(BASE_MAP[bases & 0x0f])
            readString.append(BASE_MAP[bases>>4])

        return (''.join(refString), ''.join(readString))

    def getSubreadMap (self):
        '''Create and cache a lookup table which maps a ZMW number to its mapped subreads'''

        # Entries in AlnInfo/AlnIndex can be sorted in various ways
        # (e.g., by alignment position), so the subread entries for a
        # given ZMW may not be contiguous. Here we construct a list,
        # indexed by ZMW#, of lists containing indexes into
        # AlnInfo/AlnIndex for all subreads for the ZMW.

        if self._subreadMap is None:

            logger.debug("creating subread map")
            
            ishape = self._index.shape
            self._subreadMap = [None] * (self._maxHole+1)

            readGroups = self.getReadGroups()          # set of read group IDs for this movie

            kept   = 0
            movieID = self._movieID

            for ix in xrange(ishape[0]):               # for all alignments

####                if self._index[ix,1] in readGroups:    # if this is an alignment for the current movie
                if self._index[ix,2] == movieID:    # if this is an alignment for the current movie

                    kept += 1

                    hole = self._index[ix,7]
                    if self._subreadMap[hole] == None:
                        self._subreadMap[hole] = [ix]
                    else:
                        self._subreadMap[hole].append(ix)

            logger.debug("kept %d of %d subread entries for this movie" % (kept, ishape[0]))

        return self._subreadMap

    def getReadGroups (self):
        '''Find read groups in /AlnGroup/Path for this movie.'''

        # Create and cache a dict whose key is ReadGroupId and whose
        # value is h5 group containing the AlnArray and other datasets
        # for that read group. See the comments in getAlignmentAsDict
        # for further confusion on the issue.

        if self._readGroups is None:

            logger.debug("creating ReadGroup list")

            self._readGroups = dict()

            movie = self._movieName
            path  = self._top['AlnGroup/Path']
            ID    = self._top['AlnGroup/ID']

            for ix in xrange(len(path)):
                logger.debug("path %d, ID %d: %s" % (ix, ID[ix], path[ix]))
####                if path[ix].endswith(movie):
####                    self._readGroups[ID[ix]] = self._top[path[ix]]
                self._readGroups[ID[ix]] = self._top[path[ix]]
                    
            logger.debug("kept %d of %d ReadGroups for this movie" % (len(self._readGroups), len(path)))

        return self._readGroups


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
