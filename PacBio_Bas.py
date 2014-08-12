#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Author: Tom Skelly (thomas.skelly@nih.gov)

# Print detailed data from *.bas.h5 file and (optionally) the *.cmp.h5
# alignments file on a region-by-region basis. If a cmp.h5 file is
# supplied, output includes the alignment outcome for each region.

# CCS reads are also printed if they exist in the *.bas.h5 file, or if
# (starting with SMRTanalysis release 2.1.0) the location of the
# separate *.ccs.files where they reside is supplied as an additional
# parameter (see --help).

# NOTE: For help in interpreting the various line types printed by
# this script, do "PacBio_Bas.py --linehelp". Or see the lineHelp
# method below.

import sys
import optparse
import math          # for log10
import H5BasFile
import H5CmpFile
import SWAligner
from tt_log import logger

# Defaults for command line parameters -- see getParms()

DEF_SCORE_THRESHOLD  = 750
DEF_HQ_LENGTH        = 50
DEF_ADAPTER_LENGTH   = 45

def main ():

    logger.debug("%s starting" % sys.argv[0])

    opt, args = getParms()

    basFilename = args[0]
    logger.debug("bas file: %s" % basFilename)
    bf = H5BasFile.BasFile (basFilename, CCSDir=opt.ccs)

    if not opt.nocons:
        if not bf.hasConsensus():
            logger.warning('no ccs data found: turning on --nocons')
            opt.nocons = True

    cmp = None                      # no cmp file?

    if opt.aln is not None:         # was a subread cmp.h5 file specified?

        cmpFilename = opt.aln
        logger.debug("cmp file: %s" % cmpFilename)
        cf  = H5CmpFile.CmpFile (fileName=cmpFilename)
        cmp = H5CmpFile.CmpMovie (cmpObject=cf,
                                  movieName=bf.movieName(),
                                  maxHole=bf.maxZMW())

    cmpCCS = None

    if opt.alnccs is not None:      # was a CCS cmp.h5 file specified?

        cmpCCSFilename = opt.alnccs
        logger.debug("CCS cmp file: %s" % cmpCCSFilename)
        cfCCS  = H5CmpFile.CmpFile (fileName=cmpCCSFilename)
        cmpCCS = H5CmpFile.CmpMovie (cmpObject=cfCCS,
                                     movieName=bf.movieName(),
                                     maxHole=bf.maxZMW())

    aln = SWAligner.Aligner()           # we'll use this in the loop below for finding adapters
    aln.setRead (H5BasFile.ADAPTER)     # adapter sequence is query
    minAdapterScore = opt.adapter * aln.getPenalties()[0] / 2

    print "   ZMW     b/s stat prod tp  start end+1    len  aln  chr st",
    print "     from         to   off  astart  aend+1   mm  ins del    Q"
    print

    for hole in bf.holeNumbers():          # main loop!

        numBases = bf.readLen(hole)
        zStat    = bf.holeStatusStr(hole)  # this is a string, not a number
        zProd    = bf.productivity(hole)

        numGoodInserts = 0

        HQStart, HQEnd, HQScore = bf.HQregion(hole)[2:5]

        for region in bf.holeRegions(hole):

            regionHole, regionType, start, end, score = region

            inHQ = end > HQStart and start < HQEnd             # does region overlap HQ?

            regionDuration = float(bf.elapsedFrames(hole, start, end)) / H5BasFile.frameRate
            regionBps = (end-start) / regionDuration if regionDuration > 0 else 0

            if regionType == 0:                                # an adapter region?

                if not opt.noadapt:                            # if we are printing adapter lines
                    print "%6d  %6.3f  %-5s  %d"  % (hole, regionBps, zStat, zProd),    # these appear in every line
                    flag = 'A ' if inHQ else 'a '
                    print "%-2s  %5d %5d"  % (flag, start, end)

            elif regionType == 2:                              # a HQ region?

                print "%6d  %6.3f  %-5s  %d"  % (hole, regionBps, zStat, zProd),        # these appear in every line

                if zProd != 1 or not bf.isSequencingZMW(hole) or (HQEnd-HQStart) < opt.length:
                    flag = 'h '
                else:
                    flag = 'H+' if score >= opt.score else 'H '

                print "%-2s  %5d %5d"  % (flag, start, end),

                readDuration = float(bf.elapsedFrames(hole)) / H5BasFile.frameRate
                readBps = numBases/readDuration if readDuration > 0 else 0
                print "            score: %3d  HQ: %5d  read: %5d  dur: %8.3f  b/sec: %6.3f" \
                    % (score, HQEnd-HQStart, numBases, readDuration, readBps)

            elif regionType == 1:                              # a subread?

                print "%6d  %6.3f  %-5s  %d"  % (hole, regionBps, zStat, zProd),        # these appear in every line

                insSize = end - start

                align = None
                if cmp is not None:                                       # if a cmp.h5 was supplied
                    align = cmp.getAlignmentByPosition (hole, start, end) # alignment record for this region

                if align is not None:                                     # if the region aligned

                    numGoodInserts += 1

                    flag = 'I+'
                    print "%-2s  %5d %5d  %5d"  % (flag, start, end, insSize),

                    rStart, rEnd = align['rStart'], align['rEnd']   # fetch once, used many times
                    alnLen = rEnd-rStart
                    nMM, nIns, nDel = align['nMM'], align['nIns'], align['nDel']

                    print "%5d  %2d %1s  %9d  %9d  %4d  %6d  %6d  %3d %4d %3d %4.1f" % \
                        (alnLen,                               # length of aligned portion of read
                         align['contig'],                      # chr/contig id (see H5CmpFile)
                         '-' if align['RCRefStrand'] else '+', # strand
                         align['tStart'], align['tEnd'],       # reference offset of start/end of alignment
                         rStart-start,                         # offset of alignment start into insert
                         rStart, rEnd,
                         nMM, nIns, nDel,                      # # of mismatches, insertions, deletions
                         getQ (align)),                        # read quality Q score for insert

                elif insSize > opt.adapter * 2:                # if it's a non-descript, non-aligned region
                                                               # TODO: Make the '2' a parameter
                    numGoodInserts += 1

                    flag = 'I ' if inHQ else 'i '
                    print "%-2s  %5d %5d  %5d"  % (flag, start, end, insSize),

                elif insSize < opt.adapter:                    # if it's too short to be an adapter

                    flag = 'Is' if inHQ else 'is'
                    print "%-2s  %5d %5d  %5d"  % (flag, start, end, insSize),

                else:                                          # see if it's really an adapter that wasn't called

                    sequence = bf.getSequence(hole, start, end)
                    aln.setRef (sequence)
                    alnScore = aln.fillMatrix()                # align it to adapter

                    flag = 'ia' if alnScore >= minAdapterScore else 'is'
                    print "%-2s  %5d %5d  %5d  %2d"  % (flag, start, end, insSize, alnScore),

                    if alnScore >= minAdapterScore:
                        peaks = aln.peakPosits()
                        print " %2d" % len(peaks),
                        refString, readString = aln.alignmentStrings()
                        print ' ', refString                                                 # EOL here
                        print "                         i2                              ",   # new line
                        print readString,

                print

            else:
                raise ValueError ("unrecognised region type %d in ZMW %d" % (regionType, hole))

        if not opt.nocons:      # note that opt.nocons gets turned on if no CCS data is found
            if zProd == 1 and bf.isSequencingZMW(hole):
                printCCSDataForHole (bf, hole, numGoodInserts, cmpCCS)       # process consensus read passes

    logger.debug("complete")

def printCCSDataForHole (bf, hole, numGoodInserts, cmpCCS):
    '''Print consensus read info for specified hole.'''

    zStat   = bf.holeStatusStr(hole)  # this is a string, not a number
    zProd   = bf.productivity(hole)
    consLen = bf.consReadLen(hole)

    for passDict in bf.holeConsensusPasses(hole):

        start    = passDict["PassStartBase"]
        numBases = passDict["PassNumBases"]
        end      = start + numBases

        before = 'B' if passDict["AdapterHitBefore"] else ' '
        after  = 'A' if passDict["AdapterHitAfter"]  else ' '
        dir    = '-' if passDict["PassDirection"]    else '+'

        print "%6d          %-5s  %d" % (hole, zStat, zProd),
        print "c   %5d %5d  %5d        %s%s %s" % (start, end, numBases, before, after, dir)

    if consLen > 0:                 # if there is a CCS read

        numP = bf.numConsensusPasses(hole)
        avgQ = bf.getConsensusAverageQuality(hole)

        align = None
        if cmpCCS is not None:

            alignments = [ a for a in cmpCCS.getAlignmentsForHole (hole) ]     # get all alignments from generator
            numAlign = len(alignments)
            if numAlign > 1:
                raise RuntimeError ('more than one CCS alignment found for hole %d' % hole)
            if numAlign == 1:
                align = alignments[0]

        if align is None:
            print "%6d          %-5s  %d" % (hole, zStat, zProd),
            print "CC               %5d  %2d  %2d" % (consLen, numP, avgQ)
        else:

            print "%6d          %-5s  %d" % (hole, zStat, zProd),
            print "C+               %5d  %2d  %2d" % (consLen, numP, avgQ),

            rStart, rEnd = align['rStart'], align['rEnd']
            alnLen = rEnd-rStart
            nMM, nIns, nDel = align['nMM'], align['nIns'], align['nDel']

            print "%5d  %2d %1s  %9d  %9d  %6d  %6d  %3d %4d %3d %4.1f" % \
                (alnLen,                               # length of aligned portion of read
                 align['contig'],                      # chr/contig id (see H5CmpFile)
                 '-' if align['RCRefStrand'] else '+', # strand
                 align['tStart'], align['tEnd'],       # reference offset of start/end of alignment
                 rStart, rEnd,                         # CCS read offset of start/end of alignment
                 nMM, nIns, nDel,                      # # of mismatches, insertions, deletions
                 getQ (align))                         # read quality Q score for insert

    elif numGoodInserts > 2:                           # why didn't this guy get a CCS read?
        print "%6d          %-5s  %d C?" % (hole, zStat, zProd)

    return

def getQ (align):
    '''Compute alignment quality as a Q score.'''

    # Compute the Q score for the aligned insert, according to the
    # (rather strict) PacBio formulation where each mismatched,
    # inserted or deleted base is an error.

    rStart, rEnd = align['rStart'], align['rEnd']
    alnLen = rEnd-rStart
    nMM, nIns, nDel = align['nMM'], align['nIns'], align['nDel']
    totErrors = nMM + nIns + nDel
    if totErrors == 0:
        Q = 40.0
    elif alnLen <= 0:
        Q = 0.0
    else:
        Q = -10.0 * math.log10 (float(totErrors) / float(alnLen))

    return Q

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <bas_file>')

    parser.add_option ('--linehelp', action='store_true', help='show line-specific help and exit')
    parser.add_option ('--filehelp', action='store_true', help='show help for file-related parameters, and exit')
    parser.add_option ('--ccs',                  help='directory containing ccs.h5 files for CCS reads, post-2.1.0')
    parser.add_option ('--aln',                  help='cmp.h5 file for subread alignments')
    parser.add_option ('--alnccs',               help='cmp.h5 file for CCS alignments')
    parser.add_option ('--score',    type='int', help='minimum HQ region score (def: %default)')
    parser.add_option ('--length',   type='int', help='minimum HQ region length (def: %default)')
    parser.add_option ('--adapter',  type='int', help='expected adapter length (def: %default)')
    parser.add_option ('--noadapt',  action='store_true', help='do not print adapter lines (for brevity)')
    parser.add_option ('--nocons',   action='store_true', help='do not print consensus passes lines')

    parser.set_defaults (score=DEF_SCORE_THRESHOLD,
                         length=DEF_HQ_LENGTH,
                         adapter=DEF_ADAPTER_LENGTH)

    opt, args = parser.parse_args()

    if opt.linehelp:
        lineHelp()

    if opt.filehelp:
        fileHelp()

    if opt.linehelp or opt.filehelp:
        sys.exit()

    if len(args) > 1:
        logger.warning ('WARNING: alignments cmp.h5 file should now be specified with --aln keyword')
        opt.aln = args.pop()      # put it where it belongs

    return opt, args

def fileHelp ():

    print ("""
This script reads several types of file, as described here. The basic
output of the script includes data from the bas/bax files. CCS reads,
subread alignments, and CCS read alignments can also be included by
specifying the relevant files with command-line keywords.

      *.bas.h5 -- specified as command-line argument (no keyword)

         The top of the tree. This file used to contain the data
         produced by primary analysis on the instrument. That data has
         since been progressively federated to other files, so the
         bas.h5 file is now just an anchor for accessing the real data
         in those other files.

      *.[123].bax.h5 -- found via bas.h5 file

         The primary analysis data which used to reside in the bas.h5
         file now lives in three bax.h5 files. Each bax contains data
         for one third of the ZMWs. The bas.h5 file points to its bax
         children.

      *.[123].ccs.h5 -- found in directory specified by --ccs keyword

         Prior to release 2.1.0, CCS reads were produced by primary
         analysis processing on the instrument, and were found in the
         bax.h5 files. With 2.1.0, CCS creation has been moved to the
         offline secondary analysis step, under control of a protocol
         such as ReadOfInsert. Three ccs.h5 files are created,
         corresponding to the three bax.h5 files.

         If you want ccs data to be included in this script's output,
         specify the secondary analysis directory where the the ccs.h5
         files reside, via the --ccs command line keyword.

      aligned_reads.cmp.h5 (subreads) -- specified by --aln keyword

         A typical step in secondary analysis processing is to invoke
         the blasr aligner to align subreads to a reference. The
         output of that processing resides in a cmp.h5 file (normally
         named aligned_reads.cmp.h5).

         If you want alignment results for each subread to be included
         in this script's output, specify the cmp.h5 file, via the
         --aln command line keyword.

      aligned_reads.cmp.h5 (ccs reads) -- specified by --alnccs keyword

         CCS reads, once produced, can also be aligned to a reference,
         using secondary analysis protocols such as Resequencing_ReadsOfInsert.
         The filename for CCS alignments is also normally aligned_reads.cmp.h5.

         If you want alignment results for each CCS to be included in
         this script's output, specify the cmp.h5 file, via the
         --alnccs command line keyword. For data from release 2.1.0
         and later, you will also need to specify the ccs.h5 file
         containing those reads, via the --ccs command line keyword.
          """)

def lineHelp ():

    print ("""
 Region line type flags:

     HQ: H+  High-scoring HQ region
         H   Low-scoring
         h   non-SEQ or not prod-1 or too short

     IN: I+  Aligned insert
         I   Non-aligned insert within HQ region
         i   Insert outside HQ region, but otherwise legitimate
         ia  Short (40-100) insert which looks like an adapter
         i2  Adapter as aligned to insert (follows ia)
         is  Short (0-100) insert, not adapter-like

     AD: A   Adapter within HQ region
         a   Adapter outside HQ region

     CO: c   Consensus read pass
         CC  Consensus read length
         C?  Missing consensus read

 Insert line columns:

     ZMW     ZMW number
     b/s     bases per second for this region
     stat    ZMW status (e.g., SEQUENCING)
     prod    Productivity (0=empty 1=good, 2=overloaded)
     tp      Region line type flag, as defined above
     start   0-based offset of region in read
     end+1   End offset+1 of region in read
     len     Length of region (for inserts)
     aln     Number of aligning bases in region
     chr     Aligning chr/contig name
     st      Aligning strand
     from    1-based reference coordinate of start of alignment
     to      End of alignment
     off     Offset of start of alignment within insert
     astart  Offset of start of alignment within read
     aend+1  End offset+1 of alignment within read
     mm      Number of mismatches in alignment
     ins     Number of insertions in alignment
     del     Number of deletions in alignment
     Q       Q-score of alignment: -10*log((mm+ins+del)/aln)

 High-Quality region fields:

     score:  Score of HQ region
     HQ:     Length of HQ region
     read:   Length of read
     dur:    Duration of read in seconds
     b/sec:  bases per second for total read
          """)

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
