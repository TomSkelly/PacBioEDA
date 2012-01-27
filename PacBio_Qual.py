#!/usr/bin/env python

# Read a bas.h5 file, print quality information for a specified ZMW#.

import sys
import H5BasFile
from tt_log import logger

def main ():

    logger.debug("%s starting" % sys.argv[0])

    basfileName = sys.argv[1]
    hole        = int(sys.argv[2])

    bf = H5BasFile.BasFile (basfileName)

    call     = bf.getBasecallField("Basecall",        hole)
    delete   = bf.getBasecallField("DeletionQV",      hole)
    wuzzit   = bf.getBasecallField("DeletionTag",     hole)
    insert   = bf.getBasecallField("InsertionQV",     hole)
    prebase  = bf.getBasecallField("PreBaseFrames",   hole)
    subst    = bf.getBasecallField("SubstitutionQV",  hole)
    couldbe  = bf.getBasecallField("SubstitutionTag", hole)
    qual     = bf.getBasecallField("QualityValue",    hole)
    duration = bf.getBasecallField("WidthInFrames",   hole)

    print "Index  Call   Qd  Del   Qi   Qs  Sub    Q               DelT    Dur"
    print

    for ix in xrange(len(call)):

        print "%5d  %4s  %3d  %3s  %3d  %3d  %3s  %3d  %4d  %3d  %6.3f  %5.3f" \
            % (ix, chr(call[ix]), delete[ix], chr(wuzzit[ix]), insert[ix], subst[ix], chr(couldbe[ix]), qual[ix],
               prebase[ix], duration[ix],
               float(prebase[ix])/H5BasFile.frameRate, float(duration[ix])/H5BasFile.frameRate)

    logger.debug("complete")

if __name__ == "__main__":
    main()
