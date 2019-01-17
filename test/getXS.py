#!/usr/bin/env python

__doc__="""

%prog [options] inputFile


"""



import numpy as np
import h5py

if __name__ == "__main__":
    import optparse, os, sys
    op = optparse.OptionParser(usage=__doc__)
    op.add_option("-v", "--debug", dest="DEBUG", action="store_true", default=False, help="Turn on some debug messages")
    op.add_option("-n", dest="NEVENTS", type=int, default=10000, help="Number of events to analyse (default: %default)")
    opts, args = op.parse_args()


    # Open existing H5 for reading and writing
    f = h5py.File(args[0])

    print("XS: {}".format(sum(f["event/weight"][0:opts.NEVENTS])/sum(f["event/trials"][0:opts.NEVENTS])))

    f.close()
    sys.exit(0)
