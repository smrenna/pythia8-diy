#!/usr/bin/env python

__doc__="""

%prog [options] inputDir


Read in hierarchical mc run data just as usual and
convert into HDF5 file.

The following datasets are created:

    * runs  --- not really necessary
    * index --- NBIN bin names used to identify bits
    * params (with attribute names) --- the (NPARAMS * NRUNS) parameter array
    * values --- (NBIN * NRUNS) size array
    * errors --- (NBIN * NRUNS) size array
    * xmin   --- (NBIN * NRUNS) size array
    * xmax   --- (NBIN * NRUNS) size array
"""



import numpy as np
import h5py

if __name__ == "__main__":
    import optparse, os, sys
    op = optparse.OptionParser(usage=__doc__)
    op.add_option("-v", "--debug", dest="DEBUG", action="store_true", default=False, help="Turn on some debug messages")
    op.add_option("-c", dest="COMPRESSION", type=int, default=4, help="GZip compression level (default: %default)")
    opts, args = op.parse_args()

    if len(args)!=3:
        sys.exit(1)

    npLO =int(args[1])
    npNLO=int(args[2])

    # Open existing H5 for reading and writing
    f = h5py.File(args[0], "r+")

    # Get number of events in file
    nevents = f['index/start'].size
    f.create_dataset("event/npLO",  data=npLO *np.ones(nevents, dtype=int), compression=opts.COMPRESSION)
    f.create_dataset("event/npNLO", data=npNLO*np.ones(nevents, dtype=int), compression=opts.COMPRESSION)
    f.close()
    sys.exit(0)
