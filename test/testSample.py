#!/usr/bin/env python

__doc__="""

%prog [options] inputfiles

Concatenate hdf5 LHE files as output by Sherpa

TODO: probably add the npLO and npNLO datasets here
so that there is only one postprocessing step.

"""


def nEvent(fname):
    import h5py
    f = h5py.File(fname)
    ret = f['event/weight'].shape[0]
    f.close()
    return ret

def nPart(fname):
    import h5py
    f = h5py.File(fname)
    ret = f['particle/id'].shape[0]
    f.close()
    return ret


import numpy as np
import h5py

if __name__ == "__main__":
    import optparse, os, sys
    op = optparse.OptionParser(usage=__doc__)
    op.add_option("-v", "--debug", dest="DEBUG", action="store_true", default=False, help="Turn on some debug messages")
    op.add_option("-c", dest="COMPRESSION", type=int, default=4, help="GZip compression level (default: %default)")
    op.add_option("-o", dest="OUTFILE", default="merged.h5", help="Output file name (default: %default)")
    op.add_option("-f", "--force", dest="FORCE", action="store_true", default=False, help="Overwrite output file if it exists (default: %default)")
    op.add_option("-n", "--nevents", dest="NEVENTS", type=int, default=100, help="Desired number of events (default: %default)")
    opts, args = op.parse_args()


    # Collect dataset sizes from files
    NE = nEvent(args[0])

    if opts.NEVENTS > NE:
        print("Not enough events, input file has {}".format(NE))
        import sys
        sys.exit(1)

    NE = opts.NEVENTS

    # Open existing new H5 for writing
    f = h5py.File(opts.OUTFILE, "w")

    # This copies the structure
    # NOTE the init and procinfo things realy are attributes and
    # don't change from file to file so we just copy them from the
    # first file to the output file.
    g = h5py.File(args[0])

    NP = g["index/end"][opts.NEVENTS]

    struct = []
    g.visit(struct.append)
    ds_I = [x for x in struct if x.startswith("index")    and "/" in x]
    ds_e = [x for x in struct if x.startswith("event")    and "/" in x]
    ds_p = [x for x in struct if x.startswith("particle") and "/" in x]
    ds_i = [x for x in struct if x.startswith("init")     and "/" in x]
    ds_P = [x for x in struct if x.startswith("procInfo") and "/" in x]
    # Empty dataset but with correct size and dtype
    for _I in ds_I: f.create_dataset(_I, data=g[_I][0:NE],   dtype=g[_I].dtype, compression=opts.COMPRESSION)
    for _e in ds_e: f.create_dataset(_e, data=g[_e][0:NE],   dtype=g[_e].dtype, compression=opts.COMPRESSION)
    for _p in ds_p: f.create_dataset(_p, data=g[_p][0:NP],   dtype=g[_p].dtype, compression=opts.COMPRESSION)
    # Full dataset TODO think about makeing those attributes???
    for _i in ds_i: f.create_dataset(_i, data=g[_i], dtype=g[_i].dtype, compression=opts.COMPRESSION)
    for _P in ds_P: f.create_dataset(_P, data=g[_P], dtype=g[_P].dtype, compression=opts.COMPRESSION)
    g.close()

    f.close()
    sys.exit(0)
