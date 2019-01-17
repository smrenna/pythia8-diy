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
    op.add_option("-f", "--FORCE", dest="FORCE", action="store_true", default=False, help="Overwrite output file if it exists (default: %default)")
    opts, args = op.parse_args()


    NP = [nPart(f)  for f in args]
    NE = [nEvent(f) for f in args]
    NPtot = sum(NP)
    NEtot = sum(NE)

    # Open existing H5 for reading and writing
    #
    f = h5py.File(opts.OUTFILE, "w")

    # This copies the structure
    g = h5py.File(args[0])
    struct = []
    g.visit(struct.append)
    ds_I = [x for x in struct if x.startswith("index")    and "/" in x]
    ds_e = [x for x in struct if x.startswith("event")    and "/" in x]
    ds_p = [x for x in struct if x.startswith("particle") and "/" in x]
    ds_i = [x for x in struct if x.startswith("init")     and "/" in x]
    ds_P = [x for x in struct if x.startswith("procInfo") and "/" in x]
    for _I in ds_I: f.create_dataset(_I, (NEtot,),   dtype=g[_I].dtype)
    for _e in ds_e: f.create_dataset(_e, (NEtot,),   dtype=g[_e].dtype)
    for _p in ds_p: f.create_dataset(_p, (NPtot,),   dtype=g[_p].dtype)
    for _i in ds_i: f.create_dataset(_i, data=g[_i], dtype=g[_i].dtype)
    for _P in ds_P: f.create_dataset(_P, data=g[_P], dtype=g[_P].dtype)
    g.close()


    for num, fname in enumerate(args):
        if num ==0:
            e_start = 0
            p_start = 0
        else:
            e_start=NE[num]
            p_start=NP[num]
        e_end = e_start + NE[num]
        p_end = p_start + NP[num]
        print ("events from {} to {}".format(e_start, e_end))
        print ("particles from {} to {}".format(p_start, p_end))
        g = h5py.File(fname)
        for _e in ds_e: f[_e][e_start:e_end] = g[_e]
        for _p in ds_p: f[_p][p_start:p_end] = g[_p]
        for _I in ds_I: f[_I][e_start:e_end] = g[_I][:]+ p_start
        g.close()



    # from IPython import embed
    # embed()
    f.close()
    sys.exit(0)

    # # Get number of events in file
    # nevents = f['index/start'].size

    # if "npLO" in f["event"].keys() and not opts.FORCE:
        # print("Error, npLO already exists. To overwrite, use -f")
        # import sys
        # sys.exit(1)
    # if "npNLO" in f["event"].keys() and not opts.FORCE:
        # print("Error, npNLO already exists. To overwrite, use -f")
        # import sys
        # sys.exit(1)

    # if "npLO" in f["event"].keys() and opts.FORCE:
        # f["event/npLO"][0:nevents] = npLO *np.ones(nevents, dtype=int)
    # if "npNLO" in f["event"].keys() and opts.FORCE:
        # f["event/npNLO"][0:nevents] = npNLO *np.ones(nevents, dtype=int)

    # if opts.FORCE:
        # f.close()
        # sys.exit(0)

    # f.create_dataset("event/npLO",  data=npLO *np.ones(nevents, dtype=int), compression=opts.COMPRESSION)
    # f.create_dataset("event/npNLO", data=npNLO*np.ones(nevents, dtype=int), compression=opts.COMPRESSION)
    # f.close()
    # sys.exit(0)
