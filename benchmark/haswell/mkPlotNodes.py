import numpy as np
import pylab

import sys
D=np.loadtxt(sys.argv[1])

X = D[:,0]
Y = D[:,1]
perfect=[Y[0]/2**i for i in range(len(X))]
# Y2 = D[:,2]
# Y3 = D[:,3]
for x in X:
    pylab.axvline(x, color="k", linewidth=0.1)
pylab.xscale("log")
pylab.yscale("log")
# pylab.axvspan(641, X.max(), facecolor='g', alpha=0.5)
pylab.plot(X,Y, label="batch size = 1000")
pylab.plot(X,perfect, "k--", label="perfect scaling")
# pylab.plot(X,Y2, label="batch size = 10000")
# pylab.plot(X,Y3, label="batch size = 100000")
pylab.xlim((X.min(), X.max()))
pylab.xlabel("# nodes")
pylab.ylabel("time/s")
pylab.margins(0.2)
pylab.subplots_adjust(bottom=0.20)
pylab.legend()

pylab.title("Scaling behaviour of processing 1 billion particles on HASWELL")

pylab.xticks(X, [str(int(x)) for x in X])#, rotation='vertical')
pylab.savefig(sys.argv[2])
