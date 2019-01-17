import numpy as np
import pylab

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

import sys
D=np.loadtxt(sys.argv[1])

X  = D[:,0]
X2 = np.array([x*68 for x in X])
Y  = D[:,1]


pylab.clf()
print X2
perfect   = [Y[0]/2**i for i in range(len(X2))]
for x in X2:
    pylab.axvline(x, color="k", linewidth=0.1)
pylab.xscale("log")
pylab.yscale("log")
pylab.plot(X2,Y, label="1ms per event, 68 ranks per node")
pylab.plot(X2,perfect, "k--", linewidth=0.1, label="Ideal scaling")
# pylab.plot(X,Y2, label="batch size = 10000")
pylab.xlim((X2.min(), X2.max()))
pylab.xlabel("\# ranks")
pylab.ylabel("time/s")
pylab.margins(0.2)
pylab.subplots_adjust(bottom=0.20)
pylab.legend()

pylab.title("Scaling behaviour of reading 1 billion particles on KNL")

pylab.xticks(X2, ["$%s$"%str(int(x)) for x in X2], rotation='vertical')
pylab.savefig("ranks"+sys.argv[2], transparent=True)




