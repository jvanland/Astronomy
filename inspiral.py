#!/usr/bin/python

from pylab import *
import os

f2 = open('inspiral.dat','w')
f3 = open('separate.dat','w')

def neg(array):
    for index, item in enumerate(array):
        if item < 0.0:
            return index

for i in [2,3,4,5]:
    for j in range(4):
        for k in range(10):

            print i, j, k+1
            name = "kozai10k_%d.%d_%d.dat" % (i, j, k+1)
            if os.path.exists(name):
                f1 = open(name,'r')
                data = [map(float, line.split()) for line in f1]
                acm = array([bit[4] for bit in data])
                pbin = array([bit[3] for bit in data])
                ebin = array([bit[2] for bit in data])
                abin = array([bit[1] for bit in data])
                time = array([bit[0] for bit in data])
                f1.close
                pmin = min(pbin)
                emin = max(ebin)
                amin = min(abin)
                amax = max(abin)
                afin = acm[-1]
                atime = time[abin.argmin()]
                ftime = time[-1]
                deltat = ftime - atime
                stime = time[neg(abin)]
                if (amin > 0.0 and amin < 10.0):
                    print >> f2, i, j, k+1, atime, ftime, pmin, emin, amin, amax, afin
                if (amin < 0.0):
                    print >> f3, i, j, k+1, stime, ftime, pmin, emin, amin, amax
f2.close
f3.close
        

