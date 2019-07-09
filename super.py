#!/usr/local/bin/python

from pylab import *
import sys

if len(sys.argv) != 2:
    print "usage: ./elements.py bin_num"
    sys.exit(2)

n = int(sys.argv[1])    

filename = "super%d.dat" % (n)
f1 = open(filename, 'r')
data = [map(float, line.split()) for line in f1] 
step = array([bit[0] for bit in data]) 
a = array([bit[2] for bit in data]) 
e = array([bit[3] for bit in data])
p = array([bit[5] for bit in data])
f1.close
time = (0.1987/(2.0*pi))*step

figure()
plot(time,a,label="c.o.m. a (AU)")
plot(time, p,label="pericenter (AU)") 
title(filename)
xlabel('time (yrs)')
#ylim(0,5)
#xlim(0,160000)
legend(loc=3)

show()
