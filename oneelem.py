#!/usr/bin/python

from pylab import *
import sys

if len(sys.argv) != 2:
    print "usage: ./plot.py input_file"
    sys.exit(2)
    
filename = sys.argv[1]
f1 = open(filename, 'r')

data = [map(float, line.split()) for line in f1] 
time = array([bit[0] for bit in data])
dist = array([bit[2] for bit in data])
sma = array([bit[3] for bit in data]) 
ecc = array([bit[4] for bit in data])
inc = (pi/180.0)*array([bit[5] for bit in data])
peri = array([sma[i]*(1.0-ecc[i]) for i in range(len(sma))])
f1.close

figure()
subplot(121)
plot(time,sma,label="semi-major axis")
plot(time, peri,label="pericenter")
plot(time, dist,label="dist")
xlabel('time (yrs)')
ylabel('distance (AU)')
legend(loc=2)

subplot(122)
plot(time,inc,label="inc (rad)")
plot(time, ecc,label="eccentricity") 
xlabel('time (yrs)')
legend(loc=5)

show()
