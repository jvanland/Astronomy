#!/usr/local/bin/python

from pylab import *
import sys

pi = 3.141592653589793

if len(sys.argv) != 2:
    print "usage: ./multistepelem.py input_file"
    sys.exit(2)
    
filename = sys.argv[1]
f1 = open(filename, 'r')
data1 = [map(float, line.split()) for line in f1] 
sma = array([bit[0] for bit in data1]) 
ecc = array([bit[1] for bit in data1])
peri = array([bit[2] for bit in data1])
mass = array([bit[3] for bit in data1])
dist = array([bit[4] for bit in data1])
energy = mass/(2.0*sma)
#inc = (180.0/pi)*array([bit[2] for bit in data1])
#inc = fmod(inc,360)
#lan = 180+(180.0/pi)*array([bit[3] for bit in data1])
#lan = fmod(lan,360)
#lop = 180+(180.0/pi)*array([bit[4] for bit in data1])
#lop = fmod(lop,360)
#energy = array([bit[5] for bit in data1])
f1.close
time = arange(0,10*len(sma),10)

fig1 = figure()
suptitle(filename)
subplot(511)
plot(time,sma)
xlabel('time'), ylabel('a')
subplot(512)
plot(time, ecc)
xlabel('time'), ylabel('e')
subplot(513)
plot(time, peri)
xlabel('time'), ylabel('peri')
subplot(514)
plot(time, mass)
xlabel('time'), ylabel('Mass')
subplot(515)
plot(time, dist)
xlabel('time'), ylabel('distance')
#subplot(616)
#plot(time, energy)
#xlabel('time'), ylabel('Energy')

show()
