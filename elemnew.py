#!/usr/local/bin/python

from pylab import *
import sys

f1 = open('1elem3931_2.dat', 'r')
data1 = [map(float, line.split()) for line in f1] 
time = array([bit[0] for bit in data1]) 
acm = array([bit[1] for bit in data1]) 
ecm = array([bit[2] for bit in data1])
icm = array([bit[3] for bit in data1])
lancm = array([bit[4] for bit in data1])
lopcm = array([bit[5] for bit in data1])
#meacm = array([bit[6] for bit in data1])
f1.close

#f2 = open('binelem.dat', 'r')
#data2 = [map(float, line.split()) for line in f2] 
#time = array([bit[0] for bit in data2]) 
#abin = array([bit[1] for bit in data2]) 
#ebin = array([bit[2] for bit in data2])
#ibin = array([bit[3] for bit in data2])
#lanbin = array([bit[4] for bit in data2])
#lopbin = array([bit[5] for bit in data2])
#meabin = array([bit[6] for bit in data2])
#f2.close

fig1 = figure()
subplot(511)
plot(time,acm)
xlabel('time'), ylabel('a')
subplot(512)
plot(time, ecm)
xlabel('time'), ylabel('e')
subplot(513)
plot(time, icm)
xlabel('time'), ylabel('inc')
subplot(514)
plot(time, lancm)
xlabel('time'), ylabel('q')
subplot(515)
plot(time, lopcm)
xlabel('time'), ylabel('d')
#subplot(616)
#plot(time, meacm)
#xlabel('time'), ylabel('M')

#fig2 = figure()
#subplot(611)
#plot(time,abin)
#xlabel('time'), ylabel('a')
#subplot(612)
#plot(time, ebin)
#xlabel('time'), ylabel('e')
#subplot(613)
#plot(time, ibin)
#xlabel('time'), ylabel('inc')
#subplot(614)
#plot(time, lanbin)
#xlabel('time'), ylabel('Omega')
#subplot(615)
#plot(time, lopbin)
#xlabel('time'), ylabel('omega')
#subplot(616)
#plot(time, meabin)
#xlabel('time'), ylabel('M')

show()
