#!/usr/local/bin/python

from pylab import *
import sys

pi = 3.141592653589793

file1 = "test.dat"
f1 = open(file1, 'r')
data1 = [map(float, line.split()) for line in f1] 
dist = array([bit[0] for bit in data1])
mass = array([bit[1] for bit in data1])
f1.close

fig1 = figure()
plot(dist,mass,'o')

show()
