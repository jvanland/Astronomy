#!/usr/local/bin/python

from pylab import *
import struct

f1 = open('OUT3', 'rb')
data = f1.read()
t = struct.unpack('14f',data[0:56])
print t
f1.close()
