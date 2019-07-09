#!/usr/local/bin/python

from pylab import *

f1 = open('escapecounter1.dat', 'r')
data1 = [map(float, line.split()) for line in f1] 
time1 = 0.031624*array([bit[0] for bit in data1]) 
f1.close
f2 = open('escapecounter2.dat', 'r')
data2 = [map(float, line.split()) for line in f2] 
time2 = 0.031624*array([bit[0] for bit in data2]) 
f2.close
f3 = open('escapecounter3.dat', 'r')
data3 = [map(float, line.split()) for line in f3] 
time3 = 0.031624*array([bit[0] for bit in data3]) 
f3.close
f4 = open('escapecounter4.dat', 'r')
data4 = [map(float, line.split()) for line in f4] 
time4 = 0.031624*array([bit[0] for bit in data4]) 
f4.close
f5 = open('escapecounter5.dat', 'r')
data5 = [map(float, line.split()) for line in f5] 
time5 = 0.031624*array([bit[0] for bit in data5]) 
f5.close
f6 = open('escapecounter6.dat', 'r')
data6 = [map(float, line.split()) for line in f6] 
time6 = 0.013176*array([bit[0] for bit in data6]) 
f6.close
f7 = open('escapecounter7.dat', 'r')
data7 = [map(float, line.split()) for line in f7] 
time7 = 0.013176*array([bit[0] for bit in data7]) 
f7.close

num1 = arange(0,6157,1)
num2 = arange(0,401,1)
num3 = arange(0,721,1)
num4 = arange(0,311,1)
num5 = arange(0,66,1)
num6 = arange(0,7868,1)
num7 = arange(0,134,1)

figure()
suptitle('Ejections for different configurations')
plot(time1,num1,label = "pmin17,pmax25")
plot(time2,num2,label = "none1")
#plot(time3,num3,label = "none2")
plot(time4,num4,label = "noneStart")
plot(time5,num5,label = "rkNoGrav")
plot(time6,num6,label = "pmin25,tmatch")
plot(time7,num7,label = "pmin25,tmatch,replace")
legend(loc=1)
show()
