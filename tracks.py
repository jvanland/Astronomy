#!/usr/bin/python

from pylab import *
import glob
import sys

def esc(array,min):
    for index, item in enumerate(array):
        if item < min:
            return index

filename = "inspiral2.dat"
f1 = open(filename, 'r')
data = array([map(float, line.split()) for line in f1])
subset1 = data[data[:,0] == 10]
subset2 = data[data[:,0] == 5]
subset3 = data[data[:,0] == 2]
subset4 = data[data[:,0] == 1]
tfin10 = (0.6/0.4)*array([bit[4] for bit in subset1])/1000.0
tfin5 = array([bit[4] for bit in subset2])/1000.0
tfin2 = array([bit[4] for bit in subset3])/1000.0
tfin1 = array([bit[4] for bit in subset4])/1000.0
afin10 = array([bit[9] for bit in subset1])
afin5 = array([bit[9] for bit in subset2])
afin2 = array([bit[9] for bit in subset3])
afin1 = array([bit[9] for bit in subset4])
f1.close

figure()
xlabel(r'$\mathrm{Time} (10^{3}\mathrm{yrs})$',fontsize=16), ylabel(r'$a_2(\mathrm{AU})$',fontsize=16)
for filename in glob.glob('elem*.dat'):
    f1 = open(filename,'r')
    data = [map(float, line.split()) for line in f1]
    sma = array([bit[1] for bit in data])
    ecc = array([bit[2] for bit in data])
    time = array([bit[0] for bit in data])/1000.0
    dist = array([bit[5] for bit in data])
    f1.close
    if time[-1] < 35000.0:
        Rmin = 5.348
    elif time[-1] < 130000.0:
        Rmin = 9.097
    elif time[-1] < 360000.0:
        Rmin = 11.454
    dex = argmax(dist < Rmin)
    #if dex > 0:
        #plot(time[dex],sma[dex],'*',ms = 30)
    n = len(sma)
    index = [ ]
    for i in range(n):
        if sma[i] < 0.0:
            index.append(i)
    sma_2 = delete(sma,index)
    time_2 = delete(time,index)
    n = len(sma_2)
    sma_3 = [ ]
    time_3 = [ ]
    for i in range(n):
        if i > 0 and i < n-1 and sma_2[i] < sma_2[i-1] and sma_2[i] < sma_2[i+1]:
            sma_3.append(sma_2[i])
            time_3.append(time_2[i])
    n = len(sma_3)
    sma_4 = [ ]
    time_4 = [ ]
    for i in range(n):
        if i > 0 and i < n-1 and sma_3[i] < sma_3[i-1] and sma_3[i] < sma_3[i+1]:
            sma_4.append(sma_3[i])
            time_4.append(time_3[i])
        
            #plot(time_2,sma_2)
    plot(time_4,sma_4,color='0.75')
    #plot(time,sma*(1-ecc),color='0.75')

    #p4 = semilogy(tfin10,afin10, 'bo',label="10k")
    #for i in range(10):
    #name = "create%d.dat" % (i+1)
    #f1 = open(name,'r')
    #data = [map(float, line.split()) for line in f1]
    #c1 = array([bit[0] for bit in data])
    #c2 = array([bit[1] for bit in data])
    #plot(c2,c1,color='0.75')
    #f1.close

    #f1 = open("cins.dat",'r')
    #data = [map(float, line.split()) for line in f1]
    #afin2 = array([bit[2] for bit in data])
    #tfin2 = array([bit[3] for bit in data])
    #f1.close

    #p4 = loglog(tfin10,afin10, 'bo',label="10k")
p1 = plot(tfin10,afin10, 'r*',markersize=10)
p2 = plot(tfin5,afin5, 'ms',markersize=5)
p3 = plot(tfin1,afin1, 'bo',markersize=7)
p4 = plot(tfin2,afin2, 'g^',markersize=7)
tight_layout()
savefig("tracks.eps")
show()
