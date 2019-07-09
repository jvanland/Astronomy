#!/usr/bin/python

from pylab import *
import sys

filename = "inspiral2.dat"
f1 = open(filename, 'r')
data = array([map(float, line.split()) for line in f1])
subset1 = data[data[:,0] == 10]
subset2 = data[data[:,0] == 5]
subset3 = data[data[:,0] == 2]
subset4 = data[data[:,0] == 1]
mbh = array([bit[0] for bit in data]) 
bnum = array([bit[1] for bit in data]) 
vnum = array([bit[2] for bit in data])
tfin10 = (0.6/0.4)*array([bit[4] for bit in subset1])
tfin5 = array([bit[4] for bit in subset2])
tfin2 = array([bit[4] for bit in subset3])
tfin1 = array([bit[4] for bit in subset4])
afin10 = array([bit[9] for bit in subset1])
afin5 = array([bit[9] for bit in subset2])
afin2 = array([bit[9] for bit in subset3])
afin1 = array([bit[9] for bit in subset4])
f1.close

N_ins10 = len(tfin10)
cum_ins10 = array([float(i+1)/100 for i in range(N_ins10)])
N_ins5 = len(tfin5)
cum_ins5 = array([float(i+1)/100 for i in range(N_ins5)])
N_ins2 = len(tfin2)
cum_ins2 = array([float(i+1)/100 for i in range(N_ins2)])
cum_ins2[-1] = cum_ins2[-2]
N_ins1 = len(tfin1)
cum_ins1 = array([float(i+1)/100 for i in range(N_ins1)])
cum_ins1[-1] = cum_ins1[-2]

#figure()
#xlabel('Merger Time (yrs)'), ylabel('Merger a (AU)')
#p1 = plot(tfin10,afin10, 'bo',label="10k")
#p2 = plot(tfin5,afin5, 'go',label="5k")
#p3 = plot(tfin2,afin2, 'bo',label="2.5k")
#p4 = plot(tfin1,afin1, 'ro',label="1k")
#legend(loc=5)

#f1 = open("cins.dat",'r')
#data = array([map(float, line.split()) for line in f1])
#tfin2 = array([bit[3] for bit in data])
#N_ins2 = len(tfin2)
#cum_ins2 = array([float(i+1)/100 for i in range(N_ins2)])

figure()
xlabel(r'$\mathrm{Merger\, Time\, (Relaxation\, Times)}$',fontsize=16), ylabel(r'$\mathrm{Cumulative\, Mergers\, (Fraction)}$',fontsize=16)
p1 = plot(sorted(tfin10/175000),cum_ins10,'r--')
p2 = plot(sorted(tfin5/56000),cum_ins5,'m')
p3 = plot(sorted(tfin1/20000),cum_ins1,'b--')
p4 = plot(sorted(tfin2/30000),cum_ins2,'g')
text(2.2,0.32,r'$\mathrm{5k}$')
text(2.2,0.26,r'$\mathrm{10k}$')
text(2.2,0.125,r'$\mathrm{2.5k}$')
text(2.2,0.085,r'$\mathrm{1k}$')
#legend(loc=2)
tight_layout()
#savefig("cumulative_1.eps")

show()
