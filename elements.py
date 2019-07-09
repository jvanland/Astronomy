#!/usr/bin/python

from pylab import *
import sys

if len(sys.argv) != 2:
    print "usage: ./elements.py infile"
    sys.exit(2)    

filename = sys.argv[1]
f1 = open(filename, 'r')
data = [map(float, line.split()) for line in f1] 
time = array([bit[0] for bit in data]) 
abin = array([bit[1] for bit in data]) 
ebin = array([bit[2] for bit in data])
pbin = array([bit[3] for bit in data])
acm = array([bit[4] for bit in data])
ecm = array([bit[5] for bit in data])
inc = array([bit[6] for bit in data])
f1.close

inc_max = [141]*len(time)
inc_min = [39]*len(time)
inc_good = [90]*len(time)
abin_zero = [0]*len(time)
close = min(pbin)
emax = max(ebin)
ctime = 94.872*pbin.argmin()
clabel = "closest approach %f AU" % (close)
#print emax

figure()
subplot(411)
plot(time/1000,abin,label="abin (AU)",linewidth=2)
plot(time/1000, abin_zero, 'k--')
ylabel(r'$a_1(\mathrm{AU})$',fontsize=16)
ylim(-0.5,1.2)
tick_params(axis='x',labelbottom='off')
subplot(414)
plot(time/1000, pbin,linewidth=2)
ylabel(r'$q_1(\mathrm{AU})$',fontsize=16)
ylim(0,1.2)
xlabel(r'$\mathrm{Time} (10^{3}\mathrm{yrs})$',fontsize=16)
subplot(413)
plot(time/1000, (1-ebin*ebin),label="ebin",linewidth=2)
yscale('log')
ylabel(r'$1-e_{1}^{2}$',fontsize=16)
ylim(0,1)
#xlabel('Time ($10^{3}$yrs)')
tick_params(axis='x',labelbottom='off')
#ylim(0,1)
subplot(412)
plot(time/1000,inc,label="inc (rad)",linewidth=2)
ylabel(r'$i_{tot}(\mathrm{Deg})$',fontsize=16)
tick_params(axis='x',labelbottom='off')
tight_layout()
plot(time/1000, inc_max, 'k--')
plot(time/1000, inc_min, 'k--')
plot(time/1000, inc_good, 'k--')
#xticks([0,0.4,0.8,1.2,1.6, 2.0, 2.4],fontsize=20),yticks(fontsize=20)
#ylim(0,5)
#xlim(0,160000)
#legend(loc=6,prop={'size':20})
tight_layout()
savefig("exampleSep.eps")

#subplot(122)
#plot(time/100000,abin,label="abin (AU)",linewidth=2)
#plot(time/100000, pbin,label="pericenter (AU)",linewidth=2)
#plot(time/100000, ebin,label="ebin",linewidth=2)
#ylim(0,1.1), xlim(15.2,15.4)
#xticks(fontsize=20),yticks(fontsize=20),xlabel('time (10^5 yrs)',fontsize=20)

figure()
plot(time/1000,acm*(1-ecm),label="c.o.m. pericenter")
plot(time/1000,acm,label="c.o.m. a")
#xlabel('time (yrs)')
#title('Binary 9991')
#legend(loc=3)
#ylim(0,1000)
#xlim(0,160000)

#f2 = open('superelem9995.dat', 'r')
#data2 = [map(float, line.split()) for line in f2]
#acm2 = array([bit[1] for bit in data2])
#pcm = array([bit[4] for bit in  data2])
#f2.close

#figure()
#plot(time,acm2)
#plot(time,pcm)
#xlim(0,10000)

show()
