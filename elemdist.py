#!/usr/bin/python

from pylab import *

f1 = open('elements.dat', 'r')
#f2 = open('elemfin3.dat', 'r')
data1 = [map(float, line.split()) for line in f1]  
a1 = array([bit[1] for bit in data1]) 
e1 = array([bit[2] for bit in data1])
#r1 = array([bit[5] for bit in data1])
f1.close
#data2 = [map(float, line.split()) for line in f2]  
#a2 = array([bit[1] for bit in data2]) 
#e2 = array([bit[2] for bit in data2])
#r2 = array([bit[5] for bit in data2])
#f2.close

abins = arange(0,10000,100)
ha1,abins = histogram(a1, bins=abins)
#ha2,abins = histogram(a2, bins=abins)
#aresidual = ha2-ha1
abin_centers = 0.5 * diff(abins) + abins[:-1]
ebins = arange(0.0,1.0,0.01)
he1,ebins = histogram(e1, bins=ebins)
#he2,ebins = histogram(e2, bins=ebins)
#eresidual = he2-he1
ebin_centers = 0.5 * diff(ebins) + ebins[:-1]
#rbins = arange(0,14000,20)
#hr1,rbins = histogram(r1, bins=rbins)
#hr2,rbins = histogram(r2, bins=rbins)
#rresidual = hr2-hr1
#rbin_centers = 0.5 * diff(rbins) + rbins[:-1]

fig1 = figure()
subplot(131)
suptitle('Initial Particle Distribution',fontsize=22)
bar(abin_centers/1000,ha1,0.1)
xlabel('a (10^3 AU)',fontsize=20), ylabel('# of stars',fontsize=20)
xticks(fontsize=20),yticks(fontsize=20)
subplot(132)
bar(ebin_centers,he1,0.01)
#title('Initial Eccentricity Distribution',fontsize=20)
xlabel('e',fontsize=20), ylabel('# of stars',fontsize=20)
xticks(fontsize=20),yticks(fontsize=20)
#subplot(133)
#bar(rbin_centers/1000,hr1,0.1)
#xlabel('r (10^3 AU)',fontsize=20), ylabel('# of stars',fontsize=20)
#xticks(fontsize=20),yticks(fontsize=20)
#fig1.tight_layout()
#fig2 = figure()
#subplot(121)
#suptitle('Final - Initial Element Distribution',fontsize=22)
#bar(abin_centers/1000,aresidual,0.1)
#title('Final - Initial Semi-major Axis Distribution',fontsize=20)
#xlabel('a (10^3 AU)',fontsize=20), ylabel('Relative # of stars',fontsize=20)
#xticks(fontsize=20),yticks(fontsize=20)
#subplot(122)
#bar(ebin_centers,eresidual,0.01)
#title('Final - Initial Eccentricity Distribution',fontsize=20)
#xlabel('e',fontsize=20), ylabel('Relative # of stars',fontsize=20)
#xticks(fontsize=20),yticks(fontsize=20)
#fig2.tight_layout()
#fig2 = figure()
#subplot(131)
#suptitle('Final Particle Distribution',fontsize=22)
#bar(abin_centers/1000,ha2,0.1)
#xlabel('a (10^3 AU)',fontsize=20), ylabel('# of stars',fontsize=20)
#xticks(fontsize=20),yticks(fontsize=20)
#subplot(132)
#bar(ebin_centers,he2,0.01)
#title('Final Eccentricity Distribution',fontsize=20)
#xlabel('e',fontsize=20), ylabel('# of stars',fontsize=20)
#xticks(fontsize=20),yticks(fontsize=20)
#subplot(133)
#bar(rbin_centers/1000,hr2,0.1)
#xlabel('r (10^3 AU)',fontsize=20), ylabel('# of stars',fontsize=20)
#xticks(fontsize=20),yticks(fontsize=20)
#fig2.tight_layout()

show()

