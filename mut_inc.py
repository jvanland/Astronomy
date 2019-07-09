#!/usr/bin/python

from pylab import *
from matplotlib import pyplot as plt

label_size = 15
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

f1 = open('mut_inc.dat', 'r')
data1 = [map(float, line.split()) for line in f1]  
i_init = array([bit[0] for bit in data1]) 
delta_1 = array([bit[1] for bit in data1])
delta_2 = array([bit[2] for bit in data1])
delta_3 = array([bit[3] for bit in data1])
f1.close

error = 0.5

fig1, ax1 = plt.subplots()
ax1.set_xscale('log')
plt.xlabel(r'$\mathrm{Initial\, Outer\, Inclination (Deg)}$',fontsize=20)
plt.ylabel(r'$\Delta i_{\mathrm{tot}} \mathrm{(Deg)}$',fontsize=20)
ax1.set_ylim(0,8)
ax1.set_xlim(0.05,51.2)
ax1.errorbar(i_init,delta_1,yerr=error,fmt='r*',markersize=10)
ax1.errorbar(i_init,delta_2,yerr=error,fmt='bo',markersize=7)
ax1.errorbar(i_init,delta_3,yerr=error,fmt='g^',markersize=7)
ax1.set_xticks(i_init)
ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

plt.tight_layout()
plt.savefig("mut_inc.eps")
plt.show()
