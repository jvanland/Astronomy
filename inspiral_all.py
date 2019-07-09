#!/usr/bin/python

from pylab import *
from scipy import interpolate
from scipy.misc import derivative
import glob
import sys
import os

f2 = open('merit.dat','w')

for num in ('1k', '2.5k', '5k', '10k'):
    name = 'kozai{}*'.format(num)
    if num is '1k':
        Mass = 1e3
    if num is '2.5k':
        Mass = 2.5e3
    if num is '5k':
        Mass = 5e3
    if num is '10k':
        Mass = 1e4

    for filename in glob.glob(name):
        print filename
        f1 = open(filename,'r')
        data = [map(float, line.split()) for line in f1]
        com_inc = array([bit[6] for bit in data])
        ecm = array([bit[5] for bit in data])
        acm = array([bit[4] for bit in data])
        mut_inc = array([bit[3] for bit in data])
        ebin = array([bit[2] for bit in data])
        abin = array([bit[1] for bit in data])
        time = array([bit[0] for bit in data])
        f1.close

        spl = interpolate.UnivariateSpline(time,com_inc)
        deriv_new = gradient(spl(time))
        h = time[2]-time[1]
        span = len(time)
        deriv = [0]*(span-1)
        t_inc = [0]*(span-1)
        t_k = [0]*(span-1)
        t_new = [0]*(span-1)
        t_try = [0]*(span-1)
        delta_inc = [0]*(span-1)
        merit = [0]*(span-1)
        for i in range(span-1):
            if abin[i] < 0.9:
                span = i
                break
            t_new[i] = time[i]/1000
            deriv[i] = (com_inc[i+1]-com_inc[i])/h
            t_inc[i] = 360.0/deriv[i]
            if acm[i] < 0:
                acm[i] = abs(acm[i])
                ecm[i] = 0.5
            t_k[i] = 1.4e6*(Mass/1e6)**(-1.0)*(acm[i]/20626.5)**(3)*(1.0-ecm[i]*ecm[i])**(1.5)
            t_try[i] = 360.0/deriv_new[i]
            if t_k[i] < 2.0*h:
                t_k[i] = 2.0*h
            init = i-int((0.5*t_k[i])/h)
            if init < 0:
                init = 0
            fin = i+int((0.5*t_k[i])/h)
            if fin > span-1:
                fin = span-1
            delta_inc[i] = abs(com_inc[fin]-com_inc[init])
            delta_inc[i] = (3.1415926/180.0)*delta_inc[i]
            merit[i] = 1.0/delta_inc[i]
            #print init, fin, delta_inc[i], time[i], t_new[i], t_k[i], h
        avg = 0
        for i in range(span-11, span-1):
            avg += merit[i]
        avg = avg/10
        print >> f2, filename, avg

        #figure()
        #plot(time,com_inc)
        
        #figure()
        #plot(t_new, merit)
        #ylim(-0.1)
        #subplot(241)
        #title('Inclination and spline')
        #plot(time/1000, com_inc)
        #plot(time/1000, spl(time))
        #ylabel('degrees')
        #subplot(242)
        #title('slope of adjacent points')
        #plot(t_new, deriv)
        #ylabel('degrees/yr')
        #subplot(246)
        #title('timescale from adjacent slope')
        #plot(t_new,t_inc)
        #ylabel('yrs/360 degrees')
        #xlabel('10^3 yrs')
        #subplot(243)
        #title('slope from spline')
        #plot(time/1000, deriv_new)
        #ylabel('degrees/yr')
        #subplot(247)
        #title('timescale from spline')
        #plot(t_new,t_try)
        #ylabel('yrs/360 degrees')
        #xlabel('10^3 yrs')
        #subplot(245)
        #title('mutual inclination')
        #plot(time/1000, mut_inc)
        #ylabel('degrees')
        #xlabel('10^3 yrs')
        #subplot(244)
        #title('kozai timescale')
        #plot(t_new, t_k)
        #ylabel('yrs/kozai cycle')
        #xlabel('10^3 yrs')
        #subplot(248)
        #plot(time/1000, (1-ebin*ebin))
        #yscale('log')
        #ylabel(r'$1-e_{1}^{2}$',fontsize=16)

        #show()

f2.close
    

        

