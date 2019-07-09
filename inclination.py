#!/usr/bin/env python

from pylab import *
import sys
import glob
from scipy import stats
from scipy.stats import norm

num = 1000
alldelta = [ ]
for filename in glob.glob('kozai10k*.dat'):
    print filename
    f1 = open(filename,'r')
    data = [map(float, line.split()) for line in f1]
    inc = array([bit[6] for bit in data])
    f1.close

    size = len(inc)-1
    deltainc = size*[0]
    for i in range(size):
        if isnan(inc[i+1]):
            break
        deltainc[i] = abs(inc[i+1]-inc[i])
    alldelta = concatenate([alldelta,deltainc])
alldelta = [x for x in alldelta if x > 0.0]
X = log(alldelta)
bins = logspace(-5,2,num)
normbins = linspace(min(X),max(X),num)
hall = histogram(alldelta,bins)[0]
hnormal = histogram(X,normbins)[0]
centers = (num-1)*[0]
normcenters = (num-1)*[0]
for i in range(num-1):
    centers[i] = (bins[i+1]+bins[i])/2.0
    normcenters[i] = (normbins[i+1]+normbins[i])/2.0
mu = mean(X)
scale = exp(mu)
variance = var(X)
sigma = sqrt(variance)
shape = sigma
mode = exp(mu-variance)
print mu
print scale
print shape
print mode

area = 0.0
normarea = 0.0
for i in range(num-1):
    area+=(bins[i+1]-bins[i])*hall[i]
    normarea+=(normbins[i+1]-normbins[i])*hnormal[i]
print area
print normarea
hnorm = hnormal/normarea
hscale = hall/area

logpdf = stats.lognorm.pdf(centers, shape, loc=0, scale=scale)
npdf = normpdf(normcenters,mu,sigma)

figure()
suptitle("10k Without Stars Influencing Interior Orbit")
xlabel(r'$\mathrm{Delta\, Inclination}\, (\mathrm{Degrees})$',fontsize=16), ylabel(r'$\mathrm{Normalized\, Frequency}$',fontsize=16)
plot(centers,hscale)
plot(centers,logpdf)
#plot(normcenters,npdf)
#plot(normcenters,hnorm)
gca().set_xscale("log")
gca().set_yscale("log")
tight_layout()
#savefig("nostars_inc_hist.eps")

show()
