#!/usr/bin/python

from pylab import *

cols = 24493
rows = 20
a_array = [[0] * cols for i in range(rows)]
e_array = [[0] * cols for i in range(rows)]
for i in range(rows):
    filename = "super%d.dat" % (i+9980)
    if (i != 15):
        f1 = open(filename, 'r')
        data = [map(float, line.split()) for line in f1]  
        a_array[i] = [bit[2] for bit in data]
        e_array[i] = [bit[3] for bit in data]
        f1.close
a_avg = [0]*cols
e_avg = [0]*cols
for i in range(cols):
    for j in range(rows):
        a_avg[i] += a_array[j][i]
        e_avg[i] += e_array[j][i]
    a_avg[i] = a_avg[i]/rows
    e_avg[i] = e_avg[i]/rows
a_avg = array(a_avg)
time = array(94.872*arange(0,cols,1.0))
subplot(121)
plot(time/1e6,a_avg/1000)
xlabel('Time (Million yrs)'), ylabel('Average Semi-major Axis (10^3 AU)')#, ylim([0,max(a_avg)/1000])
subplot(122)
plot(time/1e6,e_avg)
xlabel('Time (Million yrs)'), ylabel('Average Eccentricity')#, ylim([0,1])
show()
