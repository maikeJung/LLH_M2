import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
import numpy as np
import scipy
from scipy import stats
from scipy.stats import norm
from scipy.stats import chisquare
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
import random

### PSEUDO EXPERIMENTS ###

# read in data
def readData(filename, input_mass, bin_width):
    h = open(filename,'r')
    masses = []
    for line in h:
	    m = float(line)
	    masses.append(m)
    h.close()

    events = len(masses)

    # store values in histogram -> values_pseudo: values of bins; mass_hist_pseudo: center of the bins
    bins = np.arange(input_mass-1.5 - bin_width/2.0, input_mass+1.5, bin_width)
    values_pseudo, m = np.histogram(masses, bins=bins)
    mass_hist_pseudo = np.arange(input_mass-1.5, input_mass+1.5 - bin_width/2.0, bin_width)
    return mass_hist_pseudo, values_pseudo/float(events)
    #return mass_hist_pseudo, values_pseudo


def calcError(mass, events, distance, bin_width):
    # calculate the 1 sigma uncertainty, by checking where 68% of events are detected
    #mass_hist, values = readData("DATA_TEST/masses_"+str(distance)+"Mpc_"+str(events)+".0Events_"+str(mass)+"eV_test.txt", mass, bin_width)
    mass_hist, values = readData("DATA10000M2/masses_"+str(distance)+"Mpc_"+str(events)+".0Events_"+str(mass)+"eV_0.01noiseEvents_1e-06Noise_test.txt", mass, bin_width)
    mass_hist = mass_hist[::-1]
    values = values[::-1]
    values = np.cumsum(values)
    f = interp1d(values, mass_hist)
    return f(0.84135), f(0.5), f(0.15865)
    
    
def plotHist(mass, events, distance, bin_width):
    #mass_hist, values = readData("DATA/masses_"+str(distance)+"Mpc_"+str(events)+".0Events_"+str(mass)+"eV_test.txt", mass, bin_width)
    mass_hist, values = readData("DATA10000M2/masses_"+str(distance)+"Mpc_"+str(events)+".0Events_"+str(mass)+"eV_0.01noiseEvents_1e-06Noise_test.txt", mass, bin_width)
    # create a histogram of pseudo experiments
    fig=plt.figure()
    ax = fig.add_subplot(111)
    plt.gcf().subplots_adjust(top=0.85) # to be able to read the lables
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.15)
    axes = plt.gca()
    #axes.set_xlim([min(masses)-0.001, max(masses)+0.001])
    #plt.xlim(xmin=0-bin_width/2.0, xmax=1.0)
    axes.xaxis.set_tick_params(labelsize=17, width=1)
    axes.yaxis.set_tick_params(labelsize=17, width=1)

    plt.bar(mass_hist-bin_width/2.0, values, bin_width, color="blue", label=str(distance)+"Mpc "+str(events)+"Events \n"+str(mass)+"eV", alpha=0.6)
    #plt.bar(mass_hist_lim-bin_width/2.0, values_lim, bin_width, color="green", label=str(distance)+"Mpc "+str(events)+"Events \n"+str(mass)+"eV, less events", alpha=0.6)

    #plt.xlim(0.6,1.4)
    plt.xlabel('mass [eV]', fontsize=18)
    plt.ylabel('# of pseudo experiments', fontsize=18)
    plt.legend(loc='upper right', fontsize=14)
    #fig.savefig('PLOTS10000/'+str(mass)+'eV_'+str(distance)+'Mpc_'+str(events)+'Events.png')


def calculate_error(input_mass, events, distance):
    lower = []
    value = []
    upper = []
    for mass in input_mass:
        low, val, upp = calcError(mass, int(events), distance, 0.0001)
        upp = np.sqrt(upp)
        lower.append(low)
        value.append(val)
        upper.append(upp)
    lower = np.asarray(lower)
    value = np.asarray(value)
    upper = np.asarray(upper)
    return lower, value, upper

distance = 5.0
events = 10
input_mass = [0.0, 0.0001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
#input_mass = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.3, 1.5]
lower, value, upper = calculate_error(input_mass, events, distance)

fig=plt.figure()
plt.gcf().subplots_adjust(top=0.85) # to be able to read the lables
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)

plt.plot(-1,-1,color='#66ff33', linewidth=1.5, linestyle='-', label=r'$D = '+str(distance)+'\,$Mpc,$\  '+str(events)+'\,$events')

plt.plot(input_mass,input_mass,color='#000000', linewidth=1.5, linestyle='--', label=r'$m_{input}=m_{found}$')
plt.plot(input_mass,value,color='red', linewidth=1.5, linestyle='--', label=r'$m_{found}$')

plt.fill_between(input_mass, lower, upper,edgecolor='#66ff33', facecolor='#66ff33')

plt.title('Uncertainty of determined neutrino mass \n' + r'real detector with constant noise floor', fontsize=16)
plt.xlabel('input mass [eV]', fontsize=16)
plt.ylabel('valid mass range [eV]', fontsize=16)
plt.legend(loc='lower right', numpoints = 1, fontsize=14.5)
axes = plt.gca()
axes.set_xlim([0.0,1.501])
axes.set_ylim([0.0,1.6])
axes.xaxis.set_tick_params(labelsize=15, width=1)
axes.yaxis.set_tick_params(labelsize=15, width=1)

hit_width = 0.06
plotHist(0.0,events,distance, hit_width)
plotHist(0.01,events,distance, hit_width)
plotHist(0.05,events,distance, hit_width)
plotHist(0.1,events,distance, hit_width)
#plotHist(0.2,events,distance, hit_width)
#plotHist(0.3,events,distance, hit_width)
plotHist(0.5,events,distance, hit_width)
#plotHist(0.7,events,distance, hit_width)
#plotHist(1.5,events,distance, hit_width)

plt.show()

# test hist
'''
low_05 = [ 0.00277997,  0.03988879,  0.09033385, 0.26893396,  0.39260267,  0.62407705,
  0.94041025,  1.24416853,  1.4485688 ] 
val_05 = [ 0.15095,     0.2277,      0.3067,      0.4157,      0.4898,      0.70860625,
  1.00035,     1.296325,    1.50423333] 
up_05 = [ 0.40624465,  0.40136453,  0.43778404,  0.51130755,  0.58991315,  0.77278358,
  1.06139487,  1.35848879,  1.5508092 ]

low_10 = [ 0.00293914,  0.03787748,  0.0764421,   0.23821465,  0.37437043,  0.60293186,
  0.91583475,  1.22101918,  1.41838214] 
val_10 = [ 0.1658,    0.2433,    0.32155,   0.41405,  0.5053,    0.711425,  1.0032,
  1.297075,  1.5001  ] 
up_10 = [ 0.4124289,   0.43456126,  0.48328985,  0.5534317,   0.62465695,  0.80247333,
  1.08471305,  1.37718081,  1.58351786]

low_50 = [ 0.00326077,  0.01175886,  0.07911387,  0.112282,    0.22555814,  0.50553169,
  0.7941452,  1.07479944,  1.26424604] 
val_50 = [ 0.1823,   0.2316,   0.308,    0.4049,   0.4892,   0.69835,  0.9964,   1.2852,
  1.5069 ]
up_50 = [ 0.49787149,  0.52133197,  0.5680635,   0.60617088,  0.67793548,  0.87286724,
  1.1894096,   1.51775582,  1.75745132]

input_mass = [0.01, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.3, 1.5]

fig=plt.figure()
plt.gcf().subplots_adjust(top=0.85) # to be able to read the lables
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)

plt.plot(-1,-1,color='#66ff33', linewidth=1.5, linestyle='-', label=r'$D = 0.5\,$Mpc,$\  650\,$events')
plt.plot(-1,-1,color='#ff0000', linewidth=1.5, linestyle='-', label=r'$D = 1.0\,$Mpc,$\  160\,$events')
plt.plot(-1,-1,color='#3300ff', linewidth=1.5, linestyle='-', label=r'$D = 5.0\,$Mpc,$\  10\,$events')

plt.plot(input_mass,input_mass,color='#000000', linewidth=1.5, linestyle='--', label=r'$m_{input}=m_{found}$')
#plt.plot(input_mass,value,color='red', linewidth=1.5, linestyle='--', label=r'$m_{found}$')

plt.fill_between(input_mass, low_50, up_50,edgecolor='#3300ff', facecolor='#3300ff')
plt.fill_between(input_mass, low_10, up_10,edgecolor='#ff0000', facecolor='#ff0000')
plt.fill_between(input_mass, low_05, up_05,edgecolor='#66ff33', facecolor='#66ff33')

plt.title('Uncertainty of determined neutrino mass \n' + r'real detector with constant noise floor', fontsize=16)
plt.xlabel('input mass [eV]', fontsize=16)
plt.ylabel('valid mass range [eV]', fontsize=16)
plt.legend(loc='lower right', numpoints = 1, fontsize=14.5)
axes = plt.gca()
axes.set_xlim([0.19,1.501])
axes.set_ylim([0.19,1.6])
axes.xaxis.set_tick_params(labelsize=15, width=1)
axes.yaxis.set_tick_params(labelsize=15, width=1)
plt.show()
'''
