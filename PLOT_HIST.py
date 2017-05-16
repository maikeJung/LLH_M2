import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.patches as mpatches

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
    bins = np.arange(input_mass-2.0 - bin_width/2.0, input_mass+2.0, bin_width)
    values_pseudo, m = np.histogram(masses, bins=bins)
    mass_hist_pseudo = np.arange(input_mass-2.0, input_mass+2.0 - bin_width/2.0, bin_width)
    return mass_hist_pseudo, values_pseudo/float(events)
    #return mass_hist_pseudo, values_pseudo


def calcError(mass, events, distance, bin_width):
    # calculate the 1 sigma uncertainty, by checking where 68% of events are detected
    mass_hist, values = readData("DATA/masses2_"+str(distance)+"Mpc_"+str(events)+".0Events_"+str(mass)+"eV2.txt", mass, bin_width)
    mass_hist = mass_hist[::-1]
    values = values[::-1]
    values = np.cumsum(values)
    f = interp1d(values, mass_hist)
    #return f(0.84135), f(0.5), f(0.15865)
    # 90%C.L., 2sigma, 1sigma and mean
    return f(0.97725),f(0.9), f(0.84135), f(0.5), f(0.15865), f(0.02275)
    
    
def plotHist(mass, events, distance, bin_width):
    mass_hist, values = readData("DATA/masses2_"+str(distance)+"Mpc_"+str(events)+".0Events_"+str(mass)+"eV2.txt", mass, bin_width)
    # create a histogram of pseudo experiments
    fig=plt.figure()
    ax = fig.add_subplot(111)
    plt.gcf().subplots_adjust(top=0.85) # to be able to read the lables
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.15)
    axes.yaxis.set_tick_params(labelsize=17, width=1)

    plt.bar(mass_hist-bin_width/2.0, values, bin_width, color="blue", label=str(distance)+"Mpc "+str(events)+"Events \n"+str(mass)+"eV", alpha=0.6)

    #plt.xlim(0.6,1.4)
    plt.xlabel('mass [eV]', fontsize=18)
    plt.ylabel('# of pseudo experiments', fontsize=18)
    plt.legend(loc='upper right', fontsize=14)
    #fig.savefig('PLOTS10000/'+str(mass)+'eV_'+str(distance)+'Mpc_'+str(events)+'Events.png')


def calculate_error(input_mass, events, distance):
    lower = []
    value = []
    upper = []
    lower2sigma = []
    upper2sigma = []
    low90CL = []
    for mass in input_mass:
        low2sigma, t90CL, low, val, upp, upp2sigma = calcError(mass, int(events), distance, 0.0001)
        lower.append(low)
        value.append(val)
        upper.append(upp)
        lower2sigma.append(low2sigma)
        upper2sigma.append(upp2sigma)
        low90CL.append(t90CL)
    lower = np.asarray(lower)
    value = np.asarray(value)
    upper = np.asarray(upper)
    lower2sigma = np.asarray(lower2sigma)
    upper2sigma = np.asarray(upper2sigma)
    low90CL = np.asarray(low90CL)
    return lower2sigma, low90CL, lower, value, upper, upper2sigma

#set parameters
distance = 1.0
events = 160
input_mass = [0.0, 0.01, 0.3, 0.5]
lower2sigma, test90CL, lower, value, upper, upper2sigma = calculate_error(input_mass, events, distance)

#plot
fig=plt.figure()
plt.gcf().subplots_adjust(top=0.85) # to be able to read the lables
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)

minput, = plt.plot(input_mass,input_mass,color='#000000', linewidth=1.5, linestyle='--', label=r'$m^2_{\mathrm{input}}=m^2_{\mathrm{found}}$')
mfound, = plt.plot(input_mass,value,color='red', linewidth=1.5, linestyle='--', label=r'$<m^2_{\mathrm{found}}>$')
m90CL, = plt.plot(input_mass,test90CL,color='blue', linewidth=1.5, linestyle='-', label=r'$90\%$ C.L.')

green_patch = mpatches.Patch(color='#66ff33', label=r'$1\,\sigma$')
yellow_patch = mpatches.Patch(color='#ffff00', label=r'$2\,\sigma$')
plt.fill_between(input_mass, lower2sigma, upper2sigma,edgecolor='#ffff00', facecolor='#ffff00')
plt.fill_between(input_mass, lower, upper,edgecolor='#66ff33', facecolor='#66ff33')

plt.title('Uncertainty of determined neutrino mass \n' + r'real detector with constant noise floor', fontsize=16)
plt.xlabel(r'$m^2_{\mathrm{input}}$ [eV$^2$]', fontsize=16)
plt.ylabel(r'reconstructed mass $m^2_{\mathrm{found}}$ [eV$^2$]', fontsize=16)

plt.legend(loc='upper left', numpoints = 1, fontsize=14.5, handles=[green_patch, yellow_patch, minput, mfound, m90CL])

axes = plt.gca()
axes.set_xlim([0.0,1.501])
axes.set_ylim([0.0,1.6])
axes.xaxis.set_tick_params(labelsize=15, width=1)
axes.yaxis.set_tick_params(labelsize=15, width=1)

#plot histogram for certain configuration
hit_width = 0.06
plotHist(0.3,events,distance, hit_width)
plotHist(0.5,events,distance, hit_width)

plt.show()
