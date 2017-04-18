import spectrum
from scipy.optimize import minimize
from scipy.optimize import brent
from scipy.optimize import minimize_scalar
from scipy.optimize import fmin
from scipy.interpolate import interp1d
import numpy as np
from datetime import datetime
import argparse
import os

def llh(massi2):
    # calculate LLH
    test = spectrum.getLLH( float(massi2), distance, events, useTriggerEff, useEnergyRes, noise, eventTime, eventEnergy, noise, logTime, logTimeConv)
    #print massi2, test
    return test

def createHist(masses, bin_width):
    # store determined masses in a histogram from 0 to 5 eV
    events = len(masses)
    # store values in histogram -> values_pseudo: values of bins; mass_hist_pseudo: center of the bins
    bins = np.arange(0.0 - bin_width/2.0, 5.0, bin_width)
    values_pseudo, m = np.histogram(masses, bins=bins)
    mass_hist_pseudo = np.arange(0.0, 5.0 - bin_width/2.0, bin_width)
    return mass_hist_pseudo, values_pseudo/float(events)

def calcError(masses):
    # calculate the 1 sigma uncertainty, by checking where 68% of events are detected
    mass_hist, values = createHist(masses, 0.001)
    mass_hist = mass_hist[::-1]
    values = values[::-1]
    values = np.cumsum(values)
    f = interp1d(values, mass_hist)
    return f(0.84135), f(0.5), f(0.15865)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
		description="Simulate SN Spectrum and fit pseudo-experiments.")
    parser.add_argument("-M2", "--mass2", default=1.0, type=float,
			help="Neutrino mass squared (in eV2)")
    parser.add_argument("-D", "--distance", default=5.0, type=float,
			help="SN Distance (in Mpc)")
    parser.add_argument("-N", "--nevents", default=10.0, type=float,
			help="Number of expected events from the SN")
    parser.add_argument("--perfect-trigger", dest='triggEff',
			default=True, action='store_false',
			help="Assume fully  eff. trigger across all energies.")
    parser.add_argument("--perfect-reco", dest='energyRes', default=True,
			action='store_false', help="Assume perfect energy reco.")
    parser.add_argument("--nfits", default=1, type=int,
			help="No. of pseudo-experiments to generate and fit.")
    parser.add_argument("--noiseb", default=0.01, type=float,
			help="Noise - exponential function")
    args = parser.parse_args()
    
    #TODO also import? also parse noise?
    RESE = 600
    REST = 1000
    #noise = pow(10,-3)*(10.0/REST)
    noise = pow(10,-3)*(10.0)
    noise_events = 0.01

    mass2 = args.mass2; distance = args.distance; events = args.nevents; nfits = args.nfits
    useTriggerEff = args.triggEff; useEnergyRes = args.energyRes
    #noise = args.noiseb

    #draw events and noise events from Poisson distribution
    #neutrinoEvents = np.random.poisson(events)
    #noiseEvents = np.random.poisson(noise_events)
    #events = neutrinoEvents + noiseEvents
    
    # create array with log times
    pylogTime = np.logspace(-5.0,1.0,num=(REST+1))
    logTime = spectrum.doubleArray(REST+1)
    for i in range(len(pylogTime)):
        logTime[i] = pylogTime[i]
    #print logTime[0], logTime[1], logTime[99], logTime[100]

    # create a 2nd array with log times for the time convolution
    pylogTime2 = np.logspace(0.48,-5.0,num=(0.3*REST))
    pylogTime2 *= -1
    pylogTimeConv = np.concatenate((pylogTime2, [0], pylogTime))
    logTimeConv = spectrum.doubleArray(int(1.3*REST)+2)
    for i in range(len(pylogTimeConv)):
        logTimeConv[i] = pylogTimeConv[i]
    #print logTimeConv[45]

    # create spectrum from which the events are drawn
    # TODO: draw number of events from Poisson distribution - probably needs to be created for every event
    spectrumGen = spectrum.doubleArray( (RESE - 1) * REST )
    spectrum.createSpectrum(spectrumGen, mass2, distance, events, useEnergyRes, useTriggerEff, noise, noise_events, logTime, logTimeConv);
    #for t in range(0,REST):
    #    for e in range(1,RESE):
    #        print t,e, spectrumGen[t*(RESE-1) +e-1]

    #exit();
    # arrays in which the pseudo experiments (their time and energy) will be stored
    eventEnergy = spectrum.intArray(int(events))
    eventTime = spectrum.intArray(int(events))

    # find maximum in the spectrum - needed to draw random events from spectrum
    maxSpectrum = spectrum.findSpectrumMax(spectrumGen)

    # set seed for random number generator - and store it
    spectrum.getSeed(distance, mass2, events, noise)

    if not os.path.exists('DATA'):
        os.makedirs('DATA')
    masses = []
    for i in range(nfits):
        # create one event
        spectrum.createEventsArray(events, spectrumGen, maxSpectrum, eventTime, eventEnergy,i, logTime)
        # find the mass for which the likelihood is minimal and store it
        x_min = minimize_scalar(llh, bounds=(-2.0,5.0), method='bounded', options={'disp':1,'xatol':0.005})
        print i, x_min.nfev, x_min.x
        masses.append(x_min.x)
        with open("DATA/masses_"+str(distance)+"Mpc_"+str(events)+"Events_"+str(mass2)+"eV_"+str(noise_events)+"noiseEvents_"+str(noise)+"Noise_test.txt", "a") as myfile:
            myfile.write(str(x_min.x) + '\n')

    # calculate 1 sigma uncertainty and store
    lower, value, upper = calcError(masses)

    with open("DATA/detection_error_test.txt", "a") as myfile:
        myfile.write(str(distance) +" "+ str(events) +" "+ str(mass2) +" " +str(noise)+ " " + str(lower) +" "+ str(value)+" " + str(upper) + '\n')

    print 'DONE'
