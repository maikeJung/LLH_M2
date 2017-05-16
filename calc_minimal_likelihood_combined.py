'''
*Author: Maike Jung
*Date: 16.05.2017

*Purpose: create events for a certain SN neutrino spectrum with a certain mass and then determine the most pobable mass of these events via the method of maximum likelihood

SN - Model: Lawrence-Livermore

UNITS: mass2: eV2
       energy: MeV
       distance: Mpc
       time: s

noise of the detector: 1mHz
'''

import spectrum

from scipy.optimize import minimize_scalar
import numpy as np
import argparse
import os

def llh(massi2):
    # calculate LLH
    value = spectrum.getLLH( float(massi2), distance, events, useTriggerEff, useEnergyRes, noise, eventTime, eventEnergy, logTime, logTimeConv)
    return value

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
			help="Consider engery dependent trigger efficiency.")
    parser.add_argument("--perfect-reco", dest='energyRes', default=True,
			action='store_false', help="Take energy resolution of the detector into account.")
    parser.add_argument("--nfits", default=1, type=int,
			help="No. of pseudo-experiments to generate and fit.")
    args = parser.parse_args()
    
    #set parameters:
    mass2 = args.mass2; distance = args.distance; events = args.nevents; nfits = args.nfits
    useTriggerEff = args.triggEff; useEnergyRes = args.energyRes
    #RESE and REST need to be the same as in the header file!
    RESE = 600
    REST = 1000
    #expected noise rate 1mHz
    noise = pow(10,-3)
    
    # create array with logarithmic time bins
    pylogTime = np.logspace(-5.0,1.0,num=(REST+1))
    logTime = spectrum.doubleArray(REST+1)
    for i in range(len(pylogTime)):
        logTime[i] = pylogTime[i]

    # create a 2nd array with log times for the time convolution
    pylogTime2 = np.logspace(0.48,-5.0,num=(0.3*REST))
    pylogTime2 *= -1
    pylogTimeConv = np.concatenate((pylogTime2, [0], pylogTime))
    logTimeConv = spectrum.doubleArray(int(1.3*REST)+2)
    for i in range(len(pylogTimeConv)):
        logTimeConv[i] = pylogTimeConv[i]
        
    # create spectrum from which the events are drawn
    spectrumGen = spectrum.doubleArray( (RESE - 1) * REST )
    spectrum.createSpectrum(spectrumGen, mass2, distance, events, useEnergyRes, useTriggerEff, noise, logTime, logTimeConv);

    # arrays in which the pseudo experiments (their time and energy) will be stored
    eventEnergy = spectrum.intArray(int(events))
    eventTime = spectrum.intArray(int(events))

    # find maximum in the spectrum - needed to draw random events from spectrum
    maxSpectrum = spectrum.findSpectrumMax(spectrumGen)

    # set seed for random number generator - and store it
    spectrum.getSeed(distance, mass2, events, noise)

    if not os.path.exists('DATA'):
        os.makedirs('DATA')
    # performe likelihood calculation
    for i in range(nfits):
        # create one event
        spectrum.createEventsArray(events, spectrumGen, maxSpectrum, eventTime, eventEnergy,i, logTime)
        # find the mass for which the likelihood is minimal and store it
        x_min = minimize_scalar(llh, bounds=(-2.0,5.0), method='bounded', options={'disp':1,'xatol':0.005})
        print i, x_min.nfev, x_min.x
        with open("DATA/masses2_"+str(distance)+"Mpc_"+str(events)+"Events_"+str(mass2)+"eV2_.txt", "a") as myfile:
            myfile.write(str(x_min.x) + '\n')

    print 'DONE'
