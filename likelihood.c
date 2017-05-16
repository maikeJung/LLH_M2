/* 
*Author: Maike Jung
*Date: 15.05.2017

*Purpose: Calculate the likelihood for the random events generated with generateEvents2.c to belong to a certain mass2 spectrum

SN - Model: Lawrence-Livermore
    
UNITS: mass2: eV
       energy: MeV
       distance: Mpc
       time: s

BINNING: logarithmic, see spectrum.h
*/

#include "spectrum.h"


double getLLH(double mass2, double distance, double events, bool triggEff, bool energyRes, double noise, int *eventTime, int *eventEnergy, double *logTime, double *logTimeConv){
    /*function that calculates the log likelihood - and will then be minimized in python*/
    double llh = 0.0;
    int i;
    double spectrum[(RESE-1)*REST];

    createSpectrum(spectrum, mass2, distance, events, energyRes, triggEff, noise, logTime, logTimeConv);
    for (i = 0; i < events; i++){
        if (spectrum[eventTime[i]*(RESE-1)+eventEnergy[i]] < pow(10,-200)){
            // check that the values are not too small, this might indicate an error in the calculation
            llh += -10000000;   
            printf("event number %d e %d t %d\n",i, eventEnergy[i], eventTime[i]);
            printf("value of spectrum very small - check \n");
        }
        else llh += log(spectrum[eventTime[i]*(RESE-1)+eventEnergy[i]]);
    }
    llh*=-1;
    //printf("mass llh %f %f \n", mass2, llh);
    return llh;
}
