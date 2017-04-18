/* 
*Author: Maike Jung
*Date: 15.11.2016

*Purpose: Calculate the likelihood for the random events generated with generateEvents.c to belong to a certain mass2 spectrum

SN - Model: Lawrence-Livermore
    time spectrum is convoluted with the first hit distribution, to account for not knowing the absolute arrival times

UNITS: mass2: eV
       energy: MeV
       distance: Mpc
       time: s

BINNING: see spectrum.h

*/

#include "spectrum.h"

/// load event
void getEvent(int *eventEnergy, int *eventTime, double mass2, double distance, double events, int filenumber, double noise){
    /*load events & store energy and time in arrays*/
    char filename[sizeof "DATA_COMP/10.00Mpc_700Events_1.57eV_event_1.45eV_10.5Mpc_1000Events_real_1999111.txt"];
    sprintf(filename, "DATA_COMP/%.2fMpc_%.0fEvents_%.5feV/events_%.2feV_%.2fMpc_%.0fEvents_real_%d.txt",distance, events, mass2, mass2, distance, events, filenumber);

    //char filename[sizeof "DATA/NOISETEST/10.00Mpc_700Events_1.57eV_event_b000_1.45eV_10.5Mpc_1000Events_real_1111.txt"];
    //sprintf(filename, "DATA/NOISETEST/%.2fMpc_%.0fEvents_%.2feV_b%.3f/events_%.2feV_%.2fMpc_%.0fEvents_real_%d.txt",distance, events, mass2, noise, mass2, distance, events, filenumber);

    FILE *f = fopen(filename, "r");
    int i;
    for(i = 0; i < events; eventEnergy[i++] = 1);
    for(i = 0; i < events; eventTime[i++] = 1);
    for (i = 0; i < events; i++){
        fscanf(f, "%d %d", &eventEnergy[i], &eventTime[i]);
    }
}
///function that calculates likelihood - and will then be minimized in python
double getLLH(double mass2, double distance, double events, bool triggEff, bool energyRes, double noise, int *eventTime, int *eventEnergy, double noise_events, double *logTime, double *logTimeConv){

    double llh = 0.0;
    int i;
    user_data_t spectrum[(RESE-1)*REST];

	//double *spectrum= (double*) malloc((RESE-1) * REST * sizeof(double));
    createSpectrum(spectrum, mass2, distance, events, energyRes, triggEff, noise, noise_events, logTime, logTimeConv);
    for (i = 0; i < events; i++){
        if (spectrum[eventTime[i]*(RESE-1)+eventEnergy[i]] < pow(10,-200)){
            llh += -10000000;   
            printf("event number %d e %d t %d\n",i, eventEnergy[i], eventTime[i]);
            printf("value of spectrum very small - check \n");
        }
        else llh += log(spectrum[eventTime[i]*(RESE-1)+eventEnergy[i]]);
    }
    llh*=-1;
    return llh;
}




int main(void){
	/*set parameters*/
    /*flag for trigger efficiency*/
    bool triggEff = true;
    bool energyRes = true;
    double mass2 = 1.0;
    double distance = 1.0;
    double events = 160;  
    double noise_events = 0.01;  
    int filenumber;
    double noise = pow(10,-5);

    /*calculate uncertainty for certain configuration*/
    for (filenumber=1; filenumber<1; filenumber++){ 
        printf("evaluating file %d \n", filenumber);
        calcLLH(mass2, distance, events, triggEff, energyRes, filenumber, noise, noise_events);
    }

    printf("DONE\n");
    return 0;
}
