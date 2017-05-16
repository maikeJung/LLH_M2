/* 
*Author: Maike Jung
*Date: 15.05.2017

*Purpose: Draw random events from a certain mass2 spectrum

SN - Model: Lawrence-Livermore
 
UNITS: mass2: eV
       energy: MeV
       distance: Mpc
       time: s

BINNING: logarithmic, defined in header-file (spectrum.h)
*/

#include "spectrum.h"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>


void createEventsArray(double events, double *spectrum, double max, int *timeArray, int *energyArray, int filenumber, double *logTime){
    /*create array with the times and energies of the events, by drawing an arbitrary time and energy and an arbitrary number between 0 and the maximum of the spectrum. the event is keept if the value of the spectrum at that time/energy exceeds randCheck*/
    int eventsGenerated = 0;
    int randE, randT;
    double randCheck;

    while(eventsGenerated < events){
        randE = rand() % (RESE-1);
        randT = rand() % (REST);
        randCheck = rand()*max/RAND_MAX;

        if (spectrum[randT*(RESE-1)+randE] >= randCheck){
            timeArray[eventsGenerated] = randT;
            energyArray[eventsGenerated] = randE;
            eventsGenerated ++;
        }
    }
}

void getSeed(double distance, double mass2, double events, double noise){
    unsigned seed = time(NULL);
    FILE *f = fopen("seeds.txt", "a");
    fprintf(f, "%f %f %.0f %f %d\n", distance, mass2, events, noise, seed);
    srand( seed );
    fclose(f);
}

double findSpectrumMax(double *spectrum){
    /* find the maximum value in the spectrum - needed to draw random events */
    double max = 0.0;
    int i;
    for(i=0; i<((RESE-1)*REST); i++){
        if(spectrum[i]>max) max = spectrum[i];
    }
    return max;
}
