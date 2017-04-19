/* 
*Author: Maike Jung
*Date: 26.01.2017

*Purpose: Draw random events from a certain mass2 spectrum

SN - Model: Lawrence-Livermore
    time spectrum is convoluted with the first hit distribution, to account for not knowing the absolute arrival times

UNITS: mass2: eV
       energy: MeV
       distance: Mpc
       time: s

BINNING: defined in header-file (spectrum.h)
*/

#include "spectrum.h"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>


int getArrayIndexE(double time, double *logTime){
    //determine index of the convoluted time array
    int index = 0;
    if (time > 9.5){
        index = REST;
    }
    if (time <= 0){
        index = 0;
    }
    int i;
    for (i=0;i<REST;i++){
        if(time < logTime[i]){
            index = i;
            break;
        }
    }
    return index;
}

void createEventsArray(double events, double *spectrum, double max, int *timeArray, int *energyArray, int filenumber, double *logTime){

    /*creates events and writes them in a time and energy array*/
    /* - to store Data 
    char filename[sizeof "DATA/Events/event_log_1111.txt"];
    sprintf(filename, "DATA/Events/event_log_%d.txt", filenumber);
    FILE *f = fopen(filename, "w");
    if (f == NULL){
        printf("Error opening file!\n");
        exit(1);
    }*/

    int eventsGenerated = 0;
    int randE, randT;
    double randTtest;
    double randCheck;
    while(eventsGenerated < events){
        randE = rand() % (RESE-1);
        randT = rand() % (REST);
        //randTtest = (double)rand()/(double)(RAND_MAX/6.0);
        //randT = getArrayIndexE(randTtest, logTime);
        
        //printf("%f %d\n", randTtest, randT);
        randCheck = rand()*max/RAND_MAX;
        
        //printf("events %d %d %f \n", randE, randT, max);
        if (spectrum[randT*(RESE-1)+randE] >= randCheck){
            //printf("found %f %d\n", randTtest, randT);
            timeArray[eventsGenerated] = randT;
            energyArray[eventsGenerated] = randE;
            //fprintf(f, "%d %d\n", randE, randT);
            eventsGenerated ++;
        }
    }
    //fclose(f);
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
        //printf("val %d %f\n", i, spectrum[i]);
        if(spectrum[i]>max) max = spectrum[i];
    }
    return max;
}
