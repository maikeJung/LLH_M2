/* File needed for SWIG to build the python module "spectrum" */

%module spectrum

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);

%{
#define SWIG_FILE_WITH_INIT
#include "spectrum.h"
%}

#ifdef USE_SP
typedef float user_data_t;
#else
typedef double user_data_t;
#endif

void createSpectrum(double *spectrum, double mass2, double distance, double events, bool useEnergyRes, bool useTriggerEff, double noise, double noise_events, double *logTime, double *logTimeConv);

void getEvent(int *eventEnergy, int *eventTime, double mass2, double distance, double events, int filenumber, double noise);
double getLLH(double mass2, double distance, double events, bool triggEff, bool energyRes, double noise, int *eventTime, int *eventEnergy, double noise_events, double *logTime, double *logTimeConv);

double getLLHLogBins(double mass2, double distance, double events, bool triggEff, bool energyRes, double noise, int *eventTime, int *eventEnergy, double noise_events, double *logTime, double *logTimeConv, double *spectrumGen);

void createEventsArray(double events, double *spectrum, double max, int *timeArray, int *energyArray, int filenumber, double *logTime);
void getSeed(double distance, double mass2, double events, double noise);
double findSpectrumMax(double *spectrum);

void calcLLH(double mass, double distance, double events, bool triggEff, bool energyRes, int filenumber, double noise, double noise_events, double *logTime, double *logTimeConv, int *eventEnergy, int *eventTime);
