/* File needed for SWIG to build the python module "spectrum" */

%module spectrum

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);

%{
#define SWIG_FILE_WITH_INIT
#include "spectrum.h"
%}


void createSpectrum(double *spectrum, double mass2, double distance, double events, bool useEnergyRes, bool useTriggerEff, double noise, double *logTime, double *logTimeConv);

double getLLH(double mass2, double distance, double events, bool triggEff, bool energyRes, double noise, int *eventTime, int *eventEnergy, double *logTime, double *logTimeConv);

void createEventsArray(double events, double *spectrum, double max, int *timeArray, int *energyArray, int filenumber, double *logTime);
void getSeed(double distance, double mass2, double events, double noise);
double findSpectrumMax(double *spectrum);

