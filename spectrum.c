/*
*Author: Maike Jung, Thomas Ehrhardt
*Date: 24.01.2017

*Purpose: create the arrival time spectrum of the neutrinos, that can then be uses to 
    generate random events: generateEvents.c
    calculate the likelihood for these events:  likelihood.c
    calculate the binned-likelihood using an Asimov dataset:    binned_likelihood.c

SN - Model: Lawrence-Livermore
    time spectrum is convoluted with the first hit distribution, to account for not knowing the absolute arrival times

UNITS: mass2: eV
       energy: MeV
       distance: Mpc
       time: s

noise of the detector: 1mHz
*/

#include "spectrum.h"

double getDeltaT(double E, double mass2, double dist){
    /* time shift due to neutrino mass2 - factor of 51.4 to get the proper units*/
    double tDelta = dist*51.4635*(mass2/(E*E));
    return tDelta;
}

double getTimeDelay(double t, double E, double mass2, double dist){
    return t - getDeltaT(E, mass2, dist);
}

double LL_time_spectrum_shifted(double t, double E, double mass2, double dist){
    /* retrun the value of the LL time spectrum at the shifted time*/
    double time = getTimeDelay(t, E, mass2, dist);
    if (time <= 0){
        /* spectrum not defined below 0 - there are no events before time 0 and after 10s*/
        return 0.0;
    }
    if (time >= TMAX-0.001){
        time = TMAX-0.001;
    }
    return LL_time_spectrum(time);
}

double LLSpectrumTotal(double t, double E, double mass2, double dist){
    /* 2D arrival time/energy probability for a certain mass2/distance - normalized */
    return LL_time_spectrum_shifted(t, E, mass2, dist)*LL_energy_spectrum(E);
}

void cumSumT(double *arrayToSum, double *cumulative){
    /*calculate the cumulative sum of the time array*/
    int k, l;
    double cum;
    for (k = 0; k < REST; k++){
        cum = 0.0;
        for (l = 0; l <= k; l++){
            cum += arrayToSum[l];
        }
        cumulative[k] = cum;
    }
}

void firstHitDistWeightedArrivalTimeDist(double *arrivalTimeDist, double *cumulative, double events, double *result){
    /*calculate the 1st Hit Distribution - and normalize it*/
    int m;
    double count = 0.0;
    for (m = 0; m < REST; m++){
        result[m] = arrivalTimeDist[m]*events*pow((1 - cumulative[m]), events-1);
        count += result[m];
    }
    for (m = 0; m < REST; m++){
        result[m] = result[m]/count;
    }
}

/* calculate the probability to get the first hit after a certain amount of time */
void ProbFirstHitDist (double mass2, double dist, double events, double *result, double *logTime){
    /*arrival time distribution of all the hits (for a certain mass2) - project the E,t spectrum
    on the t axis*/

    /*this array needs to be in log space*/
    double totalArrivalTimeDist[REST];
    int i;
    double sum;
    double y, e;
    /*go over i in log coordinates or choose proper entry from log array*/
    for (i = 0; i < REST; i++){
        /* set the sum to zero for each time bin */
        sum = 0.0;
        /*Integrate over the energy part for every time bin. We move in 0.01 MeV
        steps up to 60 MeV. For each pair of time and energy, we compute the
        product of time and energy PDF ("LLSpectrumTotal"), continually 
        incrementing the sum*/
        for (e = 0.01; e < EMAX; e += 0.01) {
            y = LLSpectrumTotal((logTime[i]+logTime[i+1])/2.0, e, mass2, dist);
            //y = LLSpectrumTotal(i*STEPT, e, mass2, dist);
            sum += y * 0.01;
        }
        totalArrivalTimeDist[i] = sum*(logTime[i+1]-logTime[i]);
        //totalArrivalTimeDist[i] = sum*STEPT;
    }

    double cumulative[REST];
    cumSumT(totalArrivalTimeDist, cumulative);

    firstHitDistWeightedArrivalTimeDist(totalArrivalTimeDist, cumulative, events, result);
}

void convolveHitDistWithLLTimeSpec(double *hitDist, double *convolSpec){
    int i, j;
    double pNew;
    /*perform the convolution*/
    for (i = 0; i < REST*1.3; i++){
        pNew = 0.0;
        for (j = 0; j < REST; j++){
            if ((i-0.3*REST + j) < REST && (i-0.3*REST + j) > 0){
                pNew += hitDist[j] * LL_time_spectrum( (j+i-0.3*REST)*STEPT );
            }
        }
        convolSpec[i] = pNew;
    }
}

void convolveHitDistWithLLTimeSpecLog(double *hitDist, double *convolSpec, double *logTime, double *logTimeConv){
    int i, j;
    double t, a, da;
    double pNew;
    /*perform the convolution on log array*/
    for (i = 0; i < REST*1.3; i++){
        t = (logTimeConv[i]+logTimeConv[i+1])/2.0; // calculate the proper time
        pNew = 0.0;
        for (j = 0; j < REST; j++){
            a = (logTime[j]+logTime[j+1])/2.0;
            da = (logTime[j+1]-logTime[j]);
            //if ((t/a) < 10.0 && (t/a) > 0){
            //    pNew += hitDist[j] * LL_time_spectrum( t/a ) * da/a;
            //}
            if ((t+a) < 10.0 && (t+a) > 0){
                pNew += hitDist[j] * LL_time_spectrum( t + a ) * da;
            }
        }
        convolSpec[i] = pNew;
    }
}

/*calculate the correlation - new spectrum between -3 and 10s*/
/*this is stored in an array so newSpec[0] corresponds to a time of -3s 
and newSpec[1.3*REST-1] to 10s*/
void correlation(double mass2, double dist, double events, double *newSpec, double *logTime, double *logTimeConv){
    double hitDist[REST];
    ProbFirstHitDist(mass2, dist, events, hitDist,logTime);
    /*then the hitDist array will be in log space*/
    convolveHitDistWithLLTimeSpecLog(hitDist, newSpec, logTime, logTimeConv);
    //convolveHitDistWithLLTimeSpec(hitDist, newSpec);
}


void applyEnergyRes(int t, double *distribution, double *energySpectrum){
    /*smear the energy spectrum by convolution with a gaussian*/
    int f, g;
    double pNew;
    for (f=1; f<RESE; f+=1){
        pNew = 0.0;
        for (g=-RESE; g<RESE+1; g+=2){
            if (f-g >= 0 && f-g <= RESE-1){
                pNew += GAUSS(g*STEPE, f*STEPE)*energySpectrum[f-g];
            }
        }
        distribution[t*(RESE-1)+f-1] = pNew*STEPE;
    }
}

int getArrayIndex(double time, double *logTimeConv){
    //determine index of the convoluted time array
    int index = 0;
    if (time > 9.999){
        index = (int) 1.3*REST;
    }
    if (time < -2.9){
        index = 0;
    }
    int i;
    for (i=0;i<1.3*REST;i++){
        if(time < logTimeConv[i]){
            index = i;
            break;
        }
    }
    return index;
}

void getEnergySpec(double mass2, double dist, double *timeArray, double *distribution, double *triggerEffs, bool useEnergyRes, double *logTimeConv){
	double time, pUnsmeared;
	int t, e, arrayIndex;

    for (t=0; t<REST; t++){
        double energySpectrum[RESE];
        energySpectrum[0] = 0.0;
        for (e=1; e<RESE; e++){
            time = getTimeDelay(t*STEPT, e*STEPE, mass2, dist);
            arrayIndex = getArrayIndex(time, logTimeConv);
            //arrayIndex = (int) (time/(STEPT) + 0.3*REST);
            pUnsmeared = LL_energy_spectrum(e*STEPE)*timeArray[arrayIndex]*triggerEffs[(int) (e*STEPE*10)];
            if (!useEnergyRes){
                distribution[t*(RESE-1) +e-1] = pUnsmeared;
            }
            energySpectrum[e] = pUnsmeared;
        }
        if (useEnergyRes){
            applyEnergyRes(t, distribution, energySpectrum);
        }
    }
}

void normalize(double *distribution){
	// normalize the spectrum to 1
	int k;
	double normalize = 0;

	for (k=0; k<(RESE-1)*REST; k++){
		normalize += distribution[k]*STEPT*STEPE;
	}

	for (k=0; k<(RESE-1)*REST; k++){
        	distribution[k] *= 1.0/normalize;
    	}
}

/*generate the proper distribution*/
void generateDist(double mass2, double dist, double events, double *distribution, double *triggerEffs, bool useEnergyRes, double *logTime, double *logTimeConv){
	double timeArray[(int) (1.3*REST)];
	//correlation(mass2, dist, events, timeArray, logTime);
    correlation(mass2, dist, events, timeArray, logTime, logTimeConv);
	getEnergySpec(mass2, dist, timeArray, distribution, triggerEffs, useEnergyRes, logTimeConv);
	normalize(distribution);
}


void fillTriggerEff(double *triggerEffs, bool useTriggerEff){
    /*Read in trigger efficiency.
     Note that the trigger efficiency file for the chosen resolution needs
     to be located in the proper directory.*/
    int i;
    double triggerEns[601];
    if(useTriggerEff){
        FILE *myFile;
        myFile = fopen("trigger_efficiency_100keV_steps.txt", "r");
        for (i = 0; i < 601; i++) {
            fscanf(myFile, "%lf %lf", &triggerEns[i], &triggerEffs[i]);
        }
        fclose(myFile);
    }
    else{
        /* initialize with 1s if trigger efficiency is not considered */
        for(i = 0; i < 601 ; triggerEffs[i++] = 1.0);
    }
}

double noise(double b, double E){
    /*choose noise to be an exponential function at the trigger-threshold (that is at about 5MeV)*/
    double value;
    if (E < 5){
        value = 0.0;
    }
    else{ 
        value = b*exp((5.0-E)*b);
    }
    return value;
}

void addExpNoise(double *spectrum, double b, double events, double noise_events){ 
    /*add noise to the spectrum*/
    int t, e;
    for(t = 0; t < REST; t++){
        for(e = 1; e < RESE; e++){
            spectrum[t*(RESE-1) +e-1] = noise(b, e*STEPE)*noise_events + spectrum[t*(RESE-1) +e-1]*(events-noise_events);
            //spectrum[t*(RESE-1) +e-1] += noise(b, e*STEPE);
        }
    }
    
}

void addNoise(double *spectrum, double noise){
    int i;
    // add constant noise floor to the spectrum
    for (i=0; i<(RESE-1)*REST; i++){
        spectrum[i] += noise;
    }
}

void addNoiseLog(double *spectrum, double noise, double *logTime){
    // add constant noise floor to the spectrum - proper amount depending on size of bin
    int t, e;
    for(t = 0; t < REST; t++){
        for(e = 1; e < RESE; e++){
            spectrum[t*(RESE-1) +e-1] += noise*(logTime[t+1]-logTime[t]);
        }
    }
}

void createSpectrum(double *spectrum, double mass2, double distance, double events, bool useEnergyRes, bool useTriggerEff, double noise, double *logTime, double *logTimeConv){
    /*get trigger efficiencies as function of energy*/
    double triggerEffs[601];
    fillTriggerEff(triggerEffs, useTriggerEff);
    //create a file from the triggerEff for debugging
    /*char filename[sizeof "triggerEff_CUDA.txt"];
    sprintf(filename, "triggerEff_CUDA.txt");
    FILE *f = fopen(filename, "w+");
    int i;
    for(i=0; i< RESE+1 ; i++){
        fprintf(f, "%e\n", triggerEffs[i]);
    }
    fclose(f);*/


    /*create the spectrum from which the random events are drawn*/
    generateDist(mass2, distance, events, spectrum, triggerEffs, useEnergyRes, logTime, logTimeConv);
    /*add noise - constant or expomential*/
    addNoiseLog(spectrum, noise, logTime);
    //addNoise(spectrum, noise);
    //addExpNoise(spectrum, noise, events, noise_events);
    //normalize(spectrum);
}
