#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <utility>
#include <sys/stat.h>
#include <typeinfo>
#include <fcntl.h>
#include <cassert>

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"
#include "Plotting.hh"

#include "TDetector.hpp"
#include "TAvalanche1D.hpp"
#include "TAvalanche1D.hpp"
#include "helper_functions.hpp"
#include "TThreadsFactory.h"
#include "TResult.hpp"
#include "TConfig.hpp"

//#define PATH_MAX 0x1000
#define PI 3.14159265358979323846

using namespace std;
using namespace Garfield;

typedef unsigned int uint;
typedef unsigned long ulong;

double cm = 0.01;

pthread_mutex_t gPipeLock;
pthread_mutex_t gTrackLock;
int gPipe[2];
TResult gNullResult;

struct ThreadData{
	TDetector * detector;
	TConfig config;
	
	ThreadData (TDetector * det, TConfig conf) : detector(det), config(conf){ };
};

void * wrapperFunction(void * Arg){
	
	assert(Arg != NULL);

	//TDetector * detector = reinterpret_cast<TDetector *>(Arg);
	ThreadData* data = reinterpret_cast< ThreadData* > (Arg);
	
	TResult result;
	TAvalanche1D avalanche(data->detector, false);
	
	
	sem_post(TThreadsFactory::GetInstance()->GetInitLock());
	
	pthread_mutex_lock(&gTrackLock);
	avalanche.initialiseTrackHeed(data->config.particleName, data->config.particleMomentum, data->config.x0, data->config.theta);
	//avalanche.initialiseSingleCluster(0);
	pthread_mutex_unlock(&gTrackLock);
	
	//avalanche.disableSpaceChargeEffect();
	avalanche.simulateEvent();
	
	result = avalanche.getResultFile();
	
	pthread_mutex_lock(&gPipeLock);
    write(gPipe[1], &result, sizeof(result));
    pthread_mutex_unlock(&gPipeLock);
	
	return 0;
}

void * WriteResults(void * Arg)
{
    TResult result;
    char * outputFile = reinterpret_cast<char *>(Arg);

    int outFD;
	
    /* Open the output file */
    string outFileBinary (outputFile);
    size_t found = outFileBinary.find(".dat");
	if (found!=string::npos)
		outFileBinary.erase(found,4);
    outFileBinary += ".bin";

    outFD = open(outFileBinary.c_str(), O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    assert(outFD != -1);

    /* Read the incoming event */
    while (read(gPipe[0], &result, sizeof(TResult)) > 0)
    {
        /* End signal => stop dealing with data */
        if (memcmp(&result, &gNullResult, sizeof(TResult)) == 0)
            break;

        /* Write the event to the output file */
        write(outFD, &result, sizeof(TResult));
    }

    close(outFD);

    return 0;
}

int main(int argc, char** argv) {
	/* Read config file */
	TConfig config = readConfigFile("config/calice.xml");
    printConfig(config);
     
    char outputFile[PATH_MAX];
    if(argc > 2)
		strncpy(outputFile, argv[2], PATH_MAX - 1);
	else
		strncpy(outputFile, config.outFile.c_str(), PATH_MAX - 1);
    outputFile[PATH_MAX - 1] = '\0';
    
    pthread_t writingThread;
    void * ret;
    
    unsigned int nThreads = config.nThreads;
    unsigned long nEvents = config.nEvents;
		
    
    /* Initialize our pipe lock */
    pthread_mutex_init(&gPipeLock, 0);
    pthread_mutex_init(&gTrackLock, 0);

    /* Start our threads factory */
    TThreadsFactory::GetInstance()->SetMaxThreads(nThreads);
    
    /* Init our null event */
    memset(&gNullResult, 0, sizeof(TResult));
    
    /* Init our detector */
	MediumMagboltz* gas = new MediumMagboltz();
	switch (config.nGases){
		case (1): 
			gas->SetComposition(config.gasNames[0], config.gasPercentage[0]);
			break;
		case (2):
			gas->SetComposition(config.gasNames[0], config.gasPercentage[0], config.gasNames[1], config.gasPercentage[1]);
			break;
		case (3):
			gas->SetComposition(config.gasNames[0], config.gasPercentage[0], config.gasNames[1], config.gasPercentage[1], config.gasNames[2], config.gasPercentage[2]);
			break;
		case (4):
			gas->SetComposition(config.gasNames[0], config.gasPercentage[0], config.gasNames[1], config.gasPercentage[1], config.gasNames[2], config.gasPercentage[2], config.gasNames[3], config.gasPercentage[3]);
			break;
		case (5):
			gas->SetComposition(config.gasNames[0], config.gasPercentage[0], config.gasNames[1], config.gasPercentage[1], config.gasNames[2], config.gasPercentage[2], config.gasNames[3], config.gasPercentage[3], config.gasNames[4], config.gasPercentage[4]);
			break;
		case (6):
			gas->SetComposition(config.gasNames[0], config.gasPercentage[0], config.gasNames[1], config.gasPercentage[1], config.gasNames[2], config.gasPercentage[2], config.gasNames[3], config.gasPercentage[3], config.gasNames[4], config.gasPercentage[4], config.gasNames[5], config.gasPercentage[5]);
			break;
	}

	gas->SetTemperature(config.gasTemperature);
	gas->SetPressure(config.gasPressure);
	
	double HV = config.ElectricField;
	if(argc > 1)	HV = atof(argv[1])*1000.;
	
	if ( argc > 1 )
		cout << "Efficiency computation runs. HV=" << HV << " OutFile=" << outputFile << endl;
	
	DetectorGeometry geom;
	geom.gapWidth = 0.12;	//0.2; cm	//CALICE 0.12
	geom.resistiveLayersWidth[0] = 0.11;	//0.2;	//CALICE 0.11
	geom.resistiveLayersWidth[1] = 0.07;	//0.2;	//CALICE 0.07
	geom.relativePermittivity[0] = 7.;	//10.;	//CALICE 7
	geom.relativePermittivity[1] = 7.;	//10.;	//CALICE 7
	TDetector* detector = new TDetector(geom,500);
	detector->setGasMixture(gas);
	detector->setElectricField(HV,0.,0.);
	detector->initialiseDetector();
	detector->makeEbarTable();
	//detector->setGarfieldSeed( 2844356529 );
	
	ThreadData* data = new ThreadData(detector, config);
	
	//TAvalanche* avalanche = new TAvalanche(detector);
	//avalanche->computeClusterDensity(detector, "muon", 6.e7, 15.e9, 600);
	//avalanche->computeElectronsProduction(detector, "muon", 5.e9, 100000);
	//delete avalanche;
	//exit(0);
	
    /* Open the communication pipe */
    if (pipe(gPipe) == -1){
		goto end;
    }

    /* Start our background writing thread */
    if (pthread_create(&writingThread, 0, WriteResults, outputFile) != 0){
		goto end2;
    }

    /* Hot loop, the simulation happens here */
    for (unsigned long i = 0; i < nEvents; ++i){
		TThreadsFactory::GetInstance()->CreateThread(wrapperFunction, data);
    }

    /* Wait for all the propagations to finish */
    TThreadsFactory::GetInstance()->WaitForAllThreads();

    /* Send the end signal to writer */
    write(gPipe[1], &gNullResult, sizeof(TResult));

    /* Wait for the end of the writer */
    pthread_join(writingThread, &ret);

end2:
    close(gPipe[0]);
    close(gPipe[1]);
end:
    pthread_mutex_destroy(&gPipeLock);
    pthread_mutex_destroy(&gTrackLock);
    delete detector;

    return 0;
}
