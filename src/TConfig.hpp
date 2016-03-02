
#pragma once

#include <iostream>
#include <vector>

#include "tinyxml2.h"


using namespace std;

struct TConfig {
	double gapWidth;
	int nSteps;
	
	int nResisLayers;
	vector<double> resisLayersWidth;
	vector<double> resisLayersEpsilon;
	double ElectricField;
	
	int nGases;
	vector<string> gasNames;
	vector<double> gasPercentage;
	double gasTemperature;
	double gasPressure;
	
	string particleName;
	double particleMomentum;
	double x0;
	double theta;
	
	int nThreads;
	int nEvents;
	string outFile;
	TConfig() : nResisLayers(0), nGases(0), outFile("out/") { }
};

TConfig readConfigFile(string fileName);
void printConfig(TConfig config);
