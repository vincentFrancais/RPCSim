
#pragma once

#include <iostream>
#include <vector>

#include "tinyxml2.h"


using namespace std;

struct Config {
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
	double particleEnergy;
	double x0;
	double theta;
	
	int nThreads;
	int nEvents;
	Config() : nResisLayers(0), nGases(0) { }
};

Config readConfigFile(string fileName);
void printConfig(Config config);
