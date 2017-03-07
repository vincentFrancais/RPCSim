
#pragma once

#include <iostream>
#include <vector>

#include "tinyxml2.h"


using namespace std;

struct TConfig {
	
	TConfig (string file);
	TConfig () {};
	~TConfig () {};
	
	void print();
	
	//private:
	double gapWidth;
	int nSteps;
	
	double anodeWidth, cathodeWidth;
	double anodePermittivity, cathodePermittivity;
	
	int nResisLayers;
	vector<double> resisLayersWidth;
	vector<double> resisLayersEpsilon;
	
	double ElectricField;
	
	int EbarTableCalculationSteps;
	
	int nGases;
	vector<string> gasNames;
	vector<double> gasPercentage;
	double gasTemperature;
	double gasPressure;
	
	double threshold;
	
	string particleName;
	double particleMomentum;
	double x0;
	double theta;
	
	int nThreads;
	int nEvents;
	string outFile;
	
	string generator;
	vector<int> seeds;
	int garfieldSeed;
	int verbosityLevel;
	bool verbose;
	bool snaps;
	
	bool singleCluster;
	int n0;
	//bool randomSingleCluster;
	
	bool clusterDensity;
	bool electronProduction;
	bool noAvalanche;
	bool computeEfficiency;
	bool useUUIDSeed;
	bool onlyMult;
	bool debugOutput;

	//TConfig() : nResisLayers(0), nGases(0), outFile("out/"), generator("SFMT"), garfieldSeed(-1) { }
};

//TConfig readConfigFile(string fileName);
//void printConfig(TConfig config);
