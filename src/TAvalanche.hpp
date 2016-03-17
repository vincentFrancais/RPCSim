#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

#include "TDetector.hpp"
#include "helper_functions.hpp"
#include "TTimer.hpp"

using namespace std;

class TAvalanche {
	public:
	//TAvalanche(){};
	TAvalanche();
	
	~TAvalanche();
	
	static void computeElectronsProduction(const TDetector* det, const string& particleName, const double& P, const int& nTracks);
	static void computeClusterDensity(const TDetector* det, const string& particleName, const double& Pmin, const double& Pmax, const int& steps);
	static int getCount() {return count;}
	
	protected:
	static int count;
	//TDetector* fDet;
	//DetectorGeometry fGeometry;
	
	//RngStream* fRandRng;
	//RngStream* fRandRngCLT;
	//RngStream* fRandRngLongiDiff;
	
	//TResult fResult;
	
	TTimer fTimer;
};
