#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

#include "TDetector.hpp"
#include "TResult.hpp"
#include "helper_functions.hpp"

using namespace std;

class TAvalanche {
	public:
	TAvalanche(){};
	TAvalanche(TDetector* det, bool const& randomSeed=false);
	
	~TAvalanche();
	
	static void computeElectronsProduction(const TDetector* det, const string& particleName, const double& P, const int& nTracks);
	static void computeClusterDensity(const TDetector* det, const string& particleName, const double& Pmin, const double& Pmax, const int& steps);
	
	protected:
	TDetector* fDet;
	DetectorGeometry fGeometry;
	
	RngStream* fRandRng;
	RngStream* fRandRngCLT;
	RngStream* fRandRngLongiDiff;
	
	TResult fResult;
};
