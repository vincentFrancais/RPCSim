// TAvalanche2D.hpp
// 
// Copyright 2017 Vincent Français <francais@clermont.in2p3.fr>
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301, USA.
// 
// 


#pragma once

#include <vector>
#include <utility>
#include <map>
#include <chrono>

#include "TRandomEngineSFMT.hpp"
#include "TRandomEngineMRG.hpp"
#include "TRandomEngineMT.hpp"
#include "TRandomEngineMTDC.hpp"
#include "TRandomEngine.hpp"

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"

#include "TAvalanche.hpp"
#include "TConfig.hpp"
#include "TDetector.hpp"
#include "TResult.hpp"
#include "TAvalError.hpp"

using namespace std;

class TAvalanche2D : public TAvalanche {
	public:
	TAvalanche2D(TDetector* det, TConfig& config, const sfmt_t sfmt, const int& id);
	
	inline double indexDetectorGrid(const int& r, const int& z);
	double electricFieldStrength (const int& r, const int&z);
	void init();
	
	inline double n_moy(const double& x);
	double multiplicationRiegler(const double& x, const double& s);
	double electronMultiplication(const double& n);
	double multiplicationCLT(const double& x, const double& n);
	
	void updateParameters();
	bool propagate();
	bool avalanche();
	
	private:
	TDetector* fDetector;
	TConfig fConfig;
	//DetectorGeometry fGeometry;
	
	int fNstepz, fNstepr, fGridSize;
	double fDz, fDr;
	double fCurrentStepR,fCurrentStepZ,fCurrentStep;
	double fCurrentEta, fCurrentAlpha, fCurrentVel;
	int fTimeStep;
	
	vector<double> fElecDetectorGrid, fPosIonDetectorGrid, fNegIonDetectorGrid;
	vector<double> fVelocity, fAlpha, fEta;
	vector<double> fEz, fEr;
	
	double fThrCLT;
	
	EAvalancheStatus fAvalStatus;
	
	TRandomEngine* fRngMult;
	TRandomEngine* fRngCLT;
	TRandomEngine* fRngLongiDiff;
	
	bool fVerbose;
};
