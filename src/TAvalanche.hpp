/*
 * TAvalanche.hpp
 * 
 * Copyright 2016 Vincent Français <francais@clermont.in2p3.fr>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

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
	static int countSim;
	//TDetector* fDet;
	//DetectorGeometry fGeometry;
	
	//RngStream* fRandRng;
	//RngStream* fRandRngCLT;
	//RngStream* fRandRngLongiDiff;
	
	//TResult fResult;
	
	TTimer fTimer;
};
