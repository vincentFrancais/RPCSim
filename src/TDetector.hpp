/*
 * TDetector.hpp
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

#include <vector>
#include <utility>

#include "TConfig.hpp"

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"
#include "Random.hh"
#include "integration.hpp"

#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_math.h"

#if defined( _PYTHON ) || defined( PYTHON ) || defined (__PYTHON__)
#   ifndef PYTHON
#       define PYTHON
#   endif
#endif

#if defined(PYTHON)
#	include <Python.h>
#endif

using namespace std;
double integrand(double x, void * params);
double Ebar(double x, void * params);

class TDetector{
	
	public:
	TDetector(const TConfig& config);
	
	~TDetector();
	
	void setGasMixture();
	void setGasMixture(Garfield::MediumMagboltz* gas);
	void makeGasTable();
	void setElectricField(const double& Ex ,const double& Ey, const double& Ez);
	Garfield::MediumMagboltz* getGas();
	void initialiseDetector();
	void writeGasTransportParameters();
	
	static void printPACSData(Garfield::MediumMagboltz* gas);
	void printPACSData();
	
	void setGarfieldSeed( const int& s );
	
	double R(const double& k, const double& z, const double& zp);
	double D(const double& k);
	double SCPotential(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCFieldSimplified(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCFieldSimplified2(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	/* SCFieldSimplfied3,4,5 and 6 are used to output the terms of the equation in an easy-plotting structure
	 * I guess they need better names ... But I'm lazy right now ...... */
	double SCFieldSimplified3(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCFieldSimplified4(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCFieldSimplified5(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCFieldSimplified6(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCField(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	
	inline double RadialChargeDistribution(const double& r, const double& l);
	double computeEbar(const double& z, const double& l, const double& zp);
	#if defined(PYTHON)
		double computeEbar_Python(const double& z, const double& l, const double& zp);
	#endif
	void makeEbarTable( bool const& binary=true );
	void printEbarTable();
	
	//TConfig getConfig(void) const {return fConfig;}
	//double getGapWidth(void) const	{return fGeometry.gapWidth;}
	//const double* getResistiveLayersWidth(void) const {return fGeometry.resistiveLayersWidth;}
	//DetectorGeometry getGeometry(void) const {return fGeometry;}
	string getGasName() const;
	double getTimeStep(void) const	{return fDt;}
	double getSpaceStep(void) const	{return fDx;}
	double getDiffL(void) const	{return fDiffL;}
	double getDiffT(void) const	{return fDiffT;}
	double getVx(void) const	{return fVx;}
	double getVy(void) const	{return fVy;}
	double getVz(void) const	{return fVz;}
	double getAlpha(void) const {return fAlpha;}
	double getEta(void) const	{return fEta;}
	double getK(void) const	{return k;}
	int getEbarTableSize(void) const {return iEbarTableSize;}
	string getEbarTableHexName(void) const {return fEbarTableHexName;}
	
	bool hasEBarTable(void) const {return bHasEbarTable;}

	vector<double> getEbarVecTable(void) const {return fEbarVecTable;}
	int getNstep(void) const	{return iNstep;}
	double* getElectricField(void) {return fElectricField;}
	vector<double> getTransportParameters(double Ex, double Ey, double Ez);
	vector<double> getDiffusionCoefficients( double const& Ex, double const& Ey, double const& Ez );
	Garfield::Sensor* getSensor(void) const	{return mSensor;}

	vector<double> getEbarZarray(void) {return fEbarZarray;}
	vector<double> getEbarZparray(void) {return fEbarZparray;}
	vector<double> getEbarLarray(void) {return fEbarLarray;}
	
	void plotSC();
	
	string getUniqueTableName();
	
	private:
	int iNstep;
	TConfig fConfig;
	
	double fDt;
	double fDx;
	
	int iEbarTableSize;
	vector<double> fEbarVecTable;
	vector<double> fEbarZarray, fEbarZparray, fEbarLarray;
	
	bool bHasEbarTable;
	bool bGasLoaded;
	bool bDetectorInitialised;
	string fEbarTableHexName;
	
	Garfield::MediumMagboltz* mGas;
	string mGasTableName;
	Garfield::Sensor* mSensor;
	
	double fTemperature;
	double fPressure;
	double fElectricField[3];
	
	double fDiffL;
	double fDiffT;
	double fVx;
	double fVy;
	double fVz;
	double fiVx, fiVy, fiVz;
	double fAlpha;
	double fEta;
	double k;
};
