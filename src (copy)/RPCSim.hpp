#ifndef DEF_RPCSIM
#define DEF_RPCSIM

#include <vector>
#include <utility>

#include <TRandom3.h>
#include "RngStream.h"

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"

using namespace std;

class RPCSim{
	
	public:
	RPCSim(){};
	RPCSim(const int& nStep,const double& gasGap);
	void setGasMixture(Garfield::MediumMagboltz* gas);
	void makeGasTable();
	void setElectricField(const double& Ex ,const double& Ey, const double& Ez);
	Garfield::MediumMagboltz* getGas();
	void initialiseDetector();
	void setTrackHeed(const string& particleName, const double& momentum, const double& x0, const double& theta);
	
	void simulateEvent();
	void singleAvalanche(const double& x0);
	void computeInducedCharges();
	void computeInducedSignal();
	void computeChargeSpectra(const double& x0, const double& angle, const int& nEvents);
	void checkDetectorGrid();
	
	vector<double> getInducedSignal(){return fSignal;}
	vector<double> getInducedCharges(){return fCharges;}
	Garfield::TrackHeed* getTrackHeed(){return mTrack;}
	
	double n_moy(const double& x);
	double electron_multiplication2(const double& x, const double& s);
	double multiplication(const double& n);
	double CLT(const double& x, const double& n);
	void computeLongitudinalDiffusion();
	bool avalanche();
	
	
	
	private:
	int iNstep;
	double fGapWidth;
	double fDt;
	double fDx;
	
	Garfield::MediumMagboltz* mGas;
	string mGasTableName;
	Garfield::Sensor* mSensor;
	Garfield::TrackHeed* mTrack;
	
	double fTemperature;
	double fPressure;
	double fElectricField[3];
	
	double* fDetectorGrid;
	vector<double> fNElectrons;
	vector<double> fSignal;
	vector<double> fCharges;
	
	double fDiffL;
	double fDiffT;
	double fVx;
	double fVy;
	double fVz;
	double fAlpha;
	double fEta;
	double k;
	
	double fThrCLT;
	double fSpaceChargeLimit;
	double fLongiDiffComputeLimit;
	
	//TRandom3 fRand;
	//TRandom3 fRandCLT;
	//TRandom3 fRandLongiDiff;
	
	RngStream fRandRng;
	RngStream fRandRngCLT;
	RngStream fRandRngLongiDiff;
	
	int iDebug;
};

#endif
