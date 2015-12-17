#ifndef DEF_AVALANCHE
#define DEF_AVALANCHE

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

#include "TDetector.hpp"
#include "TResult.hpp"

using namespace std;

class TAvalanche{
	
	public:
	TAvalanche(){};
	TAvalanche(TDetector* det);
	
	~TAvalanche();
	
	void initialiseTrackHeed(TDetector* det, const string& particleName, const double& momentum, const double& x0, const double& theta);
	
	void simulateEvent();
	void initialiseSingleCluster(const double& x0, const double& n0=1);
	void computeInducedCharges();
	void computeInducedSignal();
	void computeInducedSignal2();
	void computeChargeSpectra(const double& x0, const double& angle, const int& nEvents);
	void computeClusterDensity(TDetector* det, const string& particleName, const double& Pmin, const double& Pmax, const int& steps);
	void checkDetectorGrid();
	
	vector<double> getInducedSignal() const {return fSignal;}
	vector<double> getInducedCharges() const {return fCharges;}
	TResult getResultFile() const {return fResult;}
	
	void printDetectorGrid() const;
	
	void makeResultFile();
	
	double n_moy(const double& x);
	double electron_multiplication2(const double& x, const double& s);
	double multiplication(const double& n);
	double CLT(const double& x, const double& n);
	void computeLongitudinalDiffusion();
	bool avalanche();
	
	void computeSCEffect();
	double interpolateEbar(const double& z, const double& zp, const double& l);
	
	void makeSnapshot();
	
	void testInterpolation();
	size_t getIndex3D(const size_t& i, const size_t& j, const size_t& k);
	
	private:
	int iNstep;
	double fGapWidth;
	const double* fResistiveLayersWidth;
	double fDt;
	double fDx;
	
	int iEbarTableSize;
	vector<double> fEbarTable;
	string fEbarTableHexName;
	vector<double> fEbarZarray;
	vector<double> fEbarZparray;
	vector<double> fEbarLarray;
	
	int iCurrentDetectorStep;
	int iTimeStep;
	int iThrCrossTimeStep;
	
	double fChargeThres;
	
	DetectorGeometry fGeometry;
	TDetector* fDet;
	
	int iNElectronsSize;
	vector<double> fElecDetectorGrid;
	vector<double> fPosIonDetectorGrid;
	vector<double> fNegIonDetectorGrid;
	
	double fNelecAnode;
	vector< pair<double,double> > fAnodeLValues;
	
	double* fLongiDiffFrac;
	//vector<double> fNElectrons;
	double * fNElectrons;
	vector<double> fSignal;
	vector<double> fCharges;
	vector<double> fTotalCharges;
	
	double fDiffL;
	double fDiffT;
	double* fVx;
	double fVy;
	double fVz;
	double fiVx, fiVy, fiVz;
	double* fAlpha;
	double* fEta;
	double* fE;
	double fEini;
	//double k;
	
	double fClusterDensity;
	
	vector<double> fClPosX, fClPosY, fClPosZ;
	
	double fLongiDiffSigma;
	
	double fThrCLT;
	double fSpaceChargeLimit;
	double fLongiDiffComputeLimit;
	
	TResult fResult;
	
	RngStream* fRandRng;
	RngStream* fRandRngCLT;
	RngStream* fRandRngLongiDiff;
	RngStream* fRandRngSCE;
	
	int iDebug;
	
	bool bFullLongiDiff;
	
	bool bPrintDetectorGrid;
	bool bTestPot;
	bool bEbarComputed;
	
	bool bVerbose;
	bool bThrCrossTime;
};

#endif
