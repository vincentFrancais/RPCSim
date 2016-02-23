#ifndef DEF_AVALANCHE
#define DEF_AVALANCHE

#include <vector>
#include <utility>
//#include <cstdint>

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
#include "TRNQueue.hpp"

#if defined( _DEBUG ) || defined( DEBUG ) || defined (__DEBUG__)
#   ifndef DEBUG
#       define DEBUG
#   endif
#endif

#if defined(DEBUG)
#	include <assert.h>
#	define DEBUGASSERT                assert
#else
#   define DEBUGASSERT( x )               {}
#endif

using namespace std;

enum EAvalancheStatus{
	AVAL_SUCCESS,
	AVAL_NO_ERROR,
	
	AVAL_CLT_FAIL,
	AVAL_MULT_FAIL,
	AVAL_ERROR_GRID_NOT_EMPTY,
	AVAL_ERROR_TIMESTEP_EXCEEDING_LIMIT,
	AVAL_EXPLOSIVE_BEHAVIOR
};

class TAvalanche{
	
	public:
	TAvalanche(){};
	TAvalanche(TDetector* det, bool const& randomSeed=false);
	
	~TAvalanche();
	
	void initialiseTrackHeed(const string& particleName, const double& momentum, const double& x0, const double& theta);
	
	void simulateEvent();
	void initialiseSingleCluster(const double& x0, const double& n0=1);
	void computeInducedCharges();
	void computeInducedSignal();
	void computeInducedSignal2();
	void computeChargeSpectra(const double& x0, const double& angle, const int& nEvents);
	void computeClusterDensity(const string& particleName, const double& Pmin, const double& Pmax, const int& steps);
	void computeElectronsProduction(const string& particleName, const double& P, const int& nTracks);
	void checkDetectorGrid();
	
	vector<double> getInducedSignal() const {return fSignal;}
	vector<double> getInducedCharges() const {return fCharges;}
	TResult getResultFile() const {return fResult;}
	
	void printDetectorGrid() const;
	
	void makeResultFile();
	
	double n_moy(const double& x);
	double electron_multiplication(const double& x, const double& s);
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
	uint tId;
	
	DetectorGeometry fGeometry;
	TDetector* fDet;
	
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
	
	//RNQueue fRNQueue;
	
	bool bVerbose;
	
	bool bFullLongiDiff;
	
	bool bPrintDetectorGrid;
	bool bTestPot;
	bool bEbarComputed;
	bool bSnapshots;
	
	bool bThrCrossTime;
	
	EAvalancheStatus eAvalStatus;
};

#endif
