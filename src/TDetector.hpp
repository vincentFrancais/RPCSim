#ifndef DEF_DETECTOR
#define DEF_DETECTOR

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

struct DetectorGeometry{
	/* TODO: change array to anode and cathode variables -- more readable */
	double gapWidth;
	double resistiveLayersWidth[2];
	double relativePermittivity[2];
};

class TDetector{
	
	public:
	TDetector(const DetectorGeometry& geometry, const int& nStep = 500);
	
	~TDetector();
	
	void setGasMixture(Garfield::MediumMagboltz* gas);
	void makeGasTable();
	void setElectricField(const double& Ex ,const double& Ey, const double& Ez);
	Garfield::MediumMagboltz* getGas();
	void initialiseDetector();
	void writeGasTransportParameters();
	
	void setGarfieldSeed( const int& s );
	
	double R(const double& k, const double& z, const double& zp);
	double D(const double& k);
	double SCPotential(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCFieldSimplified(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCField(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	
	double RadialChargeDistribution(const double& r, const double& l);
	double computeEbar(const double& z, const double& l, const double& zp);
	#if defined(PYTHON)
		double computeEbar_Python(const double& z, const double& l, const double& zp);
	#endif
	void makeEbarTable( bool const& binary=true );
	
	double getGapWidth(void) const	{return fGeometry.gapWidth;}
	const double* getResistiveLayersWidth(void) const {return fGeometry.resistiveLayersWidth;}
	DetectorGeometry getGeometry(void) const {return fGeometry;}
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
	double* getTransportParameters(double Ex, double Ey, double Ez);
	double* getDiffusionCoefficients( double const& Ex, double const& Ey, double const& Ez );
	Garfield::Sensor* getSensor(void) const	{return mSensor;}

	vector<double> getEbarZarray(void) {return fEbarZarray;}
	vector<double> getEbarZparray(void) {return fEbarZparray;}
	vector<double> getEbarLarray(void) {return fEbarLarray;}
	
	void plotSC();
	
	string getUniqueTableName(int const& n); //TO BE REMOVED!
	string getUniqueTableName();
	
	private:
	int iNstep;

	DetectorGeometry fGeometry;
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

#endif

