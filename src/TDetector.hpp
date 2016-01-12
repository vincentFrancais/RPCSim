#ifndef DEF_DETECTOR
#define DEF_DETECTOR

#include <vector>
#include <utility>

//#include <TRandom3.h>
#include "RngStream.h"

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"

using namespace std;
double integrand(double x, void * params);
double Ebar(double x, void * params);
double test_int(double x, void * params);

struct DetectorGeometry{
	//DetectorGeometry(double gap=0.2, double resLayer1Width=0.2, double resLayer2Width=0.2) {
		//gapWidth = gap;
		//resistiveLayersWidth[0] = resLayer1Width;
		//resistiveLayersWidth[1] = resLayer2Width;
		//}
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
	
	double R(const double& k, const double& z, const double& zp);
	double D(const double& k);
	double SCPotential(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCFieldSimplified(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCField(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double SCField2(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp);
	double test();
	
	double RadialChargeDistribution(const double& r, const double& l);
	double computeEbar(const double& z, const double& l, const double& zp);
	double computeEbar_Python(const double& z, const double& l, const double& zp);
	void makeEbarTable();
	
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

	vector<double> getEbarVecTable(void) const {return fEbarVecTable;}
	int getNstep(void) const	{return iNstep;}
	double* getElectricField(void) {return fElectricField;}
	double* getTransportParameters(double Ex, double Ey, double Ez);
	Garfield::Sensor* getSensor(void) const	{return mSensor;}

	vector<double> getEbarZarray(void) {return fEbarZarray;}
	vector<double> getEbarZparray(void) {return fEbarZparray;}
	vector<double> getEbarLarray(void) {return fEbarLarray;}
	
	void plotSC();
	
	private:
	int iNstep;
	//double fGapWidth;
	//double[2] fResistiveLayersWidth;
	DetectorGeometry fGeometry;
	double fDt;
	double fDx;
	
	int iEbarTableSize;
	vector<double> fEbarVecTable;
	vector<double> fEbarZarray, fEbarZparray, fEbarLarray;
	
	bool bEbarComputed;
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

