
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <assert.h> 
#include <utility>
#include <ncurses.h>
#include <limits>
#include <iterator>

#include "TDetector.hpp"
#include "TAvalanche.hpp"
#include "helper_functions.hpp"
#include "gsl/gsl_const_mksa.h"
#include "gsl/gsl_sf_bessel.h"
#include <gsl/gsl_integration.h>
#include "gsl/gsl_math.h"
#include <gsl/gsl_randist.h>

#define PI 3.14159265358979323846

using namespace std;
//static  TAvalanche * tgsl = 0;

TAvalanche::TAvalanche(TDetector* det){
	iNstep = det->getNstep();
	
	//const int iEbarTableSize = 50;
	//fEbarTable = new double[iEbarTableSize][iEbarTableSize][iEbarTableSize];
	
	fRandRng = new RngStream();
	fRandRngCLT = new RngStream();
	fRandRngLongiDiff = new RngStream();
	//fRandRngTransDiff = new RngStream();
	
	fGSLRng = gsl_rng_alloc (gsl_rng_taus);
	
	fGapWidth = det->getGapWidth();
	fResistiveLayersWidth = det->getResistiveLayersWidth();
	fDt = det->getTimeStep();
	fDx = det->getSpaceStep();
	
	fGeometry = det->getGeometry();
	
	fDiffL = det->getDiffL();
	fDiffT = det->getDiffT();
	//fVx = det->getVx();
	fVx = new double[iNstep];
	fVy = det->getVy();
	fVz = det->getVz();
	//fAlpha = det->getAlpha();
	fAlpha = new double[iNstep];
	fEta = new double[iNstep];
	fE = new double[iNstep];
	for(int i=0; i<iNstep; i++){
		fVx[i] = det->getVx();
		fAlpha[i] = det->getAlpha();
		fEta[i] = det->getEta();
		fE[i] = det->getElectricField()[0];
	}
	fEini = fE[0];
	//fEta = det->getEta();
	//k = det->getK();
	
	fThrCLT = 1.e4;
	fSpaceChargeLimit = 1.e19;//5.e7;
	fLongiDiffComputeLimit = 100;
	
	iNElectronsSize = 5*iNstep;
	fElecDetectorGrid = new double[iNstep];
	fPosIonDetectorGrid = new double[iNstep];
	fNegIonDetectorGrid = new double[iNstep];
	
	fLongiDiffSigma = fDiffL*sqrt(fDx);
	fLongiDiffFrac = new double[7];
	fLongiDiffFrac[0] = gaussIntegral(2000, -0.5*fDx, 0.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[1] = gaussIntegral(2000, 0.5*fDx, 1.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[2] = gaussIntegral(2000, -1.5*fDx, -0.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[3] = gaussIntegral(2000, 1.5*fDx, 2.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[4] = gaussIntegral(2000, -2.5*fDx, -1.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[5] = gaussIntegral(2000, 2.5*fDx, 3.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[6] = gaussIntegral(2000, -3.5*fDx, -2.5*fDx, 0, fLongiDiffSigma);
	
	bPrintDetectorGrid = false;
	bTestPot = false;
	bFullLongiDiff = true;
	//bEbarComputed = false;
	
	fDet = det;
	
	iEbarTableSize = det->getEbarTableSize();
	int n = iEbarTableSize;
	fEbarTable = new double(n*n*n);
	//fEbarTable[0]
	fEbarTable = det->getEbarTable();
	cout << fEbarTable[0] << endl;
	
	fEbarZarray = new double[iEbarTableSize];
	fEbarZparray = new double[iEbarTableSize];
	fEbarLarray = new double[iEbarTableSize];
	
	double cm=0.01;
	double zStep = fGapWidth * cm / (iEbarTableSize), zpStep = fGapWidth * cm / (iEbarTableSize), lStep = iNstep*fDx / iEbarTableSize;
	
	for(int i=0; i<iEbarTableSize; i++){
		fEbarZarray[i] = (i+1)*zStep;
		fEbarZparray[i] = (i+1)*zpStep;
		fEbarLarray[i] = (i+1)*lStep;
		
	}
	//testInterpolation();
	//exit(0);
}

TAvalanche::~TAvalanche(){
	delete fNElectrons;
	delete fLongiDiffFrac;
}

void TAvalanche::computeClusterDensity(TDetector* det, const string& particleName, const double& Pmin, const double& Pmax, const int& steps){
	double P = Pmin;
	
	string outFileName = "out/cluster_density_"+particleName+".dat";
	ofstream data(outFileName.c_str(), ios::out | ios::trunc);
	
	Garfield::TrackHeed* track = new Garfield::TrackHeed();
	track->SetSensor(det->getSensor());
	track->SetParticle(particleName);
	
	while(P <= Pmax){
		track->SetMomentum(P);
		track->NewTrack(0, 0, 0, 0, 0, 0, 0);
		data << P << "\t" << track->GetClusterDensity() << endl;
		P += (Pmax-Pmin)/steps;
	}
	data.close();
	delete track;
}

void TAvalanche::initialiseTrackHeed(TDetector* det, const string& particleName, const double& momentum, const double& x0, const double& theta){
	Garfield::TrackHeed* track = new Garfield::TrackHeed();
	track->SetSensor(det->getSensor());
	track->SetParticle(particleName);
	track->SetMomentum(momentum);
	
	double t0 = 0.;
	double y0 = 0, z0 = 0; 
	double dx0 = cos(theta * PI / 180.0), dy0 = sin(theta * PI / 180.0), dz0 = 0.;
	
	track->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
	
	double xc = 0., yc = 0., zc = 0., tc = 0.; 	// Cluster coordinates
	int nc = 0; 								// Number of electrons produced in a collision
	double ec = 0.; 							// Energy loss in a collision
	double dummy = 0.; 							// Dummy variable (not used at present)
	fNElectrons = new double[iNElectronsSize];
	for(int i=0; i<iNElectronsSize; i++)	fNElectrons[i] = 0;
	for(int i=0; i<iNstep; i++){
		fElecDetectorGrid[i] = 0;
		fPosIonDetectorGrid[i] = 0;
		fNegIonDetectorGrid[i] = 0;
	}
	
	while (track->GetCluster(xc, yc, zc, tc, nc, ec, dummy)){
		fClPosX.push_back(xc);
		fClPosY.push_back(yc);
		fClPosZ.push_back(zc);
		fNElectrons[0] += nc;
		fElecDetectorGrid[int(trunc(xc/fDx))] += nc;
		fPosIonDetectorGrid[int(trunc(xc/fDx))] += nc;
	}
	
	fClusterDensity = track->GetClusterDensity();
	delete track;
}

void TAvalanche::initialiseSingleCluster(const double& x0, const double& n0){
	
	//Clusters clusters = Clusters();
	fNElectrons = new double[iNElectronsSize];
	for(int i=0; i<iNElectronsSize; i++)	fNElectrons[i] = 0;
	for(int i=0; i<iNstep; i++)	fElecDetectorGrid[i] = 0;
	
	fNElectrons[0] = n0;
	fElecDetectorGrid[int(trunc(x0/fDx))] = n0;
	
}

void TAvalanche::makeResultFile(){
	fResult.Dt = fDt;
	fResult.Dx = fDx;
	//fResult.alpha = fAlpha;
	//fResult.eta = fEta;
	fResult.fDiffCoeff[0] = fDiffL;
	fResult.fDiffCoeff[1] = fDiffT;
	fResult.fGapWidth = fGapWidth;
	fResult.fInducedCharge = fCharges.back();
	fResult.fInducedSignal = vecToArray(fResult.fInducedSignal, fSignal);
	fResult.fNElectrons = new double[iNElectronsSize];
	memcpy(fResult.fNElectrons, fNElectrons, sizeof(double)*iNElectronsSize);
	//fResult.fClPosX = fClPosX;
	//fResult.fClPosY = fClPosY;
	//fResult.fClPosZ = fClPosZ;
	fResult.fSpaceChargeLimitThr = fSpaceChargeLimit;
	fResult.fThrCLT = fThrCLT;
	//fResult.fVx = fVx;
	fResult.iNstep = iNstep;
	//fResult.sMixture = 
}

void TAvalanche::simulateEvent(){
	
	//makeEbarTable();
	//PlotField();
	//testInterpolation();
	//exit(0);
	
	if( avalanche() ){
		checkDetectorGrid();
		computeInducedSignal();
		computeInducedCharges();
		ofstream data("out/electrons.dat", ios::out | ios::trunc);
		for(int i=0; i<iNElectronsSize; i++) data << fNElectrons[i] << endl;
		makeResultFile();
		data.close();
	}
	else{
		fSignal.clear();
		fCharges.clear();
		fSignal.push_back(-1);
		fCharges.push_back(-1);
	}
}

void TAvalanche::checkDetectorGrid(){
	for(int i=0; i<iNstep; i++) assert(fElecDetectorGrid[i] == 0);
}

void TAvalanche::computeInducedSignal(){
	// drift velocity in cm/ns
	double e0 = GSL_CONST_MKSA_ELECTRON_CHARGE;//1.60217657e-19; //Coulombs
	double eps = 10.;
	double glassThickness = fResistiveLayersWidth[0];//0.2; //cm
	double weightingField = eps/(2*glassThickness + fGapWidth*eps);
	
	// Flushing signal vector
	fSignal.clear();
	
	ofstream data("out/induced_signal.dat", ios::out | ios::trunc);
	for(int i=0; i < iNElectronsSize; i++){
		fSignal.push_back(weightingField * fVx[0]*1e9 * e0 * fNElectrons[i]);
		data << fSignal.back() << endl;
	}
	data.close();
}

void TAvalanche::computeInducedCharges(){
	fCharges = vector<double>(1,fSignal[0]*fDt*1.e-9);
	ofstream data("out/induced_charges.dat", ios::out | ios::trunc);
	for(uint i=1; i<fSignal.size(); i++){
		fCharges.push_back(fCharges.back() + fSignal[i] * fDt*1.e-9);
		data << fCharges.back() << endl;
	}
	data.close();
}




double TAvalanche::n_moy(const double& x){
	return exp((fAlpha[iCurrentDetectorStep]-fEta[iCurrentDetectorStep])*x);
}

double TAvalanche::electron_multiplication2(const double& x, const double& s){
	//double s = r.Rndm();
	double nm = n_moy(x);
	double alpha = fAlpha[iCurrentDetectorStep];
	double eta = fEta[iCurrentDetectorStep];
	double k = fEta[iCurrentDetectorStep]/fAlpha[iCurrentDetectorStep];
	double thr;
	
	if (alpha > eta){
		thr = k * (nm-1)/(nm-k);
		if (s<=thr)
			return 0;
		else
		{
			if (nm > 1.e5){
				double val = 1 + trunc( log((nm-k)*(1-s)  / (nm * (1-k))) / (-( (1-k)/(nm-k) + 0.5*pow((1-k)/(nm-k),2) + 0.5*pow((1-k)/(nm-k),3)) ) );
				return val;
			}
			else{
				double val = 1 + trunc(log( (nm-k)*(1-s)  / (nm * (1-k)) ) / log( 1 - (1-k)/(nm-k) ));
				if(val<=0){
					cout << "k: " << k << endl;
					cout << "nm: " << nm << endl;
					cout << "s: " << s << endl;
					cout << "val: " << val << endl;
				}
				assert(val>0);
				return val;
			}
		}	
	}
	
	else if (alpha == eta){
		thr = alpha*x / (1.+alpha*x);
		if (s<=thr)	return 0;
		else{
			double val = 1 + trunc( log( (1.-s)*(1.+alpha*x) ) / log(thr)  );
			assert(val>0);
			return val;
		}
	}
	
	else { // alpha < eta
		thr = exp(-eta*x);
		if (s >= thr)	return 0;
		else return 1;
	}
}

double TAvalanche::multiplication(const double& n){
	double nProduced = 0;
	
	if(n > fThrCLT){
		double c = CLT(fDx,n);
		if(c < 0){
			cout << c << " " << n << endl;
			cin.ignore();
		}
		return c;
	}
	
	for(int i=0; i<n; i++){
		double s = fRandRng->RandU01();
		if (s==1)	s = fRandRng->RandU01();
		nProduced += electron_multiplication2(fDx,s);
	}
	return nProduced;
}

double TAvalanche::CLT(const double& x, const double& n){
	double nm = n_moy(x);
	if(fAlpha[iCurrentDetectorStep] == 0) 	fAlpha[iCurrentDetectorStep] += numeric_limits<double>::epsilon();
	double k = fEta[iCurrentDetectorStep]/fAlpha[iCurrentDetectorStep];
	double alpha = fAlpha[iCurrentDetectorStep];
	double eta = fEta[iCurrentDetectorStep];
	double m = n*nm;
	
	double sigma;
	if (alpha > eta)	sigma = sqrt(n) * sqrt( ((1+k)/(1-k)) * nm * (nm-1) );
	else if(alpha == eta)	sigma = sqrt(n) * sqrt(2*alpha*x);
	else	sigma = sqrt(n) * sqrt( exp(-2*eta*x)*(exp(eta*x)-1) ); // alpha < eta
	
	double c = Gaus(m, sigma, fRandRngCLT);
	
	if (c > 1e10){
		ofstream dump("out/dump1.dat", ios::out | ios::trunc);
		cout << "k: " << k << endl;
		cout << "eta: " << fEta[iCurrentDetectorStep] << endl;
		cout << "alpha: " << fAlpha[iCurrentDetectorStep] << endl;
		cout << "E: " << fE[iCurrentDetectorStep] << endl;
		cout << "nm: " << nm << endl;
		cout << "m: " << m << endl;
		cout << "sigma: " << sigma << endl;
		cout << "c: " << c << endl;
		for(int i=0; i<iNstep; i++)
			dump << fE[i] << "\t" << fAlpha[i] << endl;
		dump.close();
		cin.ignore();
	}
	
	if(abs(c > 1) and !(trunc(c)>0)){
		ofstream dump("out/dump2.dat", ios::out | ios::trunc);
		cout << "k: " << k << endl;
		cout << "eta: " << fEta[iCurrentDetectorStep] << endl;
		cout << "alpha: " << fAlpha[iCurrentDetectorStep] << endl;
		cout << "E: " << fE[iCurrentDetectorStep] << endl;
		cout << "nm: " << nm << endl;
		cout << "m: " << m << endl;
		cout << "sigma: " << sigma << endl;
		cout << "c: " << c << endl;
		for(int i=0; i<iNstep; i++)
			dump << fE[i] << "\t" << fAlpha[i] << endl;
		dump.close();
	}
	
	if(abs(c < 1) and !(round(c)>=0)){
		ofstream dump("out/dump3.dat", ios::out | ios::trunc);
		cout << "k: " << k << endl;
		cout << "eta: " << fEta[iCurrentDetectorStep] << endl;
		cout << "alpha: " << fAlpha[iCurrentDetectorStep] << endl;
		cout << "E: " << fE[iCurrentDetectorStep] << endl;
		cout << "nm: " << nm << endl;
		cout << "m: " << m << endl;
		cout << "sigma: " << sigma << endl;
		cout << "c: " << c << endl;
		for(int i=0; i<iNstep; i++)
			dump << fE[i] << "\t" << fAlpha[i] << endl;
		dump.close();
	}
	
	if (abs(c>1)){
		assert(trunc(c)>0);
		return trunc(c);
	}
	else{
		assert(round(c) >= 0);
		return round(c);
	}
}

void TAvalanche::computeLongitudinalDiffusion(){
	bool computed = true;
	
	double newDetectorGrid[iNstep];
	for(int i=0; i<iNstep; i++) 	newDetectorGrid[i] = 0;
	
	//double nIni = sumArray(fElecDetectorGrid, iNstep);
	
	if(bFullLongiDiff){
		double pos, newPos;
		for(int iz=0; iz<iNstep; iz++){
			for(int n=0; n<fElecDetectorGrid[iz]; n++){
				pos = (iz+1) * fDx;
				//newPos = gsl_ran_gaussian(fGSLRng, fLongiDiffSigma) + pos;
				newPos = Gaus(pos, fLongiDiffSigma, fRandRngLongiDiff);
				newDetectorGrid[(int)trunc(newPos/fDx)]++;
				//cout << pos << " " << newPos << " " << (int)trunc(newPos/fDx) << " " << iz << " " << fElecDetectorGrid[iz] << " " << newDetectorGrid[(int)trunc(newPos/fDx)] << endl;
				//cin.ignore();
			}
		}
	}
	
	else{	/*	FIXME: DOES NOT WORK	*/
		for(int i=0; i<iNstep; i++)	newDetectorGrid[i] = 0;
		for(int i=0; i<3; i++) newDetectorGrid[i] = fElecDetectorGrid[i];
		for(int i=iNstep; i<iNstep-3; i--) newDetectorGrid[i] = fElecDetectorGrid[i];
		
		for(int iStep=0; iStep<iNstep; iStep++){
		
		if(iStep > 3 and iStep < iNstep-3){
			if(!computed)	computed = true;
				double n0 = fElecDetectorGrid[iStep];
				double n[7];
				double intPart;
				for(int i=0; i<7; i++){
					if ( modf(fLongiDiffFrac[i]*n0, &intPart) >= 0.5 ) n[i] = ceil(fLongiDiffFrac[i]*n0);
					else 	n[i] = floor(fLongiDiffFrac[i]*n0);
				}
				newDetectorGrid[iStep] += n[0];
				newDetectorGrid[iStep+1] += n[1];
				newDetectorGrid[iStep-1] += n[2];
				newDetectorGrid[iStep+2] += n[3];
				newDetectorGrid[iStep-2] += n[4];
				newDetectorGrid[iStep+3] += n[5];
				newDetectorGrid[iStep-3] += n[6];
			
				// The electrons "lost" by computing the fractions are put in the center bin
				newDetectorGrid[iStep] += n0 - sumArray(n,7);
				if(modf(sumArray(n,7), &intPart) != 0){
					for(int i=0; i<7; i++)	cout << n[i] << endl;
					cin.ignore();
				}
			}
		}
	
		//double nFin = sumArray(newDetectorGrid, iNstep);

		//cout << nIni << " " << nFin << " " << computed << " " << iDebug << " " << fNElectrons[iDebug] << endl;
		iDebug++;
		double intPart;
		if(modf(sumArray(newDetectorGrid, iNstep), &intPart) > 0.)	cin.ignore();
		//if(computed)	
	}
	
	memcpy(fElecDetectorGrid, newDetectorGrid, sizeof(double)*iNstep); //fElecDetectorGrid = newDetectorGrid;
}

bool TAvalanche::avalanche(){
	bool bSpaceChargeLimit = false;
	iTimeStep = 1;
	iDebug = 1;
	
	double nAnode = 0;
	
	if(bTestPot){
		cout << "test Pot" << endl;
		//PlotPot();
		//PlotField();
		//PlotSimplifiedField();
		return true;
	}
		
	while(true){
		iTimeStep++;
		computeSCEffect();
		//if(fNElectrons[iTimeStep-1] > fSpaceChargeLimit and !bSpaceChargeLimit)	bSpaceChargeLimit = true;
		
		if(bPrintDetectorGrid) {
			printDetectorGrid();
			cin.ignore();
			cout << "\033[2J\033[1;1H";
		}
		double copy[iNstep];
		memcpy(copy, fElecDetectorGrid, sizeof(double)*iNstep);
		
		for(iCurrentDetectorStep=0; iCurrentDetectorStep<iNstep-1; iCurrentDetectorStep++){
			double n = copy[iCurrentDetectorStep];
			double nProduced;
			//if(bSpaceChargeLimit)	{
				//nProduced = n;
				////fElecDetectorGrid[i] = 0;
			//}
			//else	
			nProduced = multiplication(n);
			
			fElecDetectorGrid[iCurrentDetectorStep+1] = nProduced;
			
			if( (nProduced - copy[iCurrentDetectorStep]) > 0 )
				fPosIonDetectorGrid[iCurrentDetectorStep+1] += nProduced - copy[iCurrentDetectorStep];
			else
				fNegIonDetectorGrid[iCurrentDetectorStep+1] += copy[iCurrentDetectorStep] - nProduced;
			
			if(iCurrentDetectorStep==0)	fElecDetectorGrid[iCurrentDetectorStep] = 0;
			fNElectrons[iTimeStep] += nProduced;
		}
		
		computeLongitudinalDiffusion();
		
		//fElecDetectorGrid[iNstep] = 0; // Electrons in the last detector step will leave the gas gap
		nAnode += fElecDetectorGrid[iNstep-1];
		if(fNElectrons[iTimeStep] == 0)	break;
		
		if(iTimeStep % 100 == 0)	makeSnapshot();
		
		cout << "time step: " << iTimeStep << "\t Nelec: " << fNElectrons[iTimeStep] << "\t" << "NelecLastBin: " << nAnode << " ";
		cout << " " << -sumArray(fPosIonDetectorGrid, iNstep)+sumArray(fElecDetectorGrid, iNstep)+sumArray(fNegIonDetectorGrid, iNstep) << endl;
		
		if( iTimeStep > iNElectronsSize )	return false;
	}
	//cout << sumArray(fNElectrons,iNstep) << " " << nAnode << endl;
	return true;
}

void TAvalanche::computeSCEffect(){
	double cm = 0.01;
	double SCEField[iNstep];
	for(int z=0; z<iNstep; z++){
		double tmp=0;
		for(int zp=0; zp<iNstep; zp++){
			tmp += (-fElecDetectorGrid[zp]) * interpolateEbar((z+1)*fDx*cm, (zp+1)*fDx*cm, iTimeStep*fDx);
		}
		SCEField[z] = tmp;
	}
	
	for(int z=0; z<iNstep; z++){
		fE[z] = fEini + SCEField[z];
		double* transportParams = fDet->getTransportParameters(fE[z],0.,0.);
		fAlpha[z] = transportParams[0];
		fEta[z] = transportParams[1];
		fVx[z] = transportParams[2];
	}
}

void TAvalanche::printDetectorGrid() const{

	for(int i=0; i<iNstep; i++)
		cout << fElecDetectorGrid[i] << " ";
	cout << endl << endl;
}

void TAvalanche::makeSnapshot(){
	cout << "Snapshot at: " << iTimeStep*fDt << endl;
	ofstream data(("out/snap-"+to_string(round(iTimeStep*fDt))+".dat").c_str(), ios::out | ios::trunc);
	for(int i=0; i<iNstep; i++){
		data << i*fDx << "\t" << fElecDetectorGrid[i] << "\t" << fPosIonDetectorGrid[i] << "\t" << fNegIonDetectorGrid[i] << "\t" << fE[i] << "\t" << fAlpha[i] << endl;
	}
	data.close();
}

void TAvalanche::testInterpolation(){
	cout << interpolateEbar(4.e-5,4.e-5,0.006) << endl;
}

size_t TAvalanche::getIndex3D(const size_t& i, const size_t& j, const size_t& k){
	int n = iEbarTableSize+1;
	return (long)i*(long)n*(long)n + (long)j*(long)n + (long)k;
}

double TAvalanche::interpolateEbar(const double& z, const double& zp, const double& l){  //x==z y==zp z==l
	//cout << z << " " << zp << " " << l << endl;
	//cout << fEbarZarray[iEbarTableSize-1] << " " << fEbarZparray[iEbarTableSize-1] << " " << fEbarLarray[iEbarTableSize-1] << endl;
	//for (int i=0; i<iEbarTableSize; i++)	cout << fEbarZarray[i] << " ";
	//cout << endl;
	/*
	assert(z >= fEbarZarray[0]);
	assert(zp >= fEbarZparray[0]);
	assert(l >= fEbarLarray[0]);
	assert(z <= fEbarZarray[iEbarTableSize-1]);
	assert(zp <= fEbarZparray[iEbarTableSize-1]);
	assert(l <= fEbarLarray[iEbarTableSize-1]);
	*/
	size_t iz0 = getLowerIndex(fEbarZarray,fEbarZarray+iEbarTableSize, z);
	size_t izp0 = getLowerIndex(fEbarZparray,fEbarZparray+iEbarTableSize, zp);
	size_t il0 = getLowerIndex(fEbarLarray,fEbarLarray+iEbarTableSize, l);
	
	//cout << iz0 << " " << izp0 << " " << il0 << endl;
	
	size_t iz1 = iz0+1;
	size_t izp1 = izp0+1;
	size_t il1 = il0+1;
	
	double z0 = fEbarZarray[iz0], zp0 = fEbarZparray[izp0], l0 = fEbarLarray[il0];
	double z1 = fEbarZarray[iz1], zp1 = fEbarZparray[izp1], l1 = fEbarLarray[il1];
	
	//cout << z0 << " " << zp0 << " " << l0 << endl;
	//cout << z1 << " " << zp1 << " " << l1 << endl;
	
	double zd = (z - z0)/(z1 - z0);
	double zpd = (zp - zp0)/(zp1 - zp0);
	double ld = (l - l0)/(l1 - l0);
	
	//interpolate along z
	double c00 = fEbarTable[getIndex3D(iz0,izp0,il0)]*(1-zd) + fEbarTable[getIndex3D(iz1,izp0,il0)]*zd;
	double c10 = fEbarTable[getIndex3D(iz0,izp1,il0)]*(1-zd) + fEbarTable[getIndex3D(iz1,izp1,il0)]*zd;
	double c01 = fEbarTable[getIndex3D(iz0,izp0,il1)]*(1-zd) + fEbarTable[getIndex3D(iz1,izp0,il1)]*zd;
	double c11 = fEbarTable[getIndex3D(iz0,izp1,il1)]*(1-zd) + fEbarTable[getIndex3D(iz1,izp1,il1)]*zd;
	
	//interpolate along zp
	double c0 = c00*(1-zpd) + c10*zpd;
	double c1 = c01*(1-zpd) + c11*zpd;
	
	//interpolate along l
	double c = c0*(1-ld) + c1*ld;
	
	return c;
}
