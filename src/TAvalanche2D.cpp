
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <limits>
#include <iterator>
#include <algorithm>

#include "TAvalanche2D.hpp"

using namespace std;


TAvalanche2D::TAvalanche2D(TDetector* det, TConfig& config, sfmt_t sfmt, const int& id) : TAvalanche() {
	fDetector = det;
	fConfig = config;
	
	fRngMult = new TRandomEngineMTDC(id,1234,4321);
	fRngCLT = new TRandomEngineMTDC(id+1,1234,4321);
	fRngLongiDiff = new TRandomEngineMTDC(id+2,1234,4321);
}

inline double TAvalanche2D::indexDetectorGrid(const int& r, const int& z) {
	return r * fNstepr + z;
}

double TAvalanche2D::electricFieldStrength (const int& r, const int&z) {
	return sqrt ( fEz.at(indexDetectorGrid(r,z))*fEz.at(indexDetectorGrid(r,z)) + fEr.at(indexDetectorGrid(r,z))*fEr.at(indexDetectorGrid(r,z)) );
}

inline double TAvalanche2D::n_moy(const double& x) {
	return exp( ( fAlpha.at(indexDetectorGrid(fCurrentStepR,fCurrentStepZ))-fEta.at(indexDetectorGrid(fCurrentStepR,fCurrentStepZ)) )*x );
}

void TAvalanche2D::init() {
	
	fGridSize = fNstepr*fNstepz;
	
	fElecDetectorGrid = vector<double> (fGridSize, 0);
	fPosIonDetectorGrid = vector<double> (fGridSize, 0);
	fNegIonDetectorGrid = vector<double> (fGridSize, 0);
	fEz = vector<double> (fGridSize, 0);
	fEr = vector<double> (fGridSize, 0);
	
	/* At first only the Z component of the electric field is non null */
	for (int i=0; i<fGridSize; i++)
		fEz.at(i) = fConfig.ElectricField;
		//fEr.at(i) = fConfig.ElectricField;
	
}

void TAvalanche2D::updateParameters() {
	vector<double> params = fDetector->getTransportParameters( electricFieldStrength(fCurrentStepR,fCurrentStepZ), 0. , 0. );
	fCurrentAlpha = params.at(0);
	fCurrentEta = params.at(1);
	fCurrentVel = params.at(3);
}

double TAvalanche2D::multiplicationCLT(const double& x, const double& n){
	double nm = n_moy(x);
	double alpha = fAlpha.at(indexDetectorGrid(fCurrentStepR,fCurrentStepZ));
	double eta = fEta.at(indexDetectorGrid(fCurrentStepR,fCurrentStepZ));
	double k = eta/alpha;
	double m = n*nm;
	
	double sigma;
	if (alpha == eta)	sigma = sqrt(n) * sqrt(2*alpha*x);
	else if (alpha < 0.01)	sigma = sqrt(n) * sqrt( exp(-2*eta*x)*(exp(eta*x)-1) ); // alpha << 0
	else	sigma = sqrt(n) * sqrt( ((1+k)/(1-k)) * nm * (nm-1) ); // alpha, eta > 0
	
	double c = Gaus(m, sigma, fRngCLT);
	
	if (c > 1e10){
		if ( fVerbose ){
			cout << "k: " << k << endl;
			cout << "eta: " << eta << endl;
			cout << "alpha: " << alpha << endl;
			cout << "E: " << electricFieldStrength(fCurrentStepR,fCurrentStepZ) << endl;
			cout << "nm: " << nm << endl;
			cout << "m: " << m << endl;
			cout << "sigma: " << sigma << endl;
			cout << "c: " << c << endl;
		}
		//cerr << "Explosive behavior detected -- stopping avalanche" <<endl;
		fAvalStatus = AVAL_EXPLOSIVE_BEHAVIOR_CLT;
		return 0;
		//cin.ignore();
	}
	
	if( ( (abs(c > 1) and !trunc(c)>0) or (abs(c < 1) and !round(c)>=0) ) and fVerbose ) {
		cout << "k: " << k << endl;
		cout << "eta: " << eta << endl;
		cout << "alpha: " << alpha << endl;
		cout << "E: " << electricFieldStrength(fCurrentStepR,fCurrentStepZ) << endl;
		cout << "nm: " << nm << endl;
		cout << "m: " << m << endl;
		cout << "sigma: " << sigma << endl;
		cout << "c: " << c << endl;
	}
	
	if (abs(c>1)){
		if( !(trunc(c)>0) )
			fAvalStatus = AVAL_CLT_FAIL;
		return trunc(c);
	}
	
	else{
		if ( !(round(c) >= 0) ) 
			fAvalStatus = AVAL_CLT_FAIL;
		return round(c);
	}
}

double TAvalanche2D::multiplicationRiegler(const double& x, const double& s) {
	double nm = n_moy(x);
	double alpha = fAlpha.at(indexDetectorGrid(fCurrentStepR,fCurrentStepZ));
	double eta = fEta.at(indexDetectorGrid(fCurrentStepR,fCurrentStepZ));
	double k = eta/alpha;
	double thr;
	
	
	if (alpha == eta){
		thr = alpha*x / (1.+alpha*x);
		if (s<=thr)	return 0;
		else{
			double val = 1 + trunc( log( (1.-s)*(1.+alpha*x) ) / log(thr)  );
				if (val < 0)
					fAvalStatus = AVAL_MULT_FAIL;
			return val;
		}
	}
	
	else if (alpha < 0.01) { // alpha << 0
		thr = exp(-eta*x);
		if (s >= thr)	return 0;
		else return 1;
	}
	
	else { //alpha, eta > 0
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
				if (val < 0)
					fAvalStatus = AVAL_MULT_FAIL;
				return val;
			}
		}	
	}
	
}

double TAvalanche2D::electronMultiplication(const double& n){
	double nProduced = 0;
	double s;
	double l;
	
	if(n > fThrCLT){
		double c = multiplicationCLT(fDz,n);
		if( fAvalStatus  == AVAL_CLT_FAIL ){
			/* CLT has failed to return an acceptable value. Try with classic multiplication */
			for(int i=0; i<n; i++){
				s = fRngMult->RandU01();
				if (s==1)	s = fRngMult->RandU01();
				nProduced += multiplicationRiegler(fDz,s);
			}
			if( fAvalStatus == AVAL_CLT_FAIL )
				fAvalStatus = AVAL_NO_ERROR;
			return nProduced;
		}
		return c;
	}
	
	for(int i=0; i<n; i++){
		double s = fRngMult->RandU01();
		if (s==1)	s = fRngMult->RandU01();
		nProduced += multiplicationRiegler(fDz,s);
	}
	return nProduced;
}

bool TAvalanche2D::propagate() {
	
}

bool TAvalanche2D::avalanche() {
	fTimeStep = 1;
	
	while (true) {
		propagate();
	}
	return true;
}
