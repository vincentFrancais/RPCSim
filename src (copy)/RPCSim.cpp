#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <assert.h> 
#include <utility>

#include "RPCSim.hpp"
#include "clusters.hpp"
#include "helper_functions.hpp"

#define PI 3.14159265358979323846


using namespace std;


RPCSim::RPCSim(const int& nStep, const double& gasGap){
	//fRand = TRandom3(0);
	//fRandCLT = TRandom3(456);
	//fRandLongiDiff = TRandom3(758);
	
	iNstep = nStep;
	fDetectorGrid = new double[iNstep];
	//fNElectrons = vector<double>(0,2*iNstep);
	
	fThrCLT = 1.e4;
	fSpaceChargeLimit = 5.e7;
	fLongiDiffComputeLimit = 500;
	
	fGapWidth = gasGap;
	
	iDebug = 0;
	
	//for(int i=0; i<100; i++)	cout << fRandRng.RandU01() << endl;
	//cin.ignore();
}

void RPCSim::setGasMixture(Garfield::MediumMagboltz* gas){
	mGas = gas;
	int nComponents = mGas->GetNumberOfComponents();
	vector< pair<string,double> > composition;
	for(int i=0; i<nComponents; i++) {
		double fraction;
		string label;
		mGas->GetComponent(i,label,fraction);
		composition.push_back( make_pair(label,fraction) );
	}
	
	for(vector< pair<string,double> >::iterator it = composition.begin(); it != composition.end(); it++)	mGasTableName += it->first + "-" + to_string(it->second) + "_";
	mGasTableName += "temp-" + to_string(mGas->GetTemperature()) + "_pres-" + to_string(mGas->GetPressure()) + ".gas";
	
	if(!file_exist("gastables/"+mGasTableName))	makeGasTable();
	mGas->LoadGasFile("gastables/"+mGasTableName);
	//gas->EnableDebugging();
	
	fTemperature = mGas->GetTemperature();
	fPressure = mGas->GetPressure();
	
}

void RPCSim::setElectricField(const double& Ex, const double& Ey, const double& Ez){
	fElectricField[0] = Ex;
	fElectricField[1] = Ey;
	fElectricField[2] = Ez;
}

Garfield::MediumMagboltz* RPCSim::getGas(){
	return mGas;
}

void RPCSim::makeGasTable(){
	cout << "Generating gas table -- could take some time" << endl;
	Garfield::MediumMagboltz* gas = mGas;
	gas->SetFieldGrid(100., 120000., 70, false, 0., 0., 1, 0., 0., 1);
	gas->EnableDebugging();
	gas->Initialise();  
	gas->DisableDebugging();
	//const double lP = 0.0;
	//const double rP = 0.0;
	//gas->EnablePenningTransfer(rP, lP, "c2h2f4");
	gas->GenerateGasTable(4, true);
	
	gas->WriteGasFile("gastables/"+mGasTableName);
}

void RPCSim::initialiseDetector(){
	
	fDx = fGapWidth/iNstep;
	double cx = fGapWidth/2., cy = 0, cz = 0, lx = fGapWidth/2., ly = 10, lz = 10;
	// Detector geometry
	// Gap [cm]
	Garfield::SolidBox* box = new Garfield::SolidBox(cx, cy, cz, lx, ly, lz);
	Garfield::GeometrySimple* geo = new Garfield::GeometrySimple();
	geo->AddSolid(box, mGas);

	// Make a component
	Garfield::ComponentConstant* comp = new Garfield::ComponentConstant();
	comp->SetGeometry(geo);
	comp->SetElectricField(fElectricField[0], fElectricField[1], fElectricField[2]);
	
	// Make a sensor
	mSensor = new Garfield::Sensor();
	mSensor->AddComponent(comp);
	
	// Get transport parameters
	mGas->ElectronTownsend(fElectricField[0],fElectricField[1],fElectricField[2],0.,0.,0., fAlpha);
	mGas->ElectronAttachment(fElectricField[0],fElectricField[1],fElectricField[2],0.,0.,0., fEta);
	mGas->ElectronVelocity(fElectricField[0],fElectricField[1],fElectricField[2],0.,0.,0.,fVx,fVy,fVz);
	mGas->ElectronDiffusion(fElectricField[0],fElectricField[1],fElectricField[2],0.,0.,0.,fDiffL,fDiffT);
	
	// Multiply by -1 in order to have the electron drift velocity going forward on the x-axis
	fVx *= -1., fVy *= -1, fVz *= -1;

	k = fEta/fAlpha;

	fDt = fDx/fVx; //ns
	
	assert( (fAlpha>fEta) or (fAlpha>0) );
	
	cout << "Transport parameters:" << endl;
	cout << "\talpha: " << fAlpha << " cm-1" << endl;
	cout << "\teta: " << fEta << " cm-1" << endl;
	cout << "\tE: " << fElectricField[0] << " kV/cm" << endl;
	cout << "\tDrift velocity: (" << fVx << "," << fVy << "," << fVz << ")" << endl;
	cout << "\tDiffusion coefficient: (" << fDiffL << ", " << fDiffT << ")" << endl;
	cout << "\tDt: " << fDt << endl;
}

void RPCSim::setTrackHeed(const string& particleName, const double& momentum, const double& x0, const double& theta){
	mTrack = new Garfield::TrackHeed();
	mTrack->SetSensor(mSensor);
	mTrack->SetParticle(particleName);
	mTrack->SetMomentum(momentum);
	
	double t0 = 0.;
	double y0 = 0, z0 = 0; 
	double dx0 = cos(theta * PI / 180.0), dy0 = sin(theta * PI / 180.0), dz0 = 0.;
	
	mTrack->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
	
	double xc = 0., yc = 0., zc = 0., tc = 0.; 	// Cluster coordinates
	int nc = 0; 								// Number of electrons produced in a collision
	double ec = 0.; 							// Energy loss in a collision
	double dummy = 0.; 							// Dummy variable (not used at present)
	fNElectrons = vector<double>(1.2*iNstep,0);
	//double detectorGrid[Nstep];
	for(int i=0; i<iNstep; i++)	fDetectorGrid[i] = 0;
	
	while (mTrack->GetCluster(xc, yc, zc, tc, nc, ec, dummy)){
		//clusters.addCluster(int(trunc(xc/fDx)), nc);
		fNElectrons[0] += nc;
		fDetectorGrid[int(trunc(xc/fDx))] += nc;
	}
	
	delete mTrack;
}

void RPCSim::simulateEvent(){
	//Garfield::TrackHeed* track = new Garfield::TrackHeed();
	//track->SetSensor(mSensor);
	//track->SetParticle(particleName);
	//track->SetMomentum(momentum);
	
	//double t0 = 0.;
	//double dx0 = cos(angle * PI / 180.0), dy0 = sin(angle * PI / 180.0), dz0 = 0.;
	
	//mTrack->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
	
	////Clusters clusters = Clusters();
	//double xc = 0., yc = 0., zc = 0., tc = 0.; 	// Cluster coordinates
	//int nc = 0; 								// Number of electrons produced in a collision
	//double ec = 0.; 							// Energy loss in a collision
	//double dummy = 0.; 							// Dummy variable (not used at present)
	//fNElectrons = vector<double>(1.2*iNstep,0);
	////double detectorGrid[Nstep];
	//for(int i=0; i<iNstep; i++)	fDetectorGrid[i] = 0;
	
	//while (mTrack->GetCluster(xc, yc, zc, tc, nc, ec, dummy)){
		////clusters.addCluster(int(trunc(xc/fDx)), nc);
		//fNElectrons[0] += nc;
		//fDetectorGrid[int(trunc(xc/fDx))] += nc;
	//}
	
	//cout << "Clusters: " << clusters.size() << endl;
	//cout << "Electrons: " << fNElectrons[0] << endl;
	//cout << "space charge limit: " << fSpaceChargeLimit << endl;
	
	avalanche();
	checkDetectorGrid();
	//cout << "====================" << endl;
	//for(uint i=0; i<fNElectrons.size(); i++)	cout << fNElectrons[i] << endl;
	computeInducedSignal();
	computeInducedCharges();
	ofstream data("out/electrons.dat", ios::out | ios::trunc);
	for(uint i=0; i<fNElectrons.size(); i++) data << fNElectrons[i] << endl;
	data.close();
}

void RPCSim::checkDetectorGrid(){
	for(int i=0; i<iNstep; i++) assert(fDetectorGrid[i] == 0);
}

void RPCSim::computeInducedSignal(){
	// drift velocity in cm/ns
	double e0 = 1.60217657e-19; //Coulombs
	double eps = 10.;
	double glassThickness = 0.2; //cm
	double weightingField = eps/(2*glassThickness + fGapWidth*eps);
	
	// Flushing signal vector
	fSignal.clear();
	
	ofstream data("out/induced_signal.dat", ios::out | ios::trunc);
	for(uint i=0; i < fNElectrons.size(); i++){
		fSignal.push_back(weightingField * fVx*1e9 * e0 * fNElectrons[i]);
		data << fSignal.back() << endl;
	}
	data.close();
}

void RPCSim::computeInducedCharges(){
	fCharges = vector<double>(1,fSignal[0]*fDt*1.e-9);
	ofstream data("out/induced_charges.dat", ios::out | ios::trunc);
	for(uint i=1; i<fSignal.size(); i++){
		fCharges.push_back(fCharges.back() + fSignal[i] * fDt*1.e-9);
		data << fCharges.back() << endl;
	}
	data.close();
}

void RPCSim::computeChargeSpectra(const double& x0, const double& angle, const int& nEvents){
	ofstream data("out/chargeSpectra.dat", ios::out | ios::trunc);
	
	//Garfield::TrackHeed* track = new Garfield::TrackHeed();
	//track->SetSensor(mSensor);
	//track->SetParticle(particleName);
	//track->SetMomentum(momentum);
	
	for(int j=0; j<nEvents; j++){
		if (j % 1 == 0) std::cout << j << "/" << nEvents << "\n";
		double y0 = 0., z0 = 0., t0 = 0.;
		//double theta = 0; //degree
		double dx0 = cos(angle * PI / 180.0), dy0 = sin(angle * PI / 180.0), dz0 = 0.;
		
		mTrack->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
		
		double xc = 0., yc = 0., zc = 0., tc = 0.; 	// Cluster coordinates
		int nc = 0; 								// Number of electrons produced in a collision
		double ec = 0.; 							// Energy loss in a collision
		double dummy = 0.; 							// Dummy variable (not used at present)
		fNElectrons = vector<double>(iNstep,0);
		fDetectorGrid = new double[iNstep];
		for(int i=0; i<iNstep; i++)	fDetectorGrid[i] = 0;
	
		while (mTrack->GetCluster(xc, yc, zc, tc, nc, ec, dummy)){
			fNElectrons[0] += nc;
			fDetectorGrid[int(trunc(xc/fDx))] += nc;
		}
	
		//cout << "Clusters: " << clusters.size() << endl;
		//cout << "Electrons: " << fNElectrons[0] << endl;
		//cout << "space charge limit: " << fSpaceChargeLimit << endl;
	
		avalanche();
		computeInducedSignal();
		computeInducedCharges();
		
		data << fCharges.back() << endl;
	}
	data.close();

}







double RPCSim::n_moy(const double& x){
	return exp((fAlpha-fEta)*x);
}

double RPCSim::electron_multiplication2(const double& x, const double& s){
	//double s = r.Rndm();
	double nm = n_moy(x);
	double thr = k * (nm-1)/(nm-k);
	
	if (s<=thr)
		return 0;
	else
	{
		if (nm > 1.e5)
		{
			double val = 1 + trunc( log((nm-k)*(1-s)  / (nm * (1-k))) / (-( (1-k)/(nm-k) + 0.5*pow((1-k)/(nm-k),2) + 0.5*pow((1-k)/(nm-k),3)) ) );
			return val;
		}
		else
		{
			double val = 1 + trunc(log( (nm-k)*(1-s)  / (nm * (1-k)) ) / log( 1 - (1-k)/(nm-k) ));
			assert(val>0);
			return val;
		}
	}
}

double RPCSim::multiplication(const double& n){
	double nProduced = 0;
	if(n > fThrCLT){
		double c = RPCSim::CLT(fDx,n);
		if(c < 0){
			cout << c << " " << n << endl;
			cin.ignore();
		}
		return c;
	}
	for(int i=0; i<n; i++){
		double s = fRandRng.RandU01();
		if (s==1)	s = fRandRng.RandU01();
		nProduced += RPCSim::electron_multiplication2(fDx,s);
	}
	return nProduced;
}

double RPCSim::CLT(const double& x, const double& n){
	//TRandom3 r(456);
	double nm = n_moy(x);
	double m = n*nm;
	double sigma = sqrt(n) * sqrt( ((1+k)/(1-k)) * nm * (nm-1) );
	double c = generateGaussianNumber(m, sigma, fRandRngCLT);	//fRand.Gaus(m,sigma);
	assert(trunc(c)>0);
	return trunc(c);
}

void RPCSim::computeLongitudinalDiffusion(){
	bool debug = false;
	double sigma = fDiffL*sqrt(fDx);
	
	double newDetectorGrid[iNstep];
	for(int i=0; i<iNstep; i++)	newDetectorGrid[i] = 0;
	
	double nIni = sumArray(fDetectorGrid, iNstep);
	
	double frac[7];
	frac[0] = gaussIntegral(200, -0.5*fDx, 0.5*fDx, 0, sigma);
	frac[1] = gaussIntegral(200, 0.5*fDx, 1.5*fDx, 0, sigma);
	frac[2] = gaussIntegral(200, -1.5*fDx, -0.5*fDx, 0, sigma);
	frac[3] = gaussIntegral(200, 1.5*fDx, 2.5*fDx, 0, sigma);
	frac[4] = gaussIntegral(200, -2.5*fDx, -1.5*fDx, 0, sigma);
	frac[5] = gaussIntegral(200, 2.5*fDx, 3.5*fDx, 0, sigma);
	frac[6] = gaussIntegral(200, -3.5*fDx, -2.5*fDx, 0, sigma);
	
	for(int iStep=0; iStep<iNstep; iStep++){
		if(fDetectorGrid[iStep] < fLongiDiffComputeLimit or iStep < 10 or iStep > iNstep-10){
			for(int i=0; i<fDetectorGrid[iStep]; i++){
				double pos = (iStep+1)*fDx - 0.5*fDx;
				double newPos = generateGaussianNumber(pos,sigma,fRandRngLongiDiff);//fRandLongiDiff.Gaus(pos,sigma);
				if(newPos<=fGapWidth or newPos>=0)	newDetectorGrid[int(trunc(newPos/fDx))] += 1;
			}
		}
		else{
			debug = true;
			double n0 = fDetectorGrid[iStep];
			double n[7];
			double intPart;
			for(int i=0; i<7; i++){
				if ( modf(frac[i]*n0, &intPart) >= 0.5 ) n[i] = ceil(frac[i]*n0);
				else 	n[i] = floor(frac[i]*n0);
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
	double nFin = sumArray(newDetectorGrid, iNstep);
	bool b;
	if(nIni == nFin) b = true;
	else 	b = false;
	//cout << nIni << " " << nFin << " " << b << " here " << debug << " " << iDebug << " " << fNElectrons[iDebug+1] << endl;
	iDebug++;
	double intPart;
	if(modf(sumArray(newDetectorGrid, iNstep), &intPart) > 0.)	cin.ignore();
	memcpy(fDetectorGrid, newDetectorGrid, sizeof(double)*iNstep); //fDetectorGrid = newDetectorGrid;
}

bool RPCSim::avalanche(){
	bool bSpaceChargeLimit = false;
	int timeStep = 1;
	while(true){
		if(fNElectrons[timeStep-1] > fSpaceChargeLimit and !bSpaceChargeLimit)	bSpaceChargeLimit = true;
		double copy[iNstep];
		memcpy(copy, fDetectorGrid, sizeof(double)*iNstep);
		for(int i=0; i<iNstep-1; i++){
			double n = copy[i];
			double nProduced;
			if(bSpaceChargeLimit)	{
				nProduced = n;
				//fDetectorGrid[i] = 0;
			}
			else	nProduced = multiplication(n);
			fDetectorGrid[i+1] = nProduced;
			fNElectrons[timeStep] += nProduced;
		}
		computeLongitudinalDiffusion();
		fDetectorGrid[iNstep] = 0; // Electrons in the last detector step will leave the gas gap
		
		if(fNElectrons[timeStep] == 0)	break;

		timeStep++;
		assert(timeStep <= static_cast<int> (fNElectrons.size()) );
		//if(timeStep > fNElectrons.size())	fNElectrons.push_back(0.);
	}

	return true;
}

void RPCSim::singleAvalanche(const double& x0){
	
	Clusters clusters = Clusters();
	fNElectrons = vector<double>(iNstep,0);
	for(int i=0; i<iNstep; i++)	fDetectorGrid[i] = 0;
	
	clusters.addCluster(int(trunc(x0/fDx)), 1);
	fNElectrons[0] = 1;
	fDetectorGrid[int(trunc(x0/fDx))] = 1;
	
	//cout << "Clusters: " << clusters.size() << endl;
	//cout << "Electrons: " << fNElectrons[0] << endl;
	//cout << "space charge limit: " << fSpaceChargeLimit << endl;
	
	avalanche();
	computeInducedSignal();
	computeInducedCharges();
	ofstream data("out/electrons.dat", ios::out | ios::trunc);
	for(int i=0; i<iNstep; i++) data << fNElectrons[i] << endl;
	data.close();
}
