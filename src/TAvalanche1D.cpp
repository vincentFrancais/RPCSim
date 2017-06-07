
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <utility>
#include <limits>
#include <iterator>
#include <random>
#include <algorithm>

#include "TAvalanche1D.hpp"
#include "helper_functions.hpp"
#include "TConstants.hpp"

using namespace std;

TAvalanche1D::TAvalanche1D(TDetector* det, TConfig& config, sfmt_t sfmt, const int& id) : TAvalanche() {
	fDet = det;
	fConfig = config;

	fRngMult = new TRandomEngineMTDC(id,1234,fConfig.globalSeed);
	fRngCLT = new TRandomEngineMTDC(id+1,1234,fConfig.globalSeed);
	
	if (fConfig.generator == "SFMT")
		fRngLongiDiff = new TRandomEngineSFMT(sfmt);
	else
		fRngLongiDiff = new TRandomEngineMTDC(id+2,1234,fConfig.globalSeed);
		
	fRngMisc = new TRandomEngineMT(getUUID());
	
	cout << "Longitudinal diffusion generator: " << fRngLongiDiff->Generator() << endl;
	cout << endl;
	
	/* Call to initialisation function */
	init();
}

TAvalanche1D::~TAvalanche1D() {
	delete fRngLongiDiff;
	delete fRngMult;
	delete fRngCLT;
}

void TAvalanche1D::init() {
	tId = gettid();
	Id = count++;
	
	iNstep = fConfig.nSteps;
	 
	fGapWidth = fConfig.gapWidth;
	fAnodeWidth = fConfig.anodeWidth;
	fAnodePermittivity = fConfig.anodePermittivity;
	fCathodePermittivity = fConfig.cathodePermittivity;
	fCathodeWidth = fConfig.cathodeWidth;
	fDt = fDet->getTimeStep();
	fDx = fDet->getSpaceStep();
	
	fDiffL = fDet->getDiffL();
	fDiffT = fDet->getDiffT();
	fVx = vector<double> (iNstep,0);
	fAlpha = vector<double> (iNstep,0);
	fEta = vector<double> (iNstep,0);
	fE = vector<double> (iNstep,0);
	for(int i=0; i<iNstep; i++){
		fVx.at(i) = fDet->getVx();
		fAlpha.at(i) = fDet->getAlpha();
		fEta.at(i) = fDet->getEta();
		fE.at(i) = fDet->getElectricField()[0];
	}
	fEini = fE.at(0);
	fVini = fVx.at(0);
	
	fThrCLT = 1.5e4;
	fSpaceChargeLimit = 5.e7;
	
	iNElectronsSize = 5*iNstep;
	
	fChargeThres = fConfig.threshold * 1.e-12;//100.e-15; //pC to Coulombs
	iThrCrossTimeStep = -1;
	
	fElecDetectorGrid = vector<double> (iNstep,0);
	fPosIonDetectorGrid = vector<double> (iNstep,0);
	fNegIonDetectorGrid = vector<double> (iNstep,0);
	fNelecAnode = 0;
	
	fLongiDiffSigma = fDiffL*sqrt(fDx);

	bThrCrossTime = false;
	bHasReachSpaceChargeLimit = false;
	bEbarComputed = fDet->hasEBarTable();
	bComputeSpaceChargeEffet = true;
	bAvalancheInitialised = false;
	
	bVerbose = fConfig.verbose;
	bSnapshots = fConfig.snaps;
	iVerbosityLevel = fConfig.verbosityLevel;	
	
	bStreamer = false;
	fStreamerThr = 4.85e8;
	
	bDummyRun = false;
	bSimUntilThr = false;
	bOnlyMultiplicationAvalanche = fConfig.onlyMult;
	
	fDebugOutputs = fConfig.debugOutput;
	
	fLongiDiffTimeLimit = 1400;
	
	if ( bEbarComputed ) {
		iEbarTableSize = fDet->getEbarTableSize();
		int n = iEbarTableSize+1;
		fEbarTable = vector<double>(n*n*n,0);
		fEbarTable = fDet->getEbarVecTable();
	
		fEbarZarray = vector<double>(n);
		fEbarZparray = vector<double>(n);
		fEbarLarray = vector<double>(n);
	
		fEbarZarray = fDet->getEbarZarray();
		fEbarZparray = fDet->getEbarZparray();
		fEbarLarray = fDet->getEbarLarray();
	}
	
	fSignal.clear();
	fCharges.clear();
	
	writeRPCParameters();
	
	eAvalStatus = AVAL_NO_ERROR; //Avalanche status to NO_ERROR at begining
	
}

void TAvalanche1D::writeRPCParameters() {
	string outFileName = "out/analytic_parameters_"+toString(fE.at(0))+"_"+toString(fConfig.particleMomentum)+".dat";
	ofstream data(outFileName, ios::out | ios::trunc);
	
	data << "#alpha\teta\tg\tEw\tQthr" << endl;
	data << fAlpha.at(0) << "\t" <<  fEta.at(0) << "\t" << fGapWidth << "\t" << fAnodePermittivity/(fCathodeWidth+fAnodeWidth + fGapWidth*fAnodePermittivity) << "\t" << fChargeThres;
	
	data.close();
}

void TAvalanche1D::computeClusterDensity(const string& particleName, const double& Pmin, const double& Pmax, const int& steps){
	double P = Pmin;
	
	string outFileName = "out/cluster_density_"+particleName+".dat";
	ofstream data(outFileName.c_str(), ios::out | ios::trunc);
	
	Garfield::TrackHeed* track = new Garfield::TrackHeed();
	track->SetSensor(fDet->getSensor());
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

void TAvalanche1D::computeElectronsProduction(const string& particleName, const double& P, const int& nTracks){

	string outFileName = "out/electron_density_"+particleName+"_"+toString(P)+".dat";
	string clSizeFileName = "out/cluster_size_"+particleName+"_"+toString(P)+".dat";
	ofstream data(outFileName.c_str(), ios::out | ios::trunc);
	ofstream clSize(clSizeFileName.c_str(), ios::out | ios::trunc);
	
	Garfield::TrackHeed* track = new Garfield::TrackHeed();
	track->SetSensor(fDet->getSensor());
	track->SetParticle(particleName);
	track->SetMomentum(P);
	
	for(int i=0; i<nTracks; i++){
		track->NewTrack(0, 0, 0, 0, 1., 0, 0);
		double xc = 0., yc = 0., zc = 0., tc = 0.;
	    int nc = 0;
	    double ec = 0.;
	    double extra = 0.;
	    double esum = 0.;
	    int nsum = 0;
	    while (track->GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
	      esum += ec;
	      nsum += nc;
	      clSize << nc << endl;
	    }
		data << esum << "\t" << nsum << endl;
	}
	
	data.close();
	delete track;
}

void TAvalanche1D::initialiseTrackHeed(){
	if( bAvalancheInitialised ) {
		cerr << "TAvalanche1D::initialiseTrackHeed -- Avalanche already initialised." << endl;
		return;
	}
	
	if (fConfig.singleCluster) {
		initialiseSingleCluster();
		return;
	}
	
	if (iVerbosityLevel >= 1)
		cout << "Initialising avalanche (Id: " << Id <<") with Heed track (thread " << Id << ") for " << fConfig.particleName << " with momentum " << toString(fConfig.particleMomentum) << " at position " << toString(fConfig.x0) << " with angle " << toString(fConfig.theta) << endl;
	
	if (fConfig.garfieldSeed != -1)
		fSeed = fConfig.garfieldSeed;
	else
		fSeed = getUUID();
	
	fDet->setGarfieldSeed(fSeed);

	Garfield::TrackHeed* track = new Garfield::TrackHeed();
	//track->EnablePhotoAbsorptionCrossSectionOutput();
	track->SetSensor(fDet->getSensor());
	track->SetParticle(fConfig.particleName);
	track->SetMomentum(fConfig.particleMomentum);

	double t0 = 0.;
	double y0 = 0, z0 = 0;
	double dx0 = cos(fConfig.theta * Constants::Pi / 180.0), dy0 = sin(fConfig.theta * Constants::Pi / 180.0), dz0 = 0.;

	track->NewTrack(fConfig.x0, y0, z0, t0, dx0, dy0, dz0);

	double xc = 0., yc = 0., zc = 0., tc = 0.; 	// Cluster coordinates
	int nc = 0; 								// Number of electrons produced in a collision
	double ec = 0.; 							// Energy loss in a collision
	double dummy = 0.; 							// Dummy variable (not used at present)
	fNElectrons = vector<double> (iNElectronsSize,0);
	
	cout << "Cluster density for " << fConfig.particleName << " with momentum " << fConfig.particleMomentum << " : " << track->GetClusterDensity() << endl;
	
	while (track->GetCluster(xc, yc, zc, tc, nc, ec, dummy)){
		fClustersX[xc] = nc;
		fClustersY[yc] = nc;
		fClustersZ[zc] = nc;
		
		fNElectrons[0] += nc;
		fElecDetectorGrid[int(trunc(xc/fDx))] += nc;
		fPosIonDetectorGrid[int(trunc(xc/fDx))] += nc;
	}
	
	fClusterDensity = track->GetClusterDensity();
	delete track;
	
	if (iVerbosityLevel >= 2)
		cout << "TAvalanche1D::initialiseTrackHeed -- " << toString( fNElectrons[0] ) << " electrons produced by heed track." << endl;
	
	bAvalancheInitialised = true;
}

void TAvalanche1D::initialiseSingleCluster(){
	if( bAvalancheInitialised ) {
		cerr << "TAvalanche1D::initialiseSingleCluster -- Avalanche already initialised." << endl;
		return;
	}

	double step = fConfig.x0;
	int n0 = fConfig.n0;
	
	if (fConfig.x0 < 0)
		step = int(trunc(fRngMisc->RandU01() * iNstep));
	
	if (fConfig.n0 < 0)
		n0 = int( 80 * fRngMisc->RandU01() + 1 );
	
	if (iVerbosityLevel >= 1)
		cout << "Initialising avalanche (Id: " << Id <<") with " << n0 << " electron(s) at step " << step << endl;
	
	fNElectrons = vector<double> (iNElectronsSize,0);
	
	fNElectrons[0] = n0;
	fElecDetectorGrid[ int(trunc(step)) ] = n0;
	
	fClustersX[step] = n0;

	
	bAvalancheInitialised = true;
}

void TAvalanche1D::makeResultFile() {
	fResult.Dt = fDt;
	fResult.Dx = fDx;
	fResult.iNstep = iNstep;
	fResult.thrCrossTimeStep = iThrCrossTimeStep;
	fResult.avalStatus = eAvalStatus;
	fResult.computeTime = fElapsed;
	fResult.nCluster = fClustersX.size();
	/*fResult.charges_size = fCharges.size();
	fResult.chargesTot_size = fTotalCharges.size();
	fResult.signal_size = fSignal.size();*/
	fResult.size = fCharges.size();
	fResult.streamer = bStreamer;
	for (uint i=0; i<fCharges.size(); i++ )
		fResult.charges[i] = fCharges[i];
	for (uint i=0; i<fTotalCharges.size(); i++ )
		fResult.chargesTot[i] = fTotalCharges[i];
	for (uint i=0; i<fSignal.size(); i++ )
		fResult.signal[i] = fSignal[i];
	for (uint i=0; i<fNegIonDetectorGrid.size(); i++)
		fResult.nions[i] = fNegIonDetectorGrid[i];
	for (uint i=0; i<fPosIonDetectorGrid.size(); i++)
		fResult.pions[i] = fPosIonDetectorGrid[i];
	for (uint i=0; i<fCharges.size(); i++)
		fResult.nelec[i] = fNElectrons[i];
		
	int index=0;
	for (map<double, int>::iterator it=fClustersX.begin(); it!=fClustersX.end(); it++) {
		fResult.clNe[index] = it->second;
		fResult.clPos[index] = it->first;
		index++;
	}
		
}

void TAvalanche1D::simulateEvent(){
	/* check if avalanche has been initialised */
	if ( !bAvalancheInitialised ){
		printError(__FILE__, toString(__LINE__), __func__, "Avalanche has not been initialised. Aborting");
		exit(0);
	}
	
	/* if Space Charge Effect is enabled check if Ebar table has been computed */
	if ( bComputeSpaceChargeEffet && !bEbarComputed && !bOnlyMultiplicationAvalanche ){
		printError(__FILE__, toString(__LINE__), __func__, "No Ebar table found. Aborting simulation");
		exit(0); 
	}
	
	cout.clear();
	
	if( avalanche() ){
		const auto elapsed = fTimer.time_elapsed();
		fElapsed = static_cast<double>( duration_cast<seconds>(elapsed).count() );
		checkDetectorGrid();
		/* debug outputs */
		if (fDebugOutputs) {
			ofstream data("out/electrons.dat", ios::out | ios::trunc);
			ofstream sigData("out/signal.dat", ios::out | ios::trunc);
			ofstream chargesData("out/charges.dat", ios::out | ios::trunc);
			ofstream chargesTotData("out/chargesTot.dat", ios::out | ios::trunc);
			ofstream sigTotData("out/signalTot.dat", ios::out | ios::trunc);
			for(int i=0; i<iNElectronsSize; i++) data << fNElectrons[i] << endl;
			for(uint i=0; i<fSignal.size(); i++)	sigData << fSignal[i] << endl;
			for(uint i=0; i<fCharges.size(); i++)	chargesData << fCharges[i] << endl;
			for(uint i=0; i<fTotalCharges.size(); i++)	chargesTotData << fTotalCharges[i] << endl;
			for(uint i=0; i<fTotalSignal.size(); i++)	sigTotData << fTotalSignal[i] << endl;
			data.close();
			sigData.close();
			chargesData.close();
			chargesTotData.close();
			sigTotData.close();
		}
		/* ========= */
		makeResultFile();
		if (iVerbosityLevel >= 1)
			cout << currentDateTime() << " - Avalanche simulation id " << Id << " (" << countSim << "nth simulation) terminated with success (" << duration_cast<seconds>(elapsed).count() << " seconds)." << endl<< endl;
	}
	else{
		const auto elapsed = fTimer.time_elapsed();
		fElapsed = static_cast<double>( duration_cast<seconds>(elapsed).count() );
		fSignal.clear();
		fCharges.clear();
		fSignal.push_back(-1);
		fCharges.push_back(-1);
		makeResultFile();
		if (iVerbosityLevel >= 1)
			cout << currentDateTime() << " - Avalanche simulation id " << Id << " (" << countSim << "nth simulation) terminated with error: " << eAvalStatus << " (" << duration_cast<seconds>(elapsed).count() << " seconds)." << endl<< endl;
	}
	countSim++;
}

void TAvalanche1D::checkDetectorGrid(){
	for(int i=0; i<iNstep; i++) {
		DEBUGASSERT(fElecDetectorGrid[i] == 0);
		if (fElecDetectorGrid[i] != 0)
			eAvalStatus = AVAL_ERROR_GRID_NOT_EMPTY;
	}
}

void TAvalanche1D::computeInducedSignal2(){
	// drift velocity in cm/ns
	double e0 = Constants::ElectronCharge; //1.60217657e-19; //Coulombs
	double eps = fAnodePermittivity; 	//10.;
	double weightingField = eps/(fCathodeWidth+fAnodeWidth + fGapWidth*eps);

	double sig = 0;
	double charges = 0;

	for(int z=0; z < iNstep; z++){
		sig += weightingField * fVx.at(z)*1e9 * e0*fElecDetectorGrid[z];
		charges += weightingField * e0*fElecDetectorGrid[z] * fDx;
	}
	
	fSignal.push_back(sig);
	fCharges.push_back(charges);
	
	if (fTotalCharges.size() == 0)	{
		fTotalCharges.push_back(charges);
		fTotalSignal.push_back(sig);
	}
	else {
		fTotalCharges.push_back(fTotalCharges.back() + charges);
		fTotalSignal.push_back(fTotalSignal.back() + sig);
	}
	
	if (fTotalCharges.back() >= fChargeThres and !bThrCrossTime ){
		iThrCrossTimeStep = iTimeStep;
		bThrCrossTime = true;
	}
}

void TAvalanche1D::computeInducedCharges(){
	/* OBSOLETE */
	fCharges = vector<double>(1,fSignal[0]*fDt*1.e-9);
	ofstream data("out/induced_charges.dat", ios::out | ios::trunc);
	for(uint i=1; i<fSignal.size(); i++){
		fCharges.push_back(fCharges.back() + fSignal[i] * fDt*1.e-9);
		data << fCharges.back() << endl;
	}
	data.close();
}


inline double TAvalanche1D::n_moy(const double& x){
	return exp( ( fAlpha.at(iCurrentDetectorStep)-fEta.at(iCurrentDetectorStep) )*x );
}

double TAvalanche1D::multiplicationRiegler(const double& x, const double& s){
	double nm = n_moy(x);
	double alpha = fAlpha.at(iCurrentDetectorStep);
	double eta = fEta.at(iCurrentDetectorStep);
	double k = eta/alpha;
	double thr;
	
	
	if (alpha == eta){
		thr = alpha*x / (1.+alpha*x);
		if (s<=thr)	return 0;
		else{
			double val = 1 + trunc( log( (1.-s)*(1.+alpha*x) ) / log(thr)  );
			DEBUGASSERT(val>0);
				if (val < 0)
					eAvalStatus = AVAL_MULT_FAIL;
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
				/* Taylor expansion */
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
				DEBUGASSERT(val>0);
				if (val < 0)
					eAvalStatus = AVAL_MULT_FAIL;
				return val;
			}
		}	
	}
	
}

double TAvalanche1D::electronMultiplication(const double& n){
	double nProduced = 0;
	double s;
	
	if(n > fThrCLT){
		double c = multiplicationCLT(fDx,n);
		if( eAvalStatus  == AVAL_CLT_FAIL or eAvalStatus == AVAL_EXPLOSIVE_BEHAVIOR_CLT ) {
			/* CLT has failed to return an acceptable value. Try with classic multiplication */
			for(int i=0; i<n; i++){
				s = fRngMult->RandU01();
				if (s==1)	s = fRngMult->RandU01();
				nProduced += multiplicationRiegler(fDx,s);
			}
			if( eAvalStatus  == AVAL_CLT_FAIL or eAvalStatus == AVAL_EXPLOSIVE_BEHAVIOR_CLT )
				eAvalStatus = AVAL_NO_ERROR;
			return nProduced;
		}
		return c;
	}
	
	for(int i=0; i<n; i++){
		s = fRngMult->RandU01();
		if (s==1)	s = fRngMult->RandU01();
		nProduced += multiplicationRiegler(fDx,s);
	}
	return nProduced;
}

double TAvalanche1D::multiplicationCLT(const double& x, const double& n){
	double nm = n_moy(x);
	double k = fEta[iCurrentDetectorStep]/fAlpha[iCurrentDetectorStep];
	double alpha = fAlpha[iCurrentDetectorStep];
	double eta = fEta[iCurrentDetectorStep];
	double m = n*nm;
	
	double sigma;
	if (alpha == eta)	sigma = sqrt(n) * sqrt(2*alpha*x);
	else if (alpha < 0.01)	sigma = sqrt(n) * sqrt( exp(-2*eta*x)*(exp(eta*x)-1) ); // alpha << 0
	else	sigma = sqrt(n) * sqrt( ((1+k)/(1-k)) * nm * (nm-1) ); // alpha, eta > 0
	
	double c = Gaus(m, sigma, fRngCLT);
	
	if (c > 1e10) {
		eAvalStatus = AVAL_EXPLOSIVE_BEHAVIOR_CLT;
		return 0;
	}
	
	if (abs(c>1)) {
		if( !(trunc(c)>0) )
			eAvalStatus = AVAL_CLT_FAIL;
		return trunc(c);
	}
	
	else {
		if ( !(round(c) >= 0) ) 
			eAvalStatus = AVAL_CLT_FAIL;
		return round(c);
	}
}

void TAvalanche1D::computeLongitudinalDiffusion() {	
	vector<double> newDetectorGrid (iNstep,0);
	double pos, newPos, sigma;
	int newPosIndex;
	double sqrtDx = sqrt(fDx);
	
	fLongiDiffTimer.start();
	
	for(int iz=0; iz<iNstep; iz++){
		if ( checkTimerExceededLimit(fLongiDiffTimer,fLongiDiffTimeLimit) or eAvalStatus == AVAL_LONGI_DIFF_TIME_LIMIT_EXCEEDED ) {
			eAvalStatus = AVAL_LONGI_DIFF_TIME_LIMIT_EXCEEDED;
			break;
		}
		
		pos = (iz) * fDx;
		sigma = fDet->getDiffusionCoefficients(fE.at(iz), 0, 0)[0] * sqrtDx;
		
		for(int n=0; n<fElecDetectorGrid.at(iz); n++){
			/*if ( checkTimerExceededLimit(fLongiDiffTimer,fLongiDiffTimeLimit) ) {
				eAvalStatus = AVAL_LONGI_DIFF_TIME_LIMIT_EXCEEDED;
				break;
			}*/
			
			newPos = Gaus(pos, sigma, fRngLongiDiff);
			newPosIndex = (int)trunc(newPos/fDx);
			
			if ( newPosIndex >= iNstep)
				 newPosIndex = iNstep-1;
			else if (newPosIndex < 0)
				newPosIndex = 0;
			
			try{
				newDetectorGrid.at(newPosIndex)++;
			}
			catch (const std::out_of_range& oor) {
				printError(__FILE__, toString(__LINE__), __func__, oor.what());
				cerr << newPosIndex << " " << pos << " " << newPos << endl;
				exit(0);
			}
			
		}
	}

	fElecDetectorGrid = newDetectorGrid;
}

void TAvalanche1D::computeLongitudinalSCEffect() {
	vector<double> copy (fElecDetectorGrid);
	vector<double> newGrid (iNstep,0.);
	double n;
	double x,trx,Px;
	bool modified = false;
	
	for (int z=0; z<iNstep; z++) {
		n = copy.at(z);
		if (n == 0)
			continue;
			
		
		x = fVx.at(z)/fVini;
		trx = trunc(x);
		Px = x - trx;
		if (Px < 0.01 or Px>0.99) {
			newGrid.at(z) = n;
			continue;
		}
		
		cout << z << "\t" << x  << "\t" << n << "\t" << Px << "\t" << z+trx << ":" << n-trunc(Px*n) << "\t" << z+1+trx << ":" << trunc(Px*n) << endl;
		
		modified = true;
		
		if (z+trx >= iNstep-1) {
			newGrid.at(iNstep-1) += n;
			//continue;
		}		
		else if (z+1+trx >= iNstep-1) {
			newGrid.at(iNstep-1) += trunc(Px * n);
			newGrid.at(z+trx) += n - trunc(Px * n);
			//continue;
		}
		else {
			cout << 1 << endl;
		//if (z+1+trx < iNstep-1 and z+trx < iNstep-1) {
			newGrid.at(z+trx) += n - trunc(Px * n);
			newGrid.at(z+1+trx) += trunc(Px * n);
		}
	}
	
	if (modified) {
		cout << sumVec(newGrid) << " " << sumVec(fElecDetectorGrid) << endl;
		assert (sumVec(newGrid) == sumVec(fElecDetectorGrid));
		fElecDetectorGrid = newGrid;
	}
	
}

bool TAvalanche1D::checkForExplosiveBehavior() {
	/* 
	 * Check if a large number of elec has been created between two steps (typically meaning
	 * that the avalanche shows explosive behavior). If there is a production superior than 
	 * 30% between two steps, the avalanche simulation is aborted.
	 */
	 
	double n0 = fNElectrons.at(iTimeStep-1);
	double n1 = fNElectrons.at(iTimeStep);
	if ( 100.*(n1-n0)/n1 > 30. )
		return true;
	else
		return false;
}

bool TAvalanche1D::propagate() {
	double n, nProduced;
	vector<double> copy (fElecDetectorGrid);	
	
	for(iCurrentDetectorStep=0; iCurrentDetectorStep<iNstep-1; iCurrentDetectorStep++){
		n = copy.at(iCurrentDetectorStep);
		
		//if (bHasReachSpaceChargeLimit and !bComputeSpaceChargeEffet)
		//	nProduced = n;
		//else 
			nProduced = electronMultiplication(n);
		
		if (eAvalStatus != AVAL_NO_ERROR)
			return false;
		
		fElecDetectorGrid.at(iCurrentDetectorStep+1) = nProduced;
		
		if( (nProduced - n) > 0 )
			fPosIonDetectorGrid.at(iCurrentDetectorStep+1) += nProduced - n;
		else
			fNegIonDetectorGrid.at(iCurrentDetectorStep+1) += n - nProduced;
		
		fNElectrons[iTimeStep] += nProduced;
	}
	
	return true;
}

bool TAvalanche1D::avalanche() {
	iTimeStep = 1;
	
	while(true){
		if (bSnapshots)
			makeSnapshot();
		
		/* First we compute the space charge electric field and update parameters */
		if ( bComputeSpaceChargeEffet and !bOnlyMultiplicationAvalanche )
			computeSCEffect();
			
		if ( !propagate() )
			return false;
			
		//computeLongitudinalSCEffect();
		
		/* Empty the first bin after the first multplication procedure to avoid infinite elec creation */
		if (iTimeStep == 1)
			fElecDetectorGrid.at(0) = 0;
		
		/* Check for explosive avalanche */
		/*if (iTimeStep > 400) {
			if( checkForExplosiveBehavior() ){
				eAvalStatus = AVAL_POTENTIAL_EXPLOSIVE_BEHAVIOR;
				return false;
			}
		}*/
		
		if ( !bOnlyMultiplicationAvalanche )
			computeLongitudinalDiffusion();
		
		if (eAvalStatus != AVAL_NO_ERROR)
			return false;
		
		if (!bStreamer and fNElectrons[iTimeStep] >= fStreamerThr)
			bStreamer = true;
			
		computeInducedSignal2();
		
		if ( iTimeStep > iNElectronsSize ) {
			eAvalStatus = AVAL_ERROR_TIMESTEP_EXCEEDING_LIMIT;
			return false;
		}
		
		if (fNElectrons.at(iTimeStep) == 0)
			break;
		
		if (iVerbosityLevel >= 2) {
			cout << "time step: " << iTimeStep << "\t Nelec: " << fNElectrons[iTimeStep] << "\t" << "NelecLastBin: " << fNelecAnode;
			cout << " " << -sumVec(fPosIonDetectorGrid)+sumVec(fElecDetectorGrid)+sumVec(fNegIonDetectorGrid) << endl;
		}
		
		/*	Electrons in the last last bin have reached the anode during the last propagation, so we empty it */
		if (fElecDetectorGrid.at(iNstep-1) > 0){
			fNelecAnode += fElecDetectorGrid.at(iNstep-1);
			fElecOnAnode.push_back( make_pair(fElecDetectorGrid.at(iNstep-1),iTimeStep*fDx) );
			fElecDetectorGrid.at(iNstep-1) = 0;
		}
		
		if (bDummyRun and iTimeStep == 500)
			break;
			
		//if (!bComputeSpaceChargeEffet and !bHasReachSpaceChargeLimit and fNElectrons[iTimeStep]>fSpaceChargeLimit)
		//	bHasReachSpaceChargeLimit = true;
		
		if (bThrCrossTime and bSimUntilThr)
			break;
		
		iTimeStep++;
	}
	eAvalStatus = AVAL_SUCCESS;
	return true;
}

void TAvalanche1D::computeSCEffect() {
	double SCEField[iNstep];
	double tmp;
	vector<double> transportParams;
	
	for(int z=0; z<iNstep; z++){
		tmp = 0;
		for(int zp=0; zp<iNstep; zp++){
			if (fElecDetectorGrid[zp]==0 and fNegIonDetectorGrid[zp]==0 and fPosIonDetectorGrid[zp]==0)
				continue;
			tmp += (fElecDetectorGrid[zp] +fNegIonDetectorGrid[zp] -fPosIonDetectorGrid[zp]) *  interpolateEbar((z)*fDx*Constants::cm, (zp)*fDx*Constants::cm, iTimeStep*fDx); 
		}
		
		for(uint i=0; i<fElecOnAnode.size(); i++)
			tmp += (fElecOnAnode[i].first) * interpolateEbar((z)*fDx*Constants::cm, fGapWidth*Constants::cm, fElecOnAnode[i].second);
			
		SCEField[z] = tmp;
	}

	for(int z=0; z<iNstep; z++){
		fE.at(z) = fEini + SCEField[z];
		transportParams = fDet->getTransportParameters(fE[z],0.,0.);
		fAlpha.at(z) = transportParams.at(0);
		fEta.at(z) = transportParams.at(1);
		fVx.at(z) = -1. * transportParams.at(2);
	}
}

void TAvalanche1D::makeSnapshot() {
	if (iVerbosityLevel >= 2)
		cout << "Snapshot at: " << iTimeStep*fDt << endl;
	
	string fileName;
	if (iTimeStep<10)
		fileName = "out/snaps/snap-000"+toString(iTimeStep)+".dat";
	else if (iTimeStep<100 and iTimeStep>9) 
		fileName = "out/snaps/snap-00"+toString(iTimeStep)+".dat";
	else if (iTimeStep<1000 and iTimeStep>99) 
		fileName = "out/snaps/snap-0"+toString(iTimeStep)+".dat";
	else
		fileName = "out/snaps/snap-"+toString(iTimeStep)+".dat";
	
	ofstream data(fileName.c_str(), ios::out | ios::trunc);
	for(int i=0; i<iNstep; i++){
		data << i*fDx << "\t" << fElecDetectorGrid[i] << "\t" << fPosIonDetectorGrid[i] << "\t" << fNegIonDetectorGrid[i] << "\t" << fE[i] << "\t" << fAlpha[i] << "\t" << fEta[i] << endl;
	}
	data.close();
}

void TAvalanche1D::testInterpolation() {
	ofstream data("out/interp30.dat", ios::out | ios::trunc);
	
	//double z = 0.001;
	double z = 0.8*0.001;
	//double zp = 0.001;
	double zp = 0.8*0.001;
	//double l;
	
	vector<double> l = linspace(0,0.14,200);
	
	
/*	for(double& val : l) {
		//cout << i << endl;
		//l = i * 0.2/20;
		data << val << "\t" << interpolateEbar(z, zp, val) << endl;
	}
	*/
	data.close();
}

inline size_t TAvalanche1D::getIndex3D(const size_t& i, const size_t& j, const size_t& k) {
	return (long)i*(long)(iEbarTableSize+1)*(long)(iEbarTableSize+1) + (long)j*(long)(iEbarTableSize+1) + (long)k;
}

double TAvalanche1D::interpolateEbar(const double& z, const double& zp, const double& l) {  //x==z y==zp z==l
	bool debug = false;

	size_t iz0 = getLowerIndex(fEbarZarray.begin(),fEbarZarray.end(), z);
	size_t izp0 = getLowerIndex(fEbarZparray.begin(),fEbarZparray.end(), zp);
	size_t il0 = getLowerIndex(fEbarLarray.begin(),fEbarLarray.end(), l);
	
	size_t iz1 = iz0+1;
	size_t izp1 = izp0+1;
	size_t il1 = il0+1;
	
	double z0,z1,zp0,zp1,l0,l1;
	
	/* FIXME: the index can be greater than the size of fEbarLarray.
	 * To avoid crash we decrement the index if necessary.
	 * This is bad, ugly and bady-bad .... but temporary ....... I guess */
	
	if (il1 >= fEbarLarray.size())
		il1--;
	z0 = fEbarZarray.at(iz0); zp0 = fEbarZparray.at(izp0); l0 = fEbarLarray.at(il0);
	z1 = fEbarZarray.at(iz1); zp1 = fEbarZparray.at(izp1); l1 = fEbarLarray.at(il1);
	
	
	double zd = (z - z0)/(z1 - z0);
	double zpd = (zp - zp0)/(zp1 - zp0);
	double ld = (l - l0)/(l1 - l0);
	double c00, c10, c01, c11;
	size_t index;
	
	/* FIXME: the index3D can be greater than the size of EbarTable.
	 * To avoid crash we use the index variable and decrement if necessary
	 * This is bad, ugly and bady-bad .... but temporary ....... I guess */
	 
	//interpolate along z
	index = getIndex3D(iz1,izp1,il1);
	if ( index >= fEbarTable.size())
		index--;
	c00 = fEbarTable.at(getIndex3D(iz0,izp0,il0))*(1-zd) + fEbarTable.at(getIndex3D(iz1,izp0,il0))*zd;
	c10 = fEbarTable.at(getIndex3D(iz0,izp1,il0))*(1-zd) + fEbarTable.at(getIndex3D(iz1,izp1,il0))*zd;
	c01 = fEbarTable.at(getIndex3D(iz0,izp0,il1))*(1-zd) + fEbarTable.at(getIndex3D(iz1,izp0,il1))*zd;
	c11 = fEbarTable.at(getIndex3D(iz0,izp1,il1))*(1-zd) + fEbarTable.at(index)*zd;


	//interpolate along zp
	double c0 = c00*(1-zpd) + c10*zpd;
	double c1 = c01*(1-zpd) + c11*zpd;
	
	//interpolate along l
	double c = c0*(1-ld) + c1*ld;
	
	if(debug){
		cout << "==========================================================" << endl;
		cout << "z:" << z << "\t zp:" << zp << "\t l:" << l << endl;
		cout << "iz0:" << iz0 << "\t" << "izp0:" << izp0 << "\t" << "il0:" << il0 << endl;
		cout << "z0:" << z0 << "\t" << "zp0:" << zp0 << "\t" << "l0:" << l0 << endl;
		cout << "iz1:" << iz1 << "\t" << "izp1:" << izp1 << "\t" << "il1:" << il1 << endl;
		cout << "z1:" << z1 << "\t" << "izp1:" << zp1 << "\t" << "l1:" << l1 << endl;
		cout << "zd:" << zd << "\t zpd:" << zpd << "\t ld:" << ld << endl;
		cout << "c0:" << c0 << "\t c1:" << c1 << endl;
		cout << "test:" << fEbarTable[getIndex3D(iz0,izp0,il0)] << " " << fEbarTable[getIndex3D(iz1,izp1,il1)] << endl; 
		cout << "==========================================================" << endl;
	}
	
	return c;
}
