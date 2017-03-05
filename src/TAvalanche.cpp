#include "TAvalanche.hpp"

using namespace std;

TAvalanche::TAvalanche () {

	fTimer.start();
}


TAvalanche::~TAvalanche() {
	
}

void TAvalanche::computeClusterDensity(const TDetector* det, const string& particleName, const double& Pmin, const double& Pmax, const int& steps){
	cout << "TAvalanche::Computing cluster density." << endl;
	double P = Pmin;
	string outFileName = "out/cluster_density_"+particleName+"_"+det->getGasName()+"_"+toString(det->getConfig().gapWidth)+".dat";
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

void TAvalanche::computeElectronsProduction(const TDetector* det, const string& particleName, const double& P, const int& nTracks){
	cout << "TAvalanche::Computing electron density and cluster size." << endl;

	string outFileName = "out/electron_density_"+particleName+"_"+toString(P)+"_"+det->getGasName()+".dat";
	string clSizeFileName = "out/cluster_size_"+particleName+"_"+toString(P)+"_"+det->getGasName()+".dat";
	string clDistribFileName = "out/cluster_distribution_"+particleName+"_"+toString(P)+"_"+det->getGasName()+".dat";
	ofstream data(outFileName.c_str(), ios::out | ios::trunc);
	ofstream clSize(clSizeFileName.c_str(), ios::out | ios::trunc);
	ofstream clDistrib(clDistribFileName.c_str(), ios::out | ios::trunc);
	
	Garfield::TrackHeed* track = new Garfield::TrackHeed();
	track->SetSensor(det->getSensor());
	track->SetParticle(particleName);
	track->SetMomentum(P);
	
	for(int i=0; i<nTracks; i++){
		track->NewTrack(0, 0, 0, 0, 1., 0., 0);
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
	      clDistrib << xc << endl;
	    }
		data << esum << "\t" << nsum << endl;
	}
	
	data.close();
	clSize.close();
	clDistrib.close();
	delete track;
}
