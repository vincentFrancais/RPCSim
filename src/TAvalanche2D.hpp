#pragma once


using namespace std;

class TAvalanche2D {
	public:
	double n_moy(const double& x);
	double electron_multiplication(const double& x, const double& s);
	double multiplication(const double& n);
	double CLT(const double& x, const double& n);
	
	private:
	TDetector fDetector;
	DetectorGeometry fGeometry;

};
