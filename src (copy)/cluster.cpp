#include "cluster.hpp"

using namespace std;

Cluster::Cluster(double x, double y, double z, double nElec) : _x(x), _y(y), _z(z), _nElec(nElec) {}
Cluster::Cluster(double pos[], double nElec) : _x(pos[0]), _y(pos[1]), _z(pos[2]), _nElec(nElec) {}

vector<double> Cluster::getPos(){
	double pos[] = {_x,_y,_z}; 
	return std::vector<double> (pos, pos + sizeof(pos) / sizeof(double) );
}

Electron::Electron(double x, double y, double z) {
	_x.push_back(x);
	_y.push_back(y);
	_z.push_back(z);
}

vector<double> Electron::getPos() {
	double pos[] = {_x.back(),_y.back(),_z.back()}; 
	return vector<double> (pos, pos + sizeof(pos) / sizeof(double) );
}

void Electron::propagate(double Dx) {
	_x.push_back( _x[_x.size()-1] + Dx );
}
