#ifndef DEF_CLUSTER
#define DEF_CLUSTER

#include <vector>

using namespace std;

class Cluster{
    public:
    
	Cluster(double x, double y, double z, double nElec);
	Cluster(double pos[], double nElec);
	
	vector<double> getPos();
	void setPosX(double x) {_x = x;}
	double getNElec(){return _nElec;}
	void setNElec(double n) {_nElec = n;}

    private:

    double _x;
    double _y;
    double _z;
    //double _pos[3] = {_x,_y,_z};
    double _nElec;
};

class Electron{
	public:
	
	Electron(double x, double y, double z);
	
	vector<double> getPos();
	void propagate(double Dx);
	
	
	private:
	
	vector<double> _x;
	vector<double> _y;
	vector<double> _z;
};

#endif
