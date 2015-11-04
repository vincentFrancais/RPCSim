#ifndef DEF_CLUSTERS
#define DEF_CLUSTERS

#include <vector>

using namespace std;

class Clusters{
    public:
    
	Clusters(){};
	void addCluster(int pos, double nElec);
	

    private:

	vector<int> fPos;
	vector<double> fNElec;
};

#endif
