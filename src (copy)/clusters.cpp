#include "clusters.hpp"

using namespace std;

void Clusters::addCluster(int pos, double nElec){
	fPos.push_back(pos);
	fNElec.push_back(nElec);
}
