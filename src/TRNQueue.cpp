#include "TRNQueue.hpp"
#include "helper_functions.hpp"

using namespace std;

RNQueue::RNQueue(){
	fSize = 5e8;
	fRandRng = new RngStream();
	fQueue = queue<double> ();
	
	//generate();
}

RNQueue::RNQueue(unsigned long const& n){
	fSize = n;
	fRandRng = new RngStream();
	fQueue = queue<double> ();
	
	generate();
}

RNQueue::~RNQueue(){
	delete fRandRng;
}


void RNQueue::generate(){
	if ( !fQueue.empty() )
		return;
	
	cout << "Populating RN queue." << endl;
	
	for(unsigned long i=0; i<fSize; i++)
		//fQueue.push( fRandRng->RandU01() );
		fQueue.push( Gaus(0,1,fRandRng) );
}

double RNQueue::next(){
	if( fQueue.empty() )
		generate();
	
	double rn = fQueue.front();
    fQueue.pop();
    
    return rn;
}
