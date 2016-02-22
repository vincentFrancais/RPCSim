#pragma once

#include <queue>

#include "RngStream.h"

using namespace std;

class RNQueue {
	public:
	RNQueue();
	RNQueue(unsigned long const& n);
	
	~RNQueue();
	
	double next();
	
	private:
	void generate();
	
	unsigned long fSize;
	RngStream* fRandRng;// = new RngStream();
	queue<double> fQueue;
};
