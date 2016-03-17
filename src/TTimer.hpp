#pragma once

#include <chrono>
#include <thread>
#include <iostream>

using namespace std::chrono;

class TTimer {
	public:
		typedef high_resolution_clock Clock;
		void start()	{ epoch = Clock::now(); }
		Clock::duration time_elapsed() const	{ return Clock::now() - epoch; }
		
	private:
		Clock::time_point epoch;
};
