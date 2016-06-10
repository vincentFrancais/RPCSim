// TGenMT.hpp
// 
// Copyright 2016 Vincent Français <francais@clermont.in2p3.fr>
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301, USA.
// 
// 

#pragma once

#include <string>

const int N = 624;
const int M = 397;
const long MATRIX_A = 0x9908b0dfUL;   /* constant vector a */
const long UPPER_MASK = 0x80000000UL; /* most significant w-r bits */
const long LOWER_MASK = 0x7fffffffUL; /* least significant r bits */

class TGenMT {
	public:
		TGenMT();
		TGenMT(std::string filename);
		TGenMT(unsigned long s);
		TGenMT(unsigned long init_key[], int key_length);
		
		
		double mtRand();
		unsigned long genrand_int32();
		double genrand_real1();
		double genrand_real2();
		void restoreStatus(const char * inFileName);
		void saveStatus(char * inFileName);
		
	
	private:
		void init_genrand(unsigned long s);
		void init_by_array(unsigned long init_key[], int key_length);
		
		unsigned long mt[N]; /* the array for the state vector  */
		int mti; /* mti==N+1 means mt[N] is not initialized */
};

