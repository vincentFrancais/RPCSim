// TRandomEngineSFMT.hpp
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

#define SFMT_MEXP 19937

#include "TRandomEngine.hpp"
#include "SFMT/SFMT.h"



class TRandomEngineSFMT : public TRandomEngine {
	public:
		TRandomEngineSFMT(sfmt_t status) : TRandomEngine(), fSFMT(status) {};
		TRandomEngineSFMT() : TRandomEngine() { sfmt_init_gen_rand(&fSFMT, 4321); }
		TRandomEngineSFMT(int seed) : TRandomEngine() { sfmt_init_gen_rand(&fSFMT, seed); }
		
		double RandU01() { return sfmt_genrand_real2( &fSFMT ); }
		std::string Generator() { return "SIMD-oriented Fast Mersenne Twister"; }
		
	private:
	sfmt_t fSFMT;
};


