// TRandomEngineMRG.hpp
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

#include "TRandomEngine.hpp"
#include "RngStream.hpp"


class TRandomEngineMRG : public TRandomEngine {
	public:
		TRandomEngineMRG() : TRandomEngine() { rng = RngStream(); }
		TRandomEngineMRG(std::string s) : TRandomEngine() { rng = RngStream(s.c_str()); } 
		
		double RandU01() { return rng.RandU01(); }
		std::string Generator() { return "MRG32k3a"; }
		
	private:
	RngStream rng;
};



