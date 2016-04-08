// TRandomEngine.hpp
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

class TRandomEngine {

	public:
	// Constructor
	TRandomEngine() {}
	// Destructor
	virtual ~TRandomEngine() {}

	// Draw a random number
	virtual double RandU01() = 0;
	// Initialise the random number generator
	//virtual void Seed(unsigned int s) = 0;
	// Give the generator type
	virtual std::string Generator() = 0;
};

