// TThreadData.hpp
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

#include "TDetector.hpp"
#include "TConfig.hpp"
#include "SFMT/SFMT.h"

/* Simple structure for passing simulation parameters to the wrapping function */
struct ThreadData{
	TDetector* detector;
	TConfig config;
	sfmt_t sfmt;
	int id;
	
	ThreadData (TDetector* det, TConfig conf, sfmt_t status, int i) : detector(det), config(conf), sfmt(status), id(i)	{ };
};
