// TResult.hpp
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


#ifndef DEF_RESULT
#define DEF_RESULT

using namespace std;

struct TResult
{
	double Dx;
	double Dt;
	int iNstep;
	int thrCrossTimeStep;
	int avalStatus;
	double computeTime;
	int streamer;
	int nCluster;
	uint size;
	double charges[2000];
	double chargesTot[2000];
	double signal[2000];
	double pions[2000];
	double nions[2000];
	double nelec[2000];
	double clPos[2000];
	double clNe[2000];
};

#endif
