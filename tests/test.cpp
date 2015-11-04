/*
 * test.cpp
 * 
 * Copyright 2015 francais <francais@clrtoport02>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "helper_functions.hpp"
#include "RngStream.h"
using namespace std;

void test1(){
	cout << "test1" << endl;
	RngStream* fRandRngLongiDiff = new RngStream();
	
	ofstream out("test1.dat", ios::out | ios::trunc);
	
	double n0 = 1000000;
	double* tab = new double[100];
	double* tabGauss = new double[100];
	double* tabTab = new double[100];
	for(int i=0; i<100; i++){
		tab[i] = 0;
		tabGauss[i] = 0;
		tabTab[i] = 0;
	}
	int iStep = 40;
	tab[iStep] = n0;
	
	double fLongiDiffSigma = 0.000170901;
	double fDx = 0.000285714;
	
	double* fLongiDiffFrac = new double[7];
	fLongiDiffFrac[0] = gaussIntegral(2000, -0.5*fDx, 0.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[1] = gaussIntegral(2000, 0.5*fDx, 1.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[2] = gaussIntegral(2000, -1.5*fDx, -0.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[3] = gaussIntegral(2000, 1.5*fDx, 2.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[4] = gaussIntegral(2000, -2.5*fDx, -1.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[5] = gaussIntegral(2000, 2.5*fDx, 3.5*fDx, 0, fLongiDiffSigma);
	fLongiDiffFrac[6] = gaussIntegral(2000, -3.5*fDx, -2.5*fDx, 0, fLongiDiffSigma);
	
	for(int i=0; i<n0; i++){
		double pos = (iStep+1)*fDx - 0.5*fDx;
		double newPos = Gaus(pos,fLongiDiffSigma,fRandRngLongiDiff);
		if(newPos > 0. and newPos < 0.2)	tabGauss[int(trunc(newPos/fDx))] += 1;
	}
	
	double n[7];
	double intPart;
	for(int i=0; i<7; i++){
		if ( modf(fLongiDiffFrac[i]*n0, &intPart) >= 0.5 ) n[i] = ceil(fLongiDiffFrac[i]*n0);
		else 	n[i] = floor(fLongiDiffFrac[i]*n0);
	}
	tabTab[iStep] += n[0];
	tabTab[iStep+1] += n[1];
	tabTab[iStep-1] += n[2];
	tabTab[iStep+2] += n[3];
	tabTab[iStep-2] += n[4];
	tabTab[iStep+3] += n[5];
	tabTab[iStep-3] += n[6];
	
	// The electrons "lost" by computing the fractions are put in the center bin
	tabTab[iStep] += n0 - sumArray(n,7);
	if(modf(sumArray(n,7), &intPart) != 0){
		for(int i=0; i<7; i++)	cout << n[i] << endl;
		cin.ignore();
	}
	
	for(int i=0; i<100; i++)	out << tabGauss[i] << " " << tabTab[i] << " " << abs(tabGauss[i]-tabTab[i]) << endl;
	out.close();
}

void test2(){
	cout << "test2" << endl;
	RngStream* fRandRngLongiDiff = new RngStream();
	
	ofstream out("test2.dat", ios::out | ios::trunc);
	
	double n0 = 1000000;
	double* tab = new double[100];
	double* tabGauss = new double[100];
	double* tabTab = new double[100];
	for(int i=0; i<100; i++){
		tab[i] = 0;
		tabGauss[i] = 0;
		tabTab[i] = 0;
	}
	int iStep = 50;
	tab[iStep] = n0;
	
	double fLongiDiffSigma = 0.000170901;
	double fDx = 0.000285714;
	
	double* fLongiDiffFrac = new double[100];
	int j = 48;
	for(int i=0; i<49; i++){
		fLongiDiffFrac[j] = gaussIntegral(2000, -(1.5+i)*fDx, -(0.5+i)*fDx, 0, fLongiDiffSigma);
		j--;
	}
	j = 50;
	for(int i=0; i<49; i++){
		fLongiDiffFrac[j] = gaussIntegral(2000, (0.5+i)*fDx, (1.5+i)*fDx, 0, fLongiDiffSigma);
		j++;
	}
	
	fLongiDiffFrac[49] = gaussIntegral(2000, -0.5*fDx, 0.5*fDx, 0, fLongiDiffSigma);
	
	for(int i=0; i<n0; i++){
		double pos = (iStep+1)*fDx - 0.5*fDx;
		double newPos = Gaus(pos,fLongiDiffSigma,fRandRngLongiDiff);
		if(newPos > 0. and newPos < 0.2)	tabGauss[int(trunc(newPos/fDx))] += 1;
	}
	
	//double n[100];
	double intPart;
	for(int i=0; i<100; i++){
		if ( modf(fLongiDiffFrac[i]*n0, &intPart) >= 0.5 ) tabTab[i] = ceil(fLongiDiffFrac[i]*n0);
		else 	tabTab[i] = floor(fLongiDiffFrac[i]*n0);
	}
	//tabTab[iStep] += n[0];
	//tabTab[iStep+1] += n[1];
	//tabTab[iStep-1] += n[2];
	//tabTab[iStep+2] += n[3];
	//tabTab[iStep-2] += n[4];
	//tabTab[iStep+3] += n[5];
	//tabTab[iStep-3] += n[6];
	
	// The electrons "lost" by computing the fractions are put in the center bin
	tabTab[iStep] += n0 - sumArray(tabTab,100);
	if(modf(sumArray(tabTab,100), &intPart) != 0){
		for(int i=0; i<100; i++)	cout << tabTab[i] << endl;
		cin.ignore();
	}
	
	for(int i=0; i<100; i++)	out << tabGauss[i] << " " << tabTab[i] << endl;
	out.close();
}

void test3(){
	RngStream* fRandRngLongiDiff = new RngStream();
	int n = 1e6;
	
	ofstream outGaus("test3Gaus.dat", ios::out | ios::trunc);
	
	for(int i=0; i<n; i++)	outGaus << Gaus(0,1,fRandRngLongiDiff) << endl;
	
	double dx = 0.1;
	double left[30] = {-5.,         -4.82758621, -4.65517241, -4.48275862, -4.31034483, -4.13793103,
 -3.96551724, -3.79310345, -3.62068966, -3.44827586, -3.27586207, -3.10344828,
 -2.93103448, -2.75862069, -2.5862069,  -2.4137931,  -2.24137931, -2.06896552,
 -1.89655172, -1.72413793, -1.55172414, -1.37931034, -1.20689655, -1.03448276,
 -0.86206897, -0.68965517, -0.51724138, -0.34482759, -0.17241379,  0.};
 
	double right[30] = {0.        ,  0.17241379,  0.34482759,  0.51724138,  0.68965517,
        0.86206897,  1.03448276,  1.20689655,  1.37931034,  1.55172414,
        1.72413793,  1.89655172,  2.06896552,  2.24137931,  2.4137931 ,
        2.5862069 ,  2.75862069,  2.93103448,  3.10344828,  3.27586207,
        3.44827586,  3.62068966,  3.79310345,  3.96551724,  4.13793103,
        4.31034483,  4.48275862,  4.65517241,  4.82758621,  5.        };
    
    ofstream outTab("test3Tab.dat", ios::out | ios::trunc);
    
	for(int i=0; i<29; i++) outTab << gaussIntegral(2000, left[i], left[i+1], 0, 1.) << endl;
	for(int i=0; i<29; i++) outTab << gaussIntegral(2000, right[i], right[i+1], 0, 1.) << endl;
	

}

int main(int argc, char **argv)
{
	test3();
	
	return 0;
}

