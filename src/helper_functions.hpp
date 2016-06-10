/*
 * helper_functions.hpp
 * 
 * Copyright 2016 Vincent Français <francais@clermont.in2p3.fr>
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


#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <utility>
#include <climits>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <sys/stat.h>
#include <numeric>
#include <ctime>
#include <functional>
#include <string>

//#include "RngStream.hpp"
#include "TRandomEngine.hpp"
#include "TRandomEngineMT.hpp"
#include "TRandomEngineMRG.hpp"
#include "TRandomEngineSFMT.hpp"
#include "TTimer.hpp"


#if defined( _PYTHON ) || defined( PYTHON ) || defined (__PYTHON__)
#   ifndef PYTHON
#       define PYTHON
#   endif
#endif

#if defined(PYTHON)
#	include <Python.h>
#endif

template<typename T>
std::string toString(const T & value)
{
	std::ostringstream oss;
	oss << value;
	return oss.str();
}

std::string currentDateTime();

uint64_t gettid();

size_t getUUID();

#if defined(PYTHON)
	int call_python_fun(std::string funName, std::vector<double> args, double& result);
#endif

std::string GetHexRepresentation(const unsigned char * Bytes, size_t Length);

bool file_exist (const std::string& name);

double* vecToArray(std::vector<double> vector);

std::vector<double> arrayToVec(double* array, size_t size);

double gauss(double x, double mean, double sigma);

double gaussIntegral(int steps, double min, double max, double mean, double sigma);

template<typename S>
S sumArray(S array[], int size){
	S sumArray = 0;
	for(int i=0; i<size; i++) 	sumArray += array[i];
	return sumArray;
}

double generateGaussianNumber(double mu, double sigma, TRandomEngine* stream);

//double Gaus(double mean, double sigma, RngStream* stream);
//double Gaus(double mean, double sigma, TRandomEngineSFMT* stream);
double Gaus(double mean, double sigma, TRandomEngine* stream);

template <typename Iterator, typename T>
std::size_t getLowerIndex(Iterator first, Iterator last, const T & val){
	Iterator lower = std::lower_bound(first,last,val);
	if (lower != first)	lower--;
	return std::distance(first, lower);
}

template<typename T>
T sumVec(std::vector<T> vec){
	return std::accumulate(vec.begin(), vec.end(), 0.0);
}

bool checkTimerExceededLimit(TTimer timer, double const& limit);

void testRNG(std::string const& rng);
