#ifndef DEF_HELPERS
#define DEF_HELPERS

#include <sstream>
#include <utility>
#include <limits>
#include <cmath>
#include <sys/stat.h>

#include "RngStream.h"


template<typename T>
std::string to_string(const T & value)
{
	std::ostringstream oss;
	oss << value;
	return oss.str();
}

bool file_exist (const std::string& name);

double gauss(double x, double mean, double sigma);

double gaussIntegral(int steps, double min, double max, double mean, double sigma);

template<typename S>
S sumArray(S array[], int size){
	S sumArray = 0;
	for(int i=0; i<size; i++) 	sumArray += array[i];
	return sumArray;
}

double generateGaussianNumber(double mu, double sigma, RngStream stream);

double Gaus(double mean, double sigma, RngStream* stream);

#endif
