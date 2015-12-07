#ifndef DEF_HELPERS
#define DEF_HELPERS

#include <vector>
#include <sstream>
#include <iomanip>
#include <utility>
#include <limits>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <sys/stat.h>
#include <numeric>

#include "RngStream.h"


template<typename T>
std::string to_string(const T & value)
{
	std::ostringstream oss;
	oss << value;
	return oss.str();
}

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

double generateGaussianNumber(double mu, double sigma, RngStream* stream);

double Gaus(double mean, double sigma, RngStream* stream);

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

#endif
