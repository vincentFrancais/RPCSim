
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <limits>
#include <iterator>
#include <algorithm>

#include "TAvalanche2D.hpp"

using namespace std;


TAvalanche2D::TAvalanche2D(TDetector* det, TConfig& config, sfmt_t sfmt, const int& id) : TAvalanche() {
	fDetector = det;
	fConfig = config;
}
