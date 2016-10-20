#include "helper_functions.hpp"
#include "TMath.h"
#define PI 3.14159265

using namespace std;

#if defined(PYTHON)
	int call_python_fun(std::string funName, std::vector<double> args, double& result){
	    //double result = 0.;
	    Py_Initialize();
	
		// Set the path to include the current directory in case the module is located there.
		PyObject *sys = PyImport_ImportModule("sys");
		PyObject *path = PyObject_GetAttrString(sys, "path");
		PyList_Append(path, PyString_FromString("./python/"));
	
		PyObject *pName, *pModule, *pFunc, *pResult = NULL;
		PyObject *pArgs, *pValue = NULL;
	
	    // Build the name object
		pName = PyString_FromString((char*)"integration");
	
		// Load the module object
	  	pModule = PyImport_Import(pName);
	    Py_DECREF(pName);
	
	    if(pModule != NULL){
			// pDict is a borrowed reference
	  		//pDict = PyModule_GetDict(pModule);
	
			// pFunc is also a borrowed reference
			//pFunc = PyDict_GetItemString(pDict, (char*)"compute_pot_correction_term");
	        pFunc = PyObject_GetAttrString(pModule, funName.c_str());
	
			if (PyCallable_Check(pFunc)){
				pArgs = PyTuple_New( args.size() );
				for(uint i=0; i<args.size(); i++){
					// r	phi	z	rp	phip	zp	eps1	eps2	eps3	p	q	g
					pValue = PyFloat_FromDouble(args[i]);
					if (!pValue) {
						Py_DECREF(pArgs);
						Py_DECREF(pModule);
	                    Py_DECREF(pFunc);
	    				Py_DECREF(pModule);
	                    //Py_DECREF(pDict);
						cerr << "Cannot convert value" << endl;
	                    //Py_Finalize();
						return 1;
					}
					PyTuple_SetItem(pArgs, i, pValue);
				}
	
	            // Call to the python function
				pResult = PyObject_CallObject(pFunc, pArgs);
	            Py_DECREF(pArgs);
	
	            if (pResult != NULL) {
	                result = PyFloat_AsDouble( pResult );
					Py_DECREF(pResult);
				}
	            else {
	                Py_DECREF(pFunc);
					Py_DECREF(pModule);
	                //Py_DECREF(pDict);
					PyErr_Print();
					cerr << "Call to pyhon function failed" << endl;
	                //Py_Finalize();
	                return 1;
	            }
	
	  		}
			else {
				if ( PyErr_Occurred() )
					PyErr_Print();
				cerr << "Cannot load function" << endl;
			}
	
	        Py_XDECREF(pFunc);
	        //Py_XDECREF(pDict);
			Py_DECREF(pModule);
		}
		else{
			cerr << "error loading python module" << endl;
			PyErr_Print();
		}
	
		//Py_Finalize();
	
		return 0;
	}
#endif


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
// From http://stackoverflow.com/a/10467633
string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

uint64_t gettid() {
    pthread_t ptid = pthread_self();
    uint64_t threadId = 0;
    memcpy(&threadId, &ptid, std::min(sizeof(threadId), sizeof(ptid)));
    return threadId;
}


bool file_exist (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

std::string GetHexRepresentation(const unsigned char * Bytes, size_t Length)
{
    std::ostringstream os;
    os.fill('0');
    os<<std::hex;
    for(const unsigned char * ptr=Bytes;ptr<Bytes+Length;ptr++)
        os<<std::setw(2)<<(unsigned int)*ptr;
    return os.str();
}

double* vecToArray(std::vector<double> vector){
	double* array = new double[vector.size()];
	for(uint i=0; i<vector.size(); i++)
		array[i] = vector[i];
	return array;
}

std::vector<double> arrayToVec(double* array, size_t size){
	std::vector<double> vec(size);
	for(uint i=0; i<size; i++)
		vec[i] = array[i];
	return vec;
}

double gauss(double x, double mean, double sigma){
	return 1/(sigma*sqrt(2*PI)) * exp( -pow((x-mean),2)/(2*sigma*sigma)  );
	//1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2) )
}

double gaussIntegral(int steps, double min, double max, double mean, double sigma){
	double val[steps];
	val[0] = min;
	val[steps-1] = max;
	for(int i=1; i<steps-1; i++)	val[i] = val[i-1] + (max-min)/steps;

	double integral = 0;

	for(int i=0; i<steps-1; i++){
		double a = val[i];
		double b = val[i+1];
		integral += gauss(a,mean,sigma) * (b - a);
	}

	return integral;
}

double generateGaussianNumber(double mu, double sigma, TRandomEngine* stream){
	// Simple Box-Muller gaussiam random number generator

	const double epsilon = std::numeric_limits<double>::min();
	const double two_pi = 2.0*3.14159265358979323846;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)	return z1 * sigma + mu;

	double u1, u2;
	do{
	   u1 = stream->RandU01();
	   u2 = stream->RandU01();
	}
	while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

double Gaus(double mean, double sigma, TRandomEngine* stream){
   // Samples a random number from the standard Normal (Gaussian) Distribution
   // with the given mean and sigma.
   // Uses the Acceptance-complement ratio from W. Hoermann and G. Derflinger
   // This is one of the fastest existing method for generating normal random variables.
   // It is a factor 2/3 faster than the polar (Box-Muller) method used in the previous
   // version of TRandom::Gaus. The speed is comparable to the Ziggurat method (from Marsaglia)
   // implemented for example in GSL and available in the MathMore library.
   //
   // REFERENCE:  - W. Hoermann and G. Derflinger (1990):
   //              The ACR Method for generating normal random variables,
   //              OR Spektrum 12 (1990), 181-185.
   //
   // Implementation taken from
   // UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien
   //
   // Implementation from ROOT TRandom, modified for using different RNG

   const double kC1 = 1.448242853;
   const double kC2 = 3.307147487;
   const double kC3 = 1.46754004;
   const double kD1 = 1.036467755;
   const double kD2 = 5.295844968;
   const double kD3 = 3.631288474;
   const double kHm = 0.483941449;
   const double kZm = 0.107981933;
   const double kHp = 4.132731354;
   const double kZp = 18.52161694;
   const double kPhln = 0.4515827053;
   const double kHm1 = 0.516058551;
   const double kHp1 = 3.132731354;
   const double kHzm = 0.375959516;
   const double kHzmp = 0.591923442;
   /*zhm 0.967882898*/

   const double kAs = 0.8853395638;
   const double kBs = 0.2452635696;
   const double kCs = 0.2770276848;
   const double kB  = 0.5029324303;
   const double kX0 = 0.4571828819;
   const double kYm = 0.187308492 ;
   const double kS  = 0.7270572718 ;
   const double kT  = 0.03895759111;

   double result;
   double rn,x,y,z;

   do {
      y = stream->RandU01();

      if (y>kHm1) {
         result = kHp*y-kHp1; break; }

      else if (y<kZm) {
         rn = kZp*y-1;
         result = (rn>0) ? (1+rn) : (-1+rn);
         break;
      }

      else if (y<kHm) {
         rn = stream->RandU01();
         rn = rn-1+rn;
         z = (rn>0) ? 2-rn : -2-rn;
         if ((kC1-y)*(kC3+TMath::Abs(z))<kC2) {
            result = z; break; }
         else {
            x = rn*rn;
            if ((y+kD1)*(kD3+x)<kD2) {
               result = rn; break; }
            else if (kHzmp-y<exp(-(z*z+kPhln)/2)) {
               result = z; break; }
            else if (y+kHzm<exp(-(x+kPhln)/2)) {
               result = rn; break; }
         }
      }

      while (1) {
         x = stream->RandU01();
         y = kYm * stream->RandU01();
         z = kX0 - kS*x - y;
         if (z>0)
            rn = 2+y/x;
         else {
            x = 1-x;
            y = kYm-y;
            rn = -(2+y/x);
         }
         if ((y-kAs+x)*(kCs+x)+kBs<0) {
            result = rn; break; }
         else if (y<x+kT)
            if (rn*rn<4*(kB-log(x))) {
               result = rn; break; }
      }
   } while(0);

   return mean + sigma * result;
}

uint getUUID() {
	uint t = 0;
	
	ifstream f;
	f.open("/dev/urandom", std::ios::binary | std::ios::in);
	
	if (!f.good()) {
		cerr << "Error reading /dev/urandom" << endl;
		exit(0);
	}
	
	f.read(reinterpret_cast<char *>(&t), sizeof(t));
	f.close();
	
	hash<uint> uinthash;
	return uinthash(t);
}

bool checkTimerExceededLimit(TTimer timer, double const& limit) {
	auto elapsed = timer.time_elapsed();
	if ( static_cast<double>( duration_cast<seconds>(elapsed).count() ) >= limit)
		return true;
	else
		return false;
}

void testRNG(string const& rng) {
	if (rng == "MT" or rng == "Mersenne-Twister" or rng == "mt") {
		ofstream mt_out("out/testRNG/mt.out", ios::out | ios::trunc);
		unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
		TRandomEngineMT mt = TRandomEngineMT(init,length);
		for (int i=0; i<1000; i++) {
			mt_out << mt.RandU01() << " ";
			if (i%5==4) mt_out << endl;
		}
		mt_out.close();
	}
	else if (rng == "MRG" or rng == "mrg" or rng == "RngStream") {
		ofstream mrg_out("out/testRNG/mrg.out", ios::out | ios::trunc);
		TRandomEngineMRG mrg = TRandomEngineMRG("g1");
		for (int i=0; i<1000; i++) {
			mrg_out << mrg.RandU01() << " ";
			if (i%5==4) mrg_out << endl;
		}
		mrg_out.close();
	}
	else if (rng == "sfmt" or rng == "SFMT") {
		ofstream sfmt_out("out/testRNG/sfmt.out", ios::out | ios::trunc);
		TRandomEngineSFMT sfmt = TRandomEngineSFMT(1234);
		for (int i=0; i<1000; i++) {
			sfmt_out << sfmt.RandU01() << " ";
			if (i%5==4) sfmt_out << endl;
		}
		sfmt_out.close();
	}
}

bool almostEquals(double a, double b, double epsilon) {
    return std::abs(a - b) < epsilon;
}
