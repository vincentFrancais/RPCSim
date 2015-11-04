#include "helper_functions.hpp"
#define PI 3.14159265

bool file_exist (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
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

double generateGaussianNumber(double mu, double sigma, RngStream stream){
	const double epsilon = std::numeric_limits<double>::min();
	const double two_pi = 2.0*3.14159265358979323846;
 
	static double z0, z1;
	static bool generate;
	generate = !generate;
 
	if (!generate)	return z1 * sigma + mu;
 
	double u1, u2;
	do{
	   //u1 = rand() * (1.0 / RAND_MAX);
	   //u2 = rand() * (1.0 / RAND_MAX);
	   u1 = stream.RandU01();
	   u2 = stream.RandU01();
	}
	while ( u1 <= epsilon );
 
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}
