#include "helper_functions.hpp"
#include "TMath.h"
#define PI 3.14159265

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

double* vecToArray(double* array, std::vector<double> vector){
	array = new double[vector.size()];
	for(uint i=0; i<vector.size(); i++)
		array[i] = vector[i];
	return array;
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

double generateGaussianNumber(double mu, double sigma, RngStream* stream){
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

double Gaus(double mean, double sigma, RngStream* stream){
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
   // Implementation from ROOT TRandom, modified for using RNGStream

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
