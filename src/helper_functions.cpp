#include "helper_functions.hpp"
#include "TMath.h"
#define PI 3.14159265

using namespace std;

/* Get current date/time, format is YYYY-MM-DD.HH:mm:ss
*  From http://stackoverflow.com/a/10467633
*/
string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

/* Return thread id (use with caution, those are not PID). */
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

/* return the hexadeximal representation of an array of char. */
std::string GetHexRepresentation(const unsigned char * Bytes, size_t Length){
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
}

/* Doesn't work */
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

/* Return an Universal Unique Id, hashed from random bits from /dev/random. */
uint getUUID() {
	uint t = 0;
	
	ifstream f;
	f.open("/dev/random", std::ios::binary | std::ios::in);
	
	if (!f.good()) {
		cerr << "Error reading /dev/random" << endl;
		exit(0);
	}
	
	f.read(reinterpret_cast<char *>(&t), sizeof(t));
	f.close();
	
	hash<uint> uinthash;
	return uinthash(t);
}

/* Check if timer as exceeded the limit. Used to control that some functions doesn't
 * too long to compute, otherwise it could mean we hit a singularity. */
bool checkTimerExceededLimit(TTimer timer, double const& limit) {
	auto elapsed = timer.time_elapsed();
	if ( static_cast<double>( duration_cast<seconds>(elapsed).count() ) >= limit)
		return true;
	else
		return false;
}


/* Create file of outputs for the different RNGs for repeatability checks. */
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

/* Equality checks for floating point numbers, classic == could give errors due to
 * rounding-up errors. */
bool almostEquals(double a, double b, double epsilon) {
    return std::abs(a - b) < epsilon;
}

double bessel_J0 (double X) {
/* **********************************************************************
      This subroutine calculates the First Kind Bessel Function of
      order 0, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
      
      From http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/tbessj_cpp.txt
*********************************************************************** */
	const double
		P1=1.0, P2=-0.1098628627E-2, P3=0.2734510407E-4,
        P4=-0.2073370639E-5, P5= 0.2093887211E-6,
		Q1=-0.1562499995E-1, Q2= 0.1430488765E-3, Q3=-0.6911147651E-5,
		Q4= 0.7621095161E-6, Q5=-0.9349451520E-7,
		R1= 57568490574.0, R2=-13362590354.0, R3=651619640.7,
		R4=-11214424.18, R5= 77392.33017, R6=-184.9052456,
		S1= 57568490411.0, S2=1029532985.0, S3=9494680.718,
		S4= 59272.64853, S5=267.8532712, S6=1.0;
		
	
	double AX,FR,FS,Z,FP,FQ,XX,Y, TMP;

      if (X==0.0) 
		return 1.0;
      
	AX = fabs(X);
	if (AX < 8.0) {
		Y = X*X;
        FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
        FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
        TMP = FR/FS;
	}
	else {
		Z = 8./AX;
        Y = Z*Z;
        XX = AX-0.785398164;
        FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
        FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
        TMP = sqrt(0.636619772/AX)*(FP*cos(XX)-Z*FQ*sin(XX));
	}
	
	return TMP;
}

/* Return a lineary spaced vector from start to end with num items.
 * Similar to numpy's linspace. 
 * From http://stackoverflow.com/questions/27028226/python-linspace-in-c */
std::vector<double> linspace(double start, double end, int num){
	double delta = (end - start) / (num - 1);
	
	std::vector<double> linspaced(num-1);
	for(int i=0; i < num-1; ++i)
		linspaced.at(i) = start + delta * i;
	
	linspaced.push_back(end);
	return linspaced;
}

void printError(std::string file, std::string line, std::string func, std::string what) {
	std::cerr << "Error -- file " << file << ", function " << func << " (l." << line << ") -- ";
	std::cerr << what << endl; 
}
