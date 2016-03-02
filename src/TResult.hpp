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
	uint charges_size;
	uint chargesTot_size;
	uint signal_size;
	double charges[1000];
	double chargesTot[1000];
	double signal[1000];
};

#endif
