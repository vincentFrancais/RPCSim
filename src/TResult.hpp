#ifndef DEF_RESULT
#define DEF_RESULT

using namespace std;

struct TResult
{
	//double* fInducedCharge;
	//vector<double> fInducedSignal;
	//double* fClPosX,fClPosY,fClPosZ;//vector<double> fClPosX,fClPosY,fClPosZ;
	//double* fNElectrons;
	//double fVx;
	//double fDiffCoeff[2];
	//double alpha;
	//double eta;
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
	//double fGapWidth;
	//string sMixture;
	//double fThrCLT;
	//double fSpaceChargeLimitThr;
    //unsigned char fId[SHA384_DIGEST_LENGTH];
    //double fPosition[kMax];
    //double fImpulsion[kMax];
    //double fEnergy;
    //unsigned int fPlans;
};

#endif
