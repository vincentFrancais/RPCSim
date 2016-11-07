
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <assert.h>
#include <fcntl.h>
#include <functional>
#include <string>

#include "TDetector.hpp"
#include "helper_functions.hpp"
#include "integration.hpp"
#include "TConstants.hpp"


//#include "gsl/gsl_sf_bessel.h"
#include <gsl/gsl_integration.h>
//#include "gsl/gsl_math.h"


using namespace std;
static TDetector* tgsl = 0;


TDetector::TDetector(const TConfig& config){
	fConfig = config;
	iNstep = fConfig.nSteps;
	//fGeometry = geometry;
	if (fConfig.EbarTableCalculationSteps != 0)
		iEbarTableSize = fConfig.EbarTableCalculationSteps;
	else
		iEbarTableSize = 15;

	bHasEbarTable = false;
	bGasLoaded = false;
	bDetectorInitialised = false;
	
	setGasMixture();
	setElectricField(fConfig.ElectricField,0.,0.);
	//initialiseDetector();
	
	// if we do not simulate avalanche, no point in computing Ebar table.
	
}

TDetector::~TDetector(){
	delete mGas;
	delete mSensor;
}


void TDetector::printPACSData(Garfield::MediumMagboltz* gas) {
	/* Print PhotoAbsorption Cross Section data to txt file heed_pacs.txt (the file name is 
	 * hardcoded in heed and cannot be changed). 
	 * Takes a pointer to a Magboltz gas as parameter.
	 * The printing is done through the function TrackHeed::SetupGas which is private. So we
	 * define a dummy TrackHeed (with simple geometry and particle) which itself will call SetupGas. */
	 
	cout << "Printing PACS data to file heed_pacs.txt" << endl;
	cout.setstate(ios_base::failbit);	// To avoid unecessary Garfield printing. 
	
	Garfield::TrackHeed* track = new Garfield::TrackHeed;
	track->EnablePhotoAbsorptionCrossSectionOutput();
	
	// Dummy geometry, component and sensor
	Garfield::SolidBox* box = new Garfield::SolidBox(0,0,0,10,10,10);
	Garfield::GeometrySimple* geo = new Garfield::GeometrySimple();
	geo->AddSolid(box, gas);

	Garfield::ComponentConstant* comp = new Garfield::ComponentConstant();
	comp->SetGeometry(geo);
	comp->SetElectricField(10,10,10);

	Garfield::Sensor* sensor = new Garfield::Sensor();
	sensor->AddComponent(comp);
	
	//Dummy track
	track->SetSensor(sensor);
	track->SetParticle("muon");
	track->SetMomentum(5.e9);
	track->NewTrack(0, 0, 0, 0, 1, 0, 0);	// <=== printing is done here!
	
	cout.clear();
	
	delete track;
	delete sensor;
	delete comp;
	delete geo;
	delete box;
}

void TDetector::printPACSData() {
	/* Print PhotoAbsorption Cross Section data to txt file heed_pacs.txt (the file name is 
	 * hardcoded in heed and cannot be changed). 
	 * Use the gas defined in the config file (initialised with TDetecor) 
	 */
	 
	cout << "Printing PACS data to file heed_pacs.txt" << endl;
	cout.setstate(ios_base::failbit);	// To avoid unecessary Garfield printing. 
	
	Garfield::TrackHeed* track = new Garfield::TrackHeed;
	track->EnablePhotoAbsorptionCrossSectionOutput(); // <=== enable the printing of PACS data
	
	//Dummy track
	track->SetSensor(mSensor);
	track->SetParticle("muon");
	track->SetMomentum(5.e9);
	track->NewTrack(0, 0, 0, 0, 1, 0, 0);	// <=== printing is done here!
	
	cout.clear();
	
	delete track;
}

string TDetector::getGasName() const {
	string gasName = "";
	int nComponents = mGas->GetNumberOfComponents();
	vector< pair<string,double> > composition;
	for(int i=0; i<nComponents; i++) {
		double fraction;
		string label;
		mGas->GetComponent(i,label,fraction);
		composition.push_back( make_pair(label,fraction) );
	}

	for(vector< pair<string,double> >::iterator it = composition.begin(); it != composition.end(); it++)	
		gasName += it->first + "-" + toString(it->second) + "_";
	
	return gasName;
}

void TDetector::setGasMixture() {
	if (bGasLoaded){
		cerr << "TDetector::setGasMixture -- Error, gas already set" << endl;
		return;
	}
	
	mGas = new Garfield::MediumMagboltz();
	
	switch (fConfig.nGases){
		case (1): 
			mGas->SetComposition(fConfig.gasNames[0], fConfig.gasPercentage[0]);
			break;
		case (2):
			mGas->SetComposition(fConfig.gasNames[0], fConfig.gasPercentage[0], fConfig.gasNames[1], fConfig.gasPercentage[1]);
			break;
		case (3):
			mGas->SetComposition(fConfig.gasNames[0], fConfig.gasPercentage[0], fConfig.gasNames[1], fConfig.gasPercentage[1], fConfig.gasNames[2], fConfig.gasPercentage[2]);
			break;
		case (4):
			mGas->SetComposition(fConfig.gasNames[0], fConfig.gasPercentage[0], fConfig.gasNames[1], fConfig.gasPercentage[1], fConfig.gasNames[2], fConfig.gasPercentage[2], fConfig.gasNames[3], fConfig.gasPercentage[3]);
			break;
		case (5):
			mGas->SetComposition(fConfig.gasNames[0], fConfig.gasPercentage[0], fConfig.gasNames[1], fConfig.gasPercentage[1], fConfig.gasNames[2], fConfig.gasPercentage[2], fConfig.gasNames[3], fConfig.gasPercentage[3], fConfig.gasNames[4], fConfig.gasPercentage[4]);
			break;
		case (6):
			mGas->SetComposition(fConfig.gasNames[0], fConfig.gasPercentage[0], fConfig.gasNames[1], fConfig.gasPercentage[1], fConfig.gasNames[2], fConfig.gasPercentage[2], fConfig.gasNames[3], fConfig.gasPercentage[3], fConfig.gasNames[4], fConfig.gasPercentage[4], fConfig.gasNames[5], fConfig.gasPercentage[5]);
			break;
	}
	mGas->SetTemperature(fConfig.gasTemperature);
	mGas->SetPressure(fConfig.gasPressure);
	
	mGasTableName = getGasName();
	mGasTableName += "temp-" + toString(mGas->GetTemperature()) + "_pres-" + toString(mGas->GetPressure()) + ".gas";
	
	cout << "\tGas table file: " << mGasTableName << endl;
	
	if(!file_exist("gastables/"+mGasTableName))	
		makeGasTable();
	mGas->LoadGasFile("gastables/"+mGasTableName);
	//gas->EnableDebugging();

	fTemperature = mGas->GetTemperature();
	fPressure = mGas->GetPressure();
	
	writeGasTransportParameters();
	
	bGasLoaded = true;
}

void TDetector::setGasMixture(Garfield::MediumMagboltz* gas){
//	if (bGasLoaded){
//		cerr << "TDetector::setGasMixture -- Error, gas already set" << endl;
//		return;
//	}
	
	mGas = gas;
	mGasTableName = getGasName();
	mGasTableName += "temp-" + toString(mGas->GetTemperature()) + "_pres-" + toString(mGas->GetPressure()) + ".gas";
	
	cout << "\tGas table file: " << mGasTableName << endl;
	
	if(!file_exist("gastables/"+mGasTableName))	
		makeGasTable();
	mGas->LoadGasFile("gastables/"+mGasTableName);
	//gas->EnableDebugging();

	fTemperature = mGas->GetTemperature();
	fPressure = mGas->GetPressure();
	
	bGasLoaded = true;
}

void TDetector::setElectricField(const double& Ex, const double& Ey, const double& Ez){
	fElectricField[0] = Ex;
	fElectricField[1] = Ey;
	fElectricField[2] = Ez;
}

Garfield::MediumMagboltz* TDetector::getGas(){
	return mGas;
}

void TDetector::makeGasTable(){
	cout << "Generating gas table -- could take some time" << endl;
	Garfield::MediumMagboltz* gas = mGas;
	gas->SetFieldGrid(100., 120000., 50, false, 0., 0., 1, 0., 0., 1);
	gas->EnableDebugging();
	gas->Initialise();
	gas->DisableDebugging();
	//const double lP = 0.0;
	//const double rP = 0.0;
	//gas->EnablePenningTransfer(rP, lP, "c2h2f4");
	gas->GenerateGasTable(4, true);

	gas->WriteGasFile("gastables/"+mGasTableName);
}

void TDetector::initialiseDetector(){

	fDx = fConfig.gapWidth/iNstep;
	double cx = fConfig.gapWidth/2., cy = 0, cz = 0, lx = fConfig.gapWidth/2., ly = 10, lz = 10;
	// Detector geometry
	// Gap [cm]
	Garfield::SolidBox* box = new Garfield::SolidBox(cx, cy, cz, lx, ly, lz);
	Garfield::GeometrySimple* geo = new Garfield::GeometrySimple();
	geo->AddSolid(box, mGas);

	// Make a component
	Garfield::ComponentConstant* comp = new Garfield::ComponentConstant();
	comp->SetGeometry(geo);
	comp->SetElectricField(fElectricField[0], fElectricField[1], fElectricField[2]);

	// Make a sensor
	mSensor = new Garfield::Sensor();
	mSensor->AddComponent(comp);

	// Get transport parameters
	mGas->ElectronTownsend(fElectricField[0],fElectricField[1],fElectricField[2],0.,0.,0., fAlpha);
	mGas->ElectronAttachment(fElectricField[0],fElectricField[1],fElectricField[2],0.,0.,0., fEta);
	mGas->ElectronVelocity(fElectricField[0],fElectricField[1],fElectricField[2],0.,0.,0.,fVx,fVy,fVz);
	mGas->IonVelocity(fElectricField[0],fElectricField[1],fElectricField[2],0.,0.,0.,fiVx,fiVy,fiVz);
	mGas->ElectronDiffusion(fElectricField[0],fElectricField[1],fElectricField[2],0.,0.,0.,fDiffL,fDiffT);

	// Multiply by -1 in order to have the electron drift velocity going forward on the x-axis
	fVx *= -1., fVy *= -1, fVz *= -1;

	k = fEta/fAlpha;

	fDt = fDx/fVx; //ns
	
	// If we simulate avalanche, check if alpha and eta are greater than 0 (otherwise no point in carrying on) ...
	if (!fConfig.noAvalanche)
		assert(  (fAlpha>0) and (fEta>0) );
	
	cout << endl;
	cout << "Transport parameters:" << endl;
	cout << "\talpha: " << fAlpha << " cm-1" << endl;
	cout << "\teta: " << fEta << " cm-1" << endl;
	cout << "\tE: " << fElectricField[0] << " kV/cm" << endl;
	cout << "\tDrift velocity: (" << fVx << "," << fVy << "," << fVz << ")" << endl;
	cout << "\tIon Drift velocity: (" << fiVx << "," << fiVy << "," << fiVz << ")" << endl;
	cout << "\tDiffusion coefficient: (" << fDiffL << ", " << fDiffT << ")" << endl;
	cout << "\tDt: " << fDt << endl;
	cout << endl;
	
	/*delete box;
	delete geo;
	delete comp;*/
	
	bDetectorInitialised = true;
	plotSC();
	//cout << (1 or 0) << endl;
	if( !(fConfig.noAvalanche or fConfig.onlyMult) )
		makeEbarTable();
}

vector<double> TDetector::getTransportParameters(double Ex, double Ey, double Ez){
	//double alpha, eta, vx, vy, vz;
	vector<double> parameters (5,-1);
	mGas->ElectronTownsend(Ex,Ey,Ez,0.,0.,0., parameters.at(0));
	mGas->ElectronAttachment(Ex,Ey,Ez,0.,0.,0., parameters.at(1));
	mGas->ElectronVelocity(Ex,Ey,Ez,0.,0.,0., parameters.at(2), parameters.at(3), parameters.at(4));

	return parameters;

}

vector<double> TDetector::getDiffusionCoefficients( double const& Ex, double const& Ey, double const& Ez ) {
	vector<double> coeff (2,-1);
	mGas->ElectronDiffusion(Ex,Ey,Ez,0.,0.,0.,coeff[0],coeff[1]);
	return coeff;
}

void TDetector::writeGasTransportParameters(){
	if (!bGasLoaded) {
		cerr << "TDetector::writeGasTransportParameters -- No gas table found/loaded." << endl;
		return;
	}
	
	ofstream data(("out/"+mGasTableName+".dat").c_str(), ios::out | ios::trunc);
	data << "#eta	alpha	v	dl	dt	E" << endl;
	double alpha, eta, vx, vy, vz, dl, dt;
	double E = 100.;
	for (int i=0; i<50; i++){
		E += (120000. - 100.)/50.;
		mGas->ElectronTownsend(E,0.,0.,0.,0.,0., alpha);
		mGas->ElectronAttachment(E,0.,0.,0.,0.,0., eta);
		mGas->ElectronVelocity(E,0,0,0.,0.,0., vx,vy,vz);
		mGas->ElectronDiffusion(E,0.,0.,0.,0.,0.,dl,dt);
		data << eta << "\t" << alpha << "\t" << vx << "\t" << dl << "\t" << dt << "\t" << E << endl;
	}
	
	data.close();
	
}

void TDetector::setGarfieldSeed( const int& s ) {
	cout.setstate(ios_base::failbit);
	Garfield::randomEngine.Seed(s);
	cout.clear();
	
	if (fConfig.verbose)
		cout << "Garfield TRandom3 seed: " << s << endl;
}

double TDetector::R(const double& k, const double& z, const double& zp){
	/*	q,g and p have to be in mm in order to have correct values (From C. Lippmann, PhD Thesis, 2003, figure 2.11)
	*	everything is defined in cm in the config */
	//double q = fGeometry.resistiveLayersWidth[0] * Constants::cm; //cathode
	double q = fConfig.cathodeWidth *10; //cathode
	double g = fConfig.gapWidth *10;
	double p = g + fConfig.anodeWidth *10; //anode
	double eps0 = Constants::VacuumPermittivity; //GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	double eps1 = fConfig.cathodePermittivity * eps0; //cathode
	double eps3 = fConfig.anodePermittivity * eps0;
	double eps2 = eps0;

    return (eps1+eps2)*(eps1+eps2) * (eps2+eps3)*(eps2+eps3) * (exp(k*(-2*p-2*q+z-zp)) + exp(k*(-2*p-2*q-z+zp)))
            - (eps1+eps2)*(eps1+eps2) * (eps2-eps3)*(eps2-eps3) * exp(k*(-4*g-2*q+z+zp))
            - 4*eps1*eps2*(eps2+eps3)*(eps2+eps3) * exp(k*(-2*q-z-zp)) - (eps1-eps2)*(eps1-eps2) * (eps2+eps3)*(eps2+eps3) * exp(k*(-2*p-z-zp))
            - (eps1*eps1 - eps2*eps2)*(eps2-eps3)*(eps2-eps3) * exp(k*(-4*g+z+zp))
            + (eps1*eps1 - eps2*eps2)*(eps2+eps3)*(eps2+eps3) * ( -exp(k*(-2*p-2*q-z-zp)) + exp(k*(-2*p+z-zp)) + exp(k*(-2*p-z+zp)) )
            - 4*(eps1*eps1 - eps2*eps2)*eps2*eps3*exp(k*(-2*p-2*q+z+zp)) - 4*(eps1+eps2)*(eps1+eps2) * eps2*eps3*exp(k*(-2*p+z+zp))
            + (eps1-eps2)*(eps1-eps2) * (eps2*eps2 - eps3*eps3) * exp(k*(-2*g-z-zp)) + 4*eps1*eps2*(eps2*eps2 - eps3*eps3)*exp(k*(2*g-2*p-2*q-z-zp))
            + (eps1+eps2)*(eps1+eps2) * (eps2*eps2 - eps3*eps3) * ( -exp(k*(-2*g-2*q+z-zp)) - exp(k*(-2*g-2*q-z+zp)) + exp(k*(-2*g-2*p-2*q+z+zp)) )
            + (eps1*eps1 - eps2*eps2)*(eps2*eps2 - eps3*eps3) * ( exp(k*(-2*g-2*q-z-zp)) - exp(k*(-2*g+z-zp)) - exp(k*(-2*g-z+zp)) + exp(k*(-2*g-2*p+z+zp)) );
}

double TDetector::D(const double& k){
	/*	q,g and p have to be in mm in order to have correct values (From C. Lippmann, PhD Thesis, 2003, figure 2.11) 
	*	everything is defined in cm in the config */
	double q = fConfig.cathodeWidth	*10; //cathode
	double g = fConfig.gapWidth *10;
	double p = g + fConfig.anodeWidth *10; //anode
	double eps0 = Constants::VacuumPermittivity;
	double eps1 = fConfig.cathodePermittivity * eps0; //cathode
	double eps3 = fConfig.anodePermittivity * eps0;
	double eps2 = eps0;

    return (eps1+eps2)*(eps2+eps3)*(1-exp(-2*k*(p+q)))
            - (eps1-eps2)*(eps2+eps3)*(exp(-2*k*p)-exp(-2*k*q))
            - (eps1+eps2)*(eps2-eps3)*(exp(-2*k*(p-g))-exp(-2*k*(q+g)))
            + (eps1-eps2)*(eps2-eps3)*(exp(-2*k*g)-exp(-2*k*(p+q-g)));
}

double integrand(double x, void * params){
	//x == kappa
	/* Everything here has to be in mm ! */ /* So it seems ... */
	double* param = reinterpret_cast<double*> (params); //[P,z,zp]
	//cout << param[0] << " " << param[1] << " " << param[2] << endl;
	//return gsl_sf_bessel_J0(x*param[0]) * tgsl->R(x,param[1],param[2])/tgsl->D(x);
	return bessel_J0(x*param[0]) * tgsl->R(x,param[1],param[2])/tgsl->D(x);
}

double TDetector::SCFieldSimplified(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp) {
	double e0 = Constants::ElectronCharge; //GSL_CONST_MKSA_ELECTRON_CHARGE;
	double g = fConfig.gapWidth * Constants::cm;
	double eps0 = Constants::VacuumPermittivity; //GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	double eps1 = fConfig.cathodePermittivity * eps0; //cathode
	double eps3 = fConfig.anodePermittivity * eps0;
	double eps2 = eps0;

	double P2 = r*r - 2*r*rp*cos(phi-phip) + rp*rp;

	return ( e0/(4*Constants::Pi*eps2) ) * (   ( (z-zp)/(pow(P2 + (z-zp)*(z-zp),1.5)) )  -  ( ((eps2-eps3)/(eps2+eps3))*(2*g-z-zp)/(pow(P2+(2*g-z-zp)*(2*g-z-zp),1.5)) )     );

}

double TDetector::SCFieldSimplified2(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp) {
	double e0 = Constants::ElectronCharge; //GSL_CONST_MKSA_ELECTRON_CHARGE;
	double g = fConfig.gapWidth * Constants::cm;
	double eps0 = Constants::VacuumPermittivity; //GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	double eps1 = fConfig.cathodePermittivity * eps0; //cathode
	double eps3 = fConfig.anodePermittivity * eps0;
	double eps2 = eps0;

	double P2 = r*r - 2*r*rp*cos(phi-phip) + rp*rp;

	return ( e0/(4*Constants::Pi*eps2) ) * (   ( (z-zp)/(pow(P2 + (z-zp)*(z-zp),1.5)) )  -  ( ((eps2-eps3)/(eps2+eps3))*(2*g-z-zp)/(pow(P2+(2*g-z-zp)*(2*g-z-zp),1.5)) ) - ( ((eps1-eps2)/(eps1+eps2))*(z+zp)/(pow(P2+(z+zp)*(z+zp),1.5)) )     );

}

double TDetector::SCField(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp) {
	double h = 0.00001 * Constants::mm;
	return -1. * (SCPotential(r,phi,z+h,rp,phip,zp)-SCPotential(r,phi,z,rp,phip,zp))/h;
	
	//return -1. * ( SCPotential(r,phi,z+h,rp,phip,zp) - SCPotential(r,phi,z-h,rp,phip,zp) ) / (2*h);
	/*double m1 = ( SCPotential(r,phi,z+h,rp,phip,zp) - SCPotential(r,phi,z-h,rp,phip,zp) ) / (2*h);
	double m2 = ( SCPotential(r,phi,z+2*h,rp,phip,zp) - SCPotential(r,phi,z-2*h,rp,phip,zp) ) / (4*h);
	double m3 = ( SCPotential(r,phi,z+3*h,rp,phip,zp) - SCPotential(r,phi,z-3*h,rp,phip,zp) ) / (6*h);

	double deriv = (1.5*m1 - 0.6*m2 + 0.1*m3);

	return -1. * deriv;*/
}

double TDetector::SCPotential(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp){
	/*	Ok! So in order to have a proper result, the parameters has to be parsed in m !	
	 *	Then P, z and zp are parsed to the integrand function for correction term in mm! */
	if( tgsl != this )
		tgsl = this;

    double P = sqrt(r*r - 2*r*rp*cos(phi-phip) + rp*rp);

  //  gsl_integration_workspace * w = gsl_integration_workspace_alloc (4000);
   // double result, error;
   
   /*	Parameters have to be converted in mm	*/
    double funParams[3] = {P*1000,z*1000,zp*1000};
   /* gsl_function F;
	F.function = &integrand;
	F.params = &funParams;
	gsl_integration_qagiu (&F, 0., 1.e-2, 1e-7, 4000, w, &result, &error);
	gsl_integration_workspace_free (w);*/
	
	Rosetta::GaussLegendreQuadrature<5> gl5;
	double corr = gl5.integrate_iu(0, integrand, funParams);
	double corrTerm = 0.;
	
	//ofstream data("out/corrTerm.dat", ios::out | std::ios::app);

    double Q = Constants::ElectronCharge;
	double eps0 = Constants::VacuumPermittivity;
    double eps1 = fConfig.cathodePermittivity * eps0; //cathode
	double eps3 = fConfig.anodePermittivity * eps0;
	double eps2 = eps0;
	double g = fConfig.gapWidth * Constants::cm;
	corrTerm = ( Q/(4*Constants::Pi*eps2) ) *  (1./((eps1+eps2)*(eps2+eps3))) * corr;
	
	double pot1,pot2,pot3 = 0.;
	pot1 = ( Q/(4*Constants::Pi*eps2) ) * (1./(sqrt(P*P+(z-zp)*(z-zp))));	/* Free charge term */
	pot2 = ( Q/(4*Constants::Pi*eps2) ) * ((eps1-eps2)/((eps1+eps2)*sqrt(P*P+(z+zp)*(z+zp))));	/* Mirror charge (-z') term */
	pot3 = ( Q/(4*Constants::Pi*eps2) ) * ((eps3-eps2)/((eps3+eps2)*sqrt(P*P+(2*g-z-zp)*(2*g-z-zp))));	/* Mirror charge (g-z') term */
    //double pot = ( Q/(4*Constants::Pi*eps2) ) * ( (1./(sqrt(P*P+(z-zp)*(z-zp)))) - ((eps1-eps2)/((eps1+eps2)*sqrt(P*P+(z+zp)*(z+zp)))) - ((eps3-eps2)/((eps3+eps2)*sqrt(P*P+(2*g-z-zp)*(2*g-z-zp)))) +
                                //  (1./((eps1+eps2)*(eps2+eps3))) * 0. );
	
	//cout << corrTerm << endl;
	
	/*	The correction term needs to be multplied by 1e3 to output correct value ... Why ? that's a mystery ...	*/
	//data << pot1 << "\t" << -pot2 << "\t" << -pot3 << "\t" << corrTerm*1000. << endl;
	//data.close();
    return pot1 - pot2 - pot3 + corrTerm*1000.;
}

inline double TDetector::RadialChargeDistribution(const double& r, const double& l){
	return ( 1./(fDiffT*fDiffT * l) ) * exp( -(r*r)/(2*fDiffT*fDiffT * l) ) ;
}

double Ebar(double x, void * params){
	// x == rp

	double* param = reinterpret_cast<double*> (params); //[z,l,zp]
	//cout << x << "\t" << param[0] << "\t" << param[1] << "\t" << param[2] << endl;
	/*	RadialChargeDistribution return a value in cm^2. SCField a value in V/m.	*/
	double res =  tgsl->RadialChargeDistribution(x,param[1]) * tgsl->SCFieldSimplified2(0.,0.,param[0],x,0.,param[2])*0.01 * x;
	return res;
}

double TDetector::computeEbar(const double& z, const double& l, const double& zp){
	/*	Z and Zp in m. L in cm.	*/
	if( tgsl != this )
		tgsl = this;
	
	bool gl = false;
	
	if (gl) {
		double funParams[3] = {z,l,zp};
		Rosetta::GaussLegendreQuadrature<90> gl5;
		//cout << z << "\t" << zp << "\t" << l << endl;
		return gl5.integrate_iu(0., Ebar, funParams);
	}
	
	else {
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
		double result, error;
		double funParams[3] = {z,l,zp};
		gsl_function F;
		F.function = &Ebar;
		F.params = &funParams;

		//gsl_integration_qag (&F, 0., 5*sqrt(sigma2), 1e-10, 1e-7, 3000, 3, w, &result, &error);

		gsl_integration_qagiu (&F, 0., 1.e-4, 1e-7, 10000, w, &result, &error);
		gsl_integration_workspace_free (w);

		return result;
	}
}

#if defined(PYTHON)
	double TDetector::computeEbar_Python(const double& z, const double& l, const double& zp){
		// Very slow !!!!
		//double mm = 1.e-3;
		
		double eps0 = Constants::VacuumPermittivity; //GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	    double eps1 = fConfig.cathodePermittivity * eps0; //cathode
		double eps3 = fConfig.anodePermittivity * eps0;
		double eps2 = eps0;
		
		double values[10] = {z, l, zp, fDiffT, eps1,eps2,eps3,4.*Constants::mm,2.*Constants::mm,2.*Constants::mm};
		vector<double> args (values, values + sizeof(values) / sizeof(values[0]) );
		double Ebar;
		call_python_fun("compute_Ebar", args, Ebar);
		
		return Ebar;
	}
#endif

void TDetector::makeEbarTable( bool const& binary ){
	// WARNING if table was written in plain text, the loading from binary format will be incorrect!!!!
	// WARNING writing and loading has to be done in the same format (binary or plain text)
	
	if( !bDetectorInitialised ){
		cerr << "Error -- TDetector::makeEbarTable -- Detector needs to be initialised first." << endl;
		exit(0);
	}
	
	if(bHasEbarTable)	return;

	int n = iEbarTableSize+1;
	int size = (n)*(n)*(n);
	
	string fileName = getUniqueTableName();
	cout << "FileName: " << fileName << endl;

	hash<string> hash;
	string hexFileName = toString( hash(fileName) );
	cout << "Hash representation: " << hexFileName << endl;
	
	/*	z and zp in m, l in cm.	*/
	double zStep = fConfig.gapWidth * Constants::cm / (n-1), zpStep = fConfig.gapWidth * Constants::cm / (n-1), lStep = iNstep*fDx / (n);
	fEbarLarray = vector<double>(n,0.);
	fEbarZarray = vector<double>(n,0.);
	fEbarZparray = vector<double>(n,0.);

	for(int i=0; i<n; i++){
		fEbarZarray[i] = (i)*zStep;
		fEbarZparray[i] = (i)*zpStep;
		fEbarLarray[i] = (i+1)*lStep;
	}

	fEbarVecTable = vector<double> (size);

	if (file_exist("EbarTables/"+hexFileName)){
		ifstream inf;
		
		if ( binary ) {
			inf.open(("EbarTables/"+hexFileName).c_str(),ios::binary);
			inf.read(reinterpret_cast<char*>(&fEbarVecTable[0]), fEbarVecTable.size() * sizeof(fEbarVecTable[0]));
		}
		else {
			inf.open(("EbarTables/"+hexFileName).c_str());
			double z,zp,l,Ebar;
			int i=0;
			while (inf >> z >> zp >> l >> Ebar){
				fEbarVecTable[i] = Ebar;
				i++;
			}
		}
		//inf.read(reinterpret_cast<char*>(&fEbarVecTable[0]), fEbarVecTable.size() * sizeof(fEbarVecTable[0]));
		
		if ( std::accumulate(fEbarVecTable.begin(),fEbarVecTable.end(),0.) == 0 ) {
			cerr << "Error while loading Ebar table." << endl;
			exit(0);
		}
		
		fEbarTableHexName = hexFileName;
		bHasEbarTable = true;
		printEbarTable();
		return;
	}

	cout << "Computing Ebar Table" << endl;

	ofstream data("out/Ebar.dat", ios::out | ios::trunc);

	int i,j,k;
	//#pragma omp parallel for private(j,k) collapse(3)
	for(i=0; i<n; i++){	//z
		for(j=0; j<n; j++){	//zp
			for(k=0; k<n; k++){	//l
				double Ebar = computeEbar( fEbarZarray[i], fEbarLarray[k], fEbarZparray[j] );
				cout << fEbarZarray[i] << "\t" << fEbarZparray[j] << "\t" << fEbarLarray[k] << "\t" << Ebar << endl;
				fEbarVecTable[ (long)i*(long)n*(long)n + (long)j*(long)n + (long)k ] = Ebar;
			}
		}
	}
	
	ofstream o;
	if (binary) {
		o.open(("EbarTables/"+hexFileName).c_str(),ios::binary);
		const char* pointer = reinterpret_cast<const char*>(&fEbarVecTable[0]);
		o.write( pointer, sizeof(fEbarVecTable[0])*fEbarVecTable.size() );
	}
	else {
		o.open(("EbarTables/"+hexFileName).c_str(), ios::out | ios::trunc);
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
				for(int k=0; k<n; k++){
					data << (i)*zStep << "\t" << (j)*zpStep << "\t" << (k+1)*lStep << "\t" << fEbarVecTable[ (long)i*(long)n*(long)n + (long)j*(long)n + (long)k ] << endl;
					o << fEbarZarray[i] << "\t" << fEbarZparray[j] << "\t" << fEbarLarray[k] << "\t" << fEbarVecTable[ (long)i*(long)n*(long)n + (long)j*(long)n + (long)k ] << endl;
				}
			}
		}
	}

	fEbarTableHexName = hexFileName;

	o.close();
	data.close();
	bHasEbarTable = true;
}

void TDetector::printEbarTable() {
	if (!bHasEbarTable) {
		cerr << "Error TDetector::printEbarTable -- No Ebar table found" << endl;
		return;
	}
	
	int n = iEbarTableSize+1;
	int size = n*n*n;
	ofstream data("out/Ebar.dat", ios::out | ios::trunc);
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			for(int k=0; k<n; k++){
				data << fEbarZarray[i] << "\t" << fEbarZparray[j] << "\t" << fEbarLarray[k] << "\t" << fEbarVecTable[ (long)i*(long)n*(long)n + (long)j*(long)n + (long)k ] << endl;
			}
		}
	}
	data.close();
}

void TDetector::plotSC(){
	if( tgsl != this )
		tgsl = this;
		
	//cout << gsl_sf_bessel_J0(0) << " " << gsl_sf_bessel_J0(1) << " " << gsl_sf_bessel_J0(2) << endl;
	//cout << BESSJ0(0) << " " << BESSJ0(1) << " " << BESSJ0(2) << endl;
	//exit(0);
	cout << "plot" << endl;
	//double mm = 0.001;
	cout << fConfig.gapWidth << endl;
	double min = -1*Constants::mm, max = 1*Constants::mm;
	double h = 0.01*Constants::mm;
	vector<double> r(1,min);
	vector<double> z(1,0);
	while(r.back() <= max)	r.push_back( r.back()+h );
	while(z.back() <= fConfig.gapWidth)	z.push_back( z.back() + h );
	ofstream data("out/plotSC.dat", ios::out | ios::trunc);
	ofstream data2("out/plotI.dat", ios::out | ios::trunc);
	ofstream data3("out/plotSC2d.dat", ios::out | ios::trunc);

//	for(uint i=0; i<r.size(); i++){
//		for(uint j=0; j<z.size(); j++) {
//				data << r[i]*1e3 << "\t" << z[j]*1e3 << "\t" << SCPotential(r[i],0,z[j],0,0,6*Constants::cm) << endl;
//		}
//	}

/*	
	double k = 0;
	double pars[3] = {0.,0.6*0.001,0.3*0.001};
	while(k<3){
		data2 << k << "\t" << integrand(k,pars) << endl;
		k+=0.01;
	}

	for(uint i=0; i<z.size(); i++)
		data3 << z[i] << "\t" << SCField(0.,0,z[i]*0.01,0,0,0.4*0.001) << "\t" << SCFieldSimplified2(0.,0,z[i]*0.01,0,0,0.4*0.001) << "\t" << SCFieldSimplified(0.,0,z[i]*0.01,0,0,0.4*0.001) << endl;
	
	data.close();
	data2.close();
	data3.close();
	*/
	//cout << computeEbar(0.00056,	0.06,	0.0012) << endl;
	makeEbarTable();
}

string TDetector::getUniqueTableName(){
	string name = "fDiffT:"+toString(fDiffT) + "-gapWidth:"+toString(fConfig.gapWidth) + "-eps1:"+toString(fConfig.cathodePermittivity) + "-eps3:"+toString(fConfig.anodePermittivity)
	+ "-cathode:"+toString(fConfig.cathodeWidth) + "-anode:"+toString(fConfig.anodeWidth)
	+ "-Nstep:"+toString(iNstep) + "-Dx:"+toString(fDx) + "-EbarTableSize:"+toString(iEbarTableSize);
	
	return name;
}
