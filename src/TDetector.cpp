
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

#include "gsl/gsl_const_mksa.h"
#include "gsl/gsl_sf_bessel.h"
#include <gsl/gsl_integration.h>
#include "gsl/gsl_math.h"


using namespace std;
static TDetector* tgsl = 0;

extern double cm;

TDetector::TDetector(const DetectorGeometry& geometry, const int& nStep){
	iNstep = nStep;
	fGeometry = geometry;
	iEbarTableSize = 75;

	bHasEbarTable = false;
	bGasLoaded = false;
	bDetectorInitialised = false;
}

TDetector::~TDetector(){
	delete mGas;
	delete mSensor;
}

void TDetector::setGasMixture(Garfield::MediumMagboltz* gas){
	mGas = gas;
	int nComponents = mGas->GetNumberOfComponents();
	vector< pair<string,double> > composition;
	for(int i=0; i<nComponents; i++) {
		double fraction;
		string label;
		mGas->GetComponent(i,label,fraction);
		composition.push_back( make_pair(label,fraction) );
	}

	for(vector< pair<string,double> >::iterator it = composition.begin(); it != composition.end(); it++)	
		mGasTableName += it->first + "-" + toString(it->second) + "_";
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

	fDx = fGeometry.gapWidth/iNstep;
	double cx = fGeometry.gapWidth/2., cy = 0, cz = 0, lx = fGeometry.gapWidth/2., ly = 10, lz = 10;
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
}

double* TDetector::getTransportParameters(double Ex, double Ey, double Ez){
	//double alpha, eta, vx, vy, vz;
	double* parameters = new double[5];
	mGas->ElectronTownsend(Ex,Ey,Ez,0.,0.,0., parameters[0]);
	mGas->ElectronAttachment(Ex,Ey,Ez,0.,0.,0., parameters[1]);
	mGas->ElectronVelocity(Ex,Ey,Ez,0.,0.,0., parameters[2], parameters[3], parameters[4]);

	return parameters;

}

double* TDetector::getDiffusionCoefficients( double const& Ex, double const& Ey, double const& Ez ) {
	double* coeff = new double[2];
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
}

double TDetector::R(const double& k, const double& z, const double& zp){
	double q = fGeometry.resistiveLayersWidth[0] * cm; //cathode
	double g = fGeometry.gapWidth * cm;
	double p = (g + fGeometry.resistiveLayersWidth[1]) * cm; //anode
	double eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	double eps1 = fGeometry.relativePermittivity[0] * eps0;
	double eps3 = eps1;
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
	double q = fGeometry.resistiveLayersWidth[0] * cm;	// cathode
	double g = fGeometry.gapWidth * cm;
	double p = (g + fGeometry.resistiveLayersWidth[1]) * cm;	// anode

	double eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	double eps1 = fGeometry.relativePermittivity[0] * eps0;
	double eps3 = eps1;
	double eps2 = eps0;

    return (eps1+eps2)*(eps2+eps3)*(1-exp(-2*k*(p+q)))
            - (eps1-eps2)*(eps2+eps3)*(exp(-2*k*p)-exp(-2*k*q))
            - (eps1+eps2)*(eps2-eps3)*(exp(-2*k*(p-g))-exp(-2*k*(q+g)))
            + (eps1-eps2)*(eps2-eps3)*(exp(-2*k*g)-exp(-2*k*(p+q-g)));
}

double integrand(double x, void * params){
	//x == k
	double* param = reinterpret_cast<double*> (params); //[P,z,zp]
	return gsl_sf_bessel_J0(x*param[0]) * tgsl->R(x,param[1],param[2])/tgsl->D(x);
}

double TDetector::SCFieldSimplified(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp) {
	double e0 = GSL_CONST_MKSA_ELECTRON_CHARGE;
	double g = fGeometry.gapWidth * cm;
	double eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	double eps1 = fGeometry.relativePermittivity[0] * eps0;
	double eps3 = eps1;
	double eps2 = eps0;

	double P2 = r*r - 2*r*rp*cos(phi-phip) + rp*rp;

	return ( e0/(4*M_PI*eps2) ) * (   ( (z-zp)/(pow(P2 + (z-zp)*(z-zp),1.5)) )  -  ( ((eps2-eps3)/(eps2+eps3))*(2*g-z-zp)/(pow(P2+(2*g-z-zp)*(2*g-z-zp),1.5)) )     );

}

double TDetector::SCField(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp) {
	double h = 0.001 * cm;

	double m1 = ( SCPotential(r,phi,z+h,rp,phip,zp) - SCPotential(r,phi,z-h,rp,phip,zp) ) / (2*h);
	double m2 = ( SCPotential(r,phi,z+2*h,rp,phip,zp) - SCPotential(r,phi,z-2*h,rp,phip,zp) ) / (4*h);
	double m3 = ( SCPotential(r,phi,z+3*h,rp,phip,zp) - SCPotential(r,phi,z-3*h,rp,phip,zp) ) / (6*h);

	double deriv = (1.5*m1 - (3./5.)*m2 + 0.1*m3);

	return -1. * deriv;
}

double TDetector::SCPotential(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp){
	if( tgsl != this )
		tgsl = this;

    double P = sqrt(r*r - 2*r*rp*cos(phi-phip) + rp*rp);

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (4000);
    double result, error;
    double funParams[3] = {P,z,zp};
    gsl_function F;
	F.function = &integrand;
	F.params = &funParams;
	gsl_integration_qagiu (&F, 0., 1.e-2, 1e-7, 4000, w, &result, &error);
	gsl_integration_workspace_free (w);

    double Q = GSL_CONST_MKSA_ELECTRON_CHARGE;
	double eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
    double eps1 = fGeometry.relativePermittivity[0] * eps0;
	double eps3 = eps1;
	double eps2 = eps0;
	double g = fGeometry.gapWidth * cm;

    double pot = ( Q/(4*M_PI*eps2) ) * ( (1./(sqrt(P*P+(z-zp)*(z-zp)))) - ((eps1-eps2)/((eps1+eps2)*sqrt(P*P+(z+zp)*(z+zp)))) - ((eps3-eps2)/((eps3+eps2)*sqrt(P*P+(2*g-z-zp)*(2*g-z-zp)))) +
                                  (1./((eps1+eps2)*(eps2+eps3))) * result );

    return pot;
}

inline double TDetector::RadialChargeDistribution(const double& r, const double& l){
	return ( 1./(fDiffT*fDiffT * l) ) * exp( -(r*r)/(2*fDiffT*fDiffT * l) ) ;
}

double Ebar(double x, void * params){
	// x == rp

	double* param = reinterpret_cast<double*> (params); //[z,l,zp]
	// RadialDistrib en cm-2, rp*drp en cm2, SCField en V/cm. rp doit etre en m dans les params de SCField
	double res =  tgsl->RadialChargeDistribution(x,param[1]) * tgsl->SCField(0.,0.,param[0],x*cm,0.,param[2])*0.01 * x;
	return res;
}

double TDetector::computeEbar(const double& z, const double& l, const double& zp){
	if( tgsl != this )
		tgsl = this;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (3000);
    double result, error;
    double funParams[3] = {z,l,zp};
    gsl_function F;
	F.function = &Ebar;
	F.params = &funParams;

	//gsl_integration_qag (&F, 0., 5*sqrt(sigma2), 1e-10, 1e-7, 3000, 3, w, &result, &error);

	gsl_integration_qagiu (&F, 0., 1.e-1, 1e-7, 3000, w, &result, &error);
	gsl_integration_workspace_free (w);

	return result;
}

#if defined(PYTHON)
	double TDetector::computeEbar_Python(const double& z, const double& l, const double& zp){
		// Very slow !!!!
		double mm = 1.e-3;
		
		double eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	    double eps1 = fGeometry.relativePermittivity[0] * eps0;
		double eps3 = eps1;
		double eps2 = eps0;
		
		double values[10] = {z, l, zp, fDiffT, eps1,eps2,eps3,4.*mm,2.*mm,2.*mm};
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

	/*unsigned char* uc = new unsigned char[fileName.size()+1];
	memcpy(uc, fileName.c_str(), fileName.size());
	uc[fileName.size()]=0;
	string hexFileName = GetHexRepresentation(uc, fileName.size());
	cout << "Hex representation: " << hexFileName << endl;
	delete uc; */
	hash<string> hash;
	string hexFileName = toString( hash(fileName) );
	cout << "Hash representation: " << hexFileName << endl;
	
	double zStep = fGeometry.gapWidth * cm / (n-1), zpStep = fGeometry.gapWidth * cm / (n-1), lStep = iNstep*fDx / (n);
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
		inf.read(reinterpret_cast<char*>(&fEbarVecTable[0]), fEbarVecTable.size() * sizeof(fEbarVecTable[0]));
		
		if ( std::accumulate(fEbarVecTable.begin(),fEbarVecTable.end(),0.) == 0 ) {
			cerr << "Error while loading Ebar table." << endl;
			exit(0);
		}
		
		fEbarTableHexName = hexFileName;
		bHasEbarTable = true;
		return;
	}

	cout << "Computing Ebar Table" << endl;

	ofstream data("out/Ebar.dat", ios::out | ios::trunc);

	int i,j,k;
	#pragma omp parallel for private(j,k) collapse(3)
	for(i=0; i<n; i++){	//z
		for(j=0; j<n; j++){	//zp
			for(k=0; k<n; k++){	//l
				double Ebar = computeEbar( fEbarZarray[i], fEbarLarray[k], fEbarZparray[j] );
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

void TDetector::plotSC(){
	cout << "plot" << endl;
	double mm = 0.001;

	double min = -1, max = 1;
	double h = 0.05;
	vector<double> r(1,min);
	vector<double> z(1,0);
	while(r.back() <= max)	r.push_back( r.back()+h );
	while(z.back() <= fGeometry.gapWidth)	z.push_back( z.back() + h );
	ofstream data("out/plotSC.dat", ios::out | ios::trunc);

	for(uint i=0; i<r.size(); i++){
		for(uint j=0; j<z.size(); j++) {
				data << r[i] << "\t" << z[j] << "\t" << SCField(r[i]*mm,0,z[j],0,0,0.5*mm) << endl;
		}
	}
	data.close();
}

string TDetector::getUniqueTableName(){
	string name = "fDiffT:"+toString(fDiffT) + "-gapWidth:"+toString(fGeometry.gapWidth) + "-eps1:"+toString(fGeometry.relativePermittivity[0]) + "-eps3:"+toString(fGeometry.relativePermittivity[1])
	+ "-cathode:"+toString(fGeometry.resistiveLayersWidth[0]) + "-anode:"+toString(fGeometry.resistiveLayersWidth[1])
	+ "-Nstep:"+toString(iNstep) + "-Dx:"+toString(fDx) + "-EbarTableSize:"+toString(iEbarTableSize);
	
	return name;
}
