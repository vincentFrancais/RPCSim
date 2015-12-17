
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <assert.h> 
#include <utility>
#include <cstring>
#include <fcntl.h>

#include "TDetector.hpp"
#include "helper_functions.hpp"
#include "gsl/gsl_const_mksa.h"
#include "gsl/gsl_sf_bessel.h"
#include <gsl/gsl_integration.h>
#include "gsl/gsl_math.h"

#define PI 3.14159265358979323846


using namespace std;
static TDetector* tgsl = 0;

TDetector::TDetector(const DetectorGeometry& geometry, const int& nStep){
	iNstep = nStep;
	//fGapWidth = gasGap;
	fGeometry = geometry;
	
	bEbarComputed = false;
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
	
	for(vector< pair<string,double> >::iterator it = composition.begin(); it != composition.end(); it++)	mGasTableName += it->first + "-" + to_string(it->second) + "_";
	mGasTableName += "temp-" + to_string(mGas->GetTemperature()) + "_pres-" + to_string(mGas->GetPressure()) + ".gas";
	
	if(!file_exist("gastables/"+mGasTableName))	makeGasTable();
	mGas->LoadGasFile("gastables/"+mGasTableName);
	//gas->EnableDebugging();
	
	fTemperature = mGas->GetTemperature();
	fPressure = mGas->GetPressure();
	
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
	gas->SetFieldGrid(100., 120000., 70, false, 0., 0., 1, 0., 0., 1);
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
	
	assert( (fAlpha>fEta) or (fAlpha>0) );
	
	cout << "Transport parameters:" << endl;
	cout << "\talpha: " << fAlpha << " cm-1" << endl;
	cout << "\teta: " << fEta << " cm-1" << endl;
	cout << "\tE: " << fElectricField[0] << " kV/cm" << endl;
	cout << "\tDrift velocity: (" << fVx << "," << fVy << "," << fVz << ")" << endl;
	cout << "\tIon Drift velocity: (" << fiVx << "," << fiVy << "," << fiVz << ")" << endl;
	cout << "\tDiffusion coefficient: (" << fDiffL << ", " << fDiffT << ")" << endl;
	cout << "\tDt: " << fDt << endl;
	
	makeEbarTable();
}

double* TDetector::getTransportParameters(double Ex, double Ey, double Ez){
	//double alpha, eta, vx, vy, vz;
	double* parameters = new double[5];
	mGas->ElectronTownsend(Ex,Ey,Ez,0.,0.,0., parameters[0]);
	mGas->ElectronAttachment(Ex,Ey,Ez,0.,0.,0., parameters[1]);
	mGas->ElectronVelocity(Ex,Ey,Ez,0.,0.,0., parameters[2], parameters[3], parameters[4]);
	
	return parameters;
	
}

double TDetector::R(const double& k, const double& z, const double& zp){
	double cm = 0.01;
	
	double q = fGeometry.resistiveLayersWidth[0] * cm;
	double g = fGeometry.gapWidth * cm;
	double p = (g + fGeometry.resistiveLayersWidth[1]) * cm;
	double eps1 = fGeometry.relativePermittivity[0];
	double eps3 = eps1;
	double eps2 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	
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
	double cm = 0.01;
	
	double q = fGeometry.resistiveLayersWidth[0] * cm;
	double g = fGeometry.gapWidth * cm;
	double p = (g + fGeometry.resistiveLayersWidth[1]) * cm;
	double eps1 = fGeometry.relativePermittivity[0];
	double eps3 = eps1;
	double eps2 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	
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

double TDetector::SCFieldSimplified(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp){
	double cm = 0.01;
	
	double e0 = GSL_CONST_MKSA_ELECTRON_CHARGE;
	double g = fGeometry.gapWidth * cm;
	double eps1 = fGeometry.relativePermittivity[0];
	double eps3 = eps1;
	double eps2 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	
	double P2 = r*r - 2*r*rp*cos(phi-phip) + rp*rp;
	
	return ( e0/(4*M_PI*eps2) ) * (   ( (z-zp)/(pow(P2 + (z-zp)*(z-zp),1.5)) )  -  ( ((eps2-eps3)/(eps2+eps3))*(2*g-z-zp)/(pow(P2+(2*g-z-zp)*(2*g-z-zp),1.5)) )     );
	
}

double TDetector::SCField(const double& r, const double& phi, const double& z, const double& rp, const double& phip, const double& zp){
	double mm = 0.001;
	double cm = 0.01;
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
	
	double cm = 0.01;
    double P = sqrt(r*r - 2*r*rp*cos(phi-phip) + rp*rp);
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (4000);
    double result, error;
    double funParams[3] = {P,z,zp};
    gsl_function F;
	F.function = &integrand;
	F.params = &funParams;
	gsl_integration_qagiu (&F, 0., 1e-4, 1e-7, 4000, w, &result, &error);
	gsl_integration_workspace_free (w);
    
    double Q = GSL_CONST_MKSA_ELECTRON_CHARGE;
    double eps1 = fGeometry.relativePermittivity[0];
	double eps3 = eps1;
	double eps2 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	double g = fGeometry.gapWidth * cm;
	
    double pot = ( Q/(4*M_PI*eps2) ) * ( (1./(sqrt(P*P+(z-zp)*(z-zp)))) - ((eps1-eps2)/((eps1+eps2)*sqrt(P*P+(z+zp)*(z+zp)))) - ((eps3-eps2)/((eps3+eps2)*sqrt(P*P+(2*g-z-zp)*(2*g-z-zp)))) + 
                                  (1./((eps1+eps2)*(eps2+eps3))) * result );
                            
    return pot;
}

double TDetector::RadialChargeDistribution(const double& r, const double& l){
	//cout << fDiffT << " " << l << " " << r << endl;
	return ( 1./(fDiffT*fDiffT * l) ) * exp( -(r*r)/(2*fDiffT*fDiffT * l) ) ;
}

double Ebar(double x, void * params){
	//x == rp
	double cm = 0.01;
	double* param = reinterpret_cast<double*> (params); //[z,l,zp]
	double res =  tgsl->RadialChargeDistribution(x,param[1]) * tgsl->SCField(0.,0.,param[0],x*cm,0.,param[2])*0.01 * x;  //RadialDistrib en cm-2, rp*drp en cm2, SCField en V/cm. rp doit etre en m dans les params de SCField
	return res;
}

double TDetector::computeEbar(const double& z, const double& l, const double& zp){
	if( tgsl != this )
		tgsl = this;
	
	double sigma2 = fDiffT*fDiffT * l;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (3000);
    double result, error;
    double funParams[3] = {z,l,zp};
    gsl_function F;
	F.function = &Ebar;
	F.params = &funParams;

	gsl_integration_qag (&F, 0., 3*sqrt(sigma2), 1e-10, 1e-7, 3000, 3, w, &result, &error);
	
	//gsl_integration_qagiu (&F, 0., 1.e-2, 1e-3, 15000, w, &result, &error);
	gsl_integration_workspace_free (w);
	
	return result;
}

void TDetector::makeEbarTable(){
	if(bEbarComputed)	return;
	
	iEbarTableSize = 100;
	int n = iEbarTableSize+1;
	int size = (n)*(n)*(n);
	
	string fileName = to_string(fDiffT)+to_string(fGeometry.gapWidth)+to_string(fGeometry.relativePermittivity[0])
	+to_string(fGeometry.relativePermittivity[1])+to_string(iNstep)+to_string(n)+to_string(fDx)+to_string(iEbarTableSize);
	
	cout << "FileName: " << fileName << endl;
	
	unsigned char* uc = new unsigned char[fileName.size()+1];
	memcpy(uc, fileName.c_str(), fileName.size());
	uc[fileName.size()]=0;
	string hexFileName = GetHexRepresentation(uc, fileName.size());
	cout << "Hex representation: " << hexFileName << endl;
	delete uc;

	//double table[size];
	
	double cm = 0.01;
	
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
		//ifstream inf(("EbarTables/"+hexFileName).c_str(),ios::binary);
		ifstream inf(("EbarTables/"+hexFileName).c_str());
		double z,zp,l,Ebar;
		int i=0;
		while (inf >> z >> zp >> l >> Ebar){
			fEbarVecTable[i] = Ebar;
			i++;
		}
		//inf.read((char*)&table,sizeof(table));
		//inf.read(reinterpret_cast<char*>(&fEbarVecTable[0]), fEbarVecTable.size() * sizeof(fEbarVecTable[0]));

		//fEbarVecTable = arrayToVec(table, size);

		fEbarTableHexName = hexFileName;
		
		bEbarComputed = true;
		return;
	}
	
	cout << "Computing Ebar Table" << endl;

	ofstream data("out/Ebar.dat", ios::out | ios::trunc);
	
	int i,j,k;
	#pragma omp parallel for private(j,k) collapse(3)
	for(i=0; i<n; i++){	//z
		for(j=0; j<n; j++){	//zp
			for(k=0; k<n; k++){	//l
				//double Ebar = computeEbar((i)*zStep,(k+1)*lStep,(j)*zpStep);
				double Ebar = computeEbar( fEbarZarray[i], fEbarLarray[k], fEbarZparray[j] );
				fEbarVecTable[ (long)i*(long)n*(long)n + (long)j*(long)n + (long)k ] = Ebar;
			}
		}
	}
	
	//ofstream o(("EbarTables/"+hexFileName).c_str(),ios::binary);
	ofstream o(("EbarTables/"+hexFileName).c_str(), ios::out | ios::trunc);
	//const char* pointer = reinterpret_cast<const char*>(&fEbarVecTable[0]);
	//o.write((char*)&fEbarVecTable[0], fEbarVecTable.size() * sizeof(fEbarVecTable[0]));
	//o.write((char*)&table,sizeof(table));
	//o.write((const char*)&fEbarVecTable, sizeof(fEbarVecTable));
	
	//fEbarVecTable = arrayToVec(table, size);
	//delete table;
	
	fEbarTableHexName = hexFileName;
 	
	for(int i=0; i<n; i++){	
		for(int j=0; j<n; j++){
			for(int k=0; k<n; k++){
				data << (i)*zStep << "\t" << (j)*zpStep << "\t" << (k+1)*lStep << "\t" << fEbarVecTable[ (long)i*(long)n*(long)n + (long)j*(long)n + (long)k ] << endl;
				o << fEbarZarray[i] << "\t" << fEbarZparray[j] << "\t" << fEbarLarray[k] << "\t" << fEbarVecTable[ (long)i*(long)n*(long)n + (long)j*(long)n + (long)k ] << endl;
			}
		}
	}
	
	o.close();
	data.close();
	bEbarComputed = true;
}

void TDetector::plotSC(){
	cout << "plot" << endl;
	double cm = 0.01;

	double min = -1, max = 1;
	vector<double> val(1,min);
	while(val.back() <= max)	val.push_back( val.back()+0.02 );
	ofstream data("out/plotSC.dat", ios::out | ios::trunc);
	
	for(uint i=0; i<val.size(); i++){
		data << val[i] << "\t" << SCFieldSimplified(0,0,0.05*cm,val[i]*cm,0,0.05*cm) << endl;
	}
	data.close();
}
