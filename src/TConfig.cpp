#include "TConfig.hpp"

using namespace std;
using namespace tinyxml2;

TConfig::TConfig (string file) {
	
	//TConfig config;
	nResisLayers = 0; 
	nGases = 0;
	outFile = "out/";
	generator = "SFMT";
	garfieldSeed = -1;
	//threshold = 0.; //pC
	verbose = false;
	verbosityLevel = 0;
	snaps = false;
	clusterDensity = false;
	electronProduction = false;
	noAvalanche = false;
	//computeEfficiency = false;
	onlyMult = false;
	debugOutput = false;
	
	XMLDocument doc;
	if ( doc.LoadFile( file.c_str() ) != XML_SUCCESS ){
		cerr << "TConfig -- Error opening file " << file << endl;
		exit(0);
	}
	
	XMLHandle docHandle( &doc );
	
	/*  Detector configuration */
	if(  docHandle.FirstChildElement("Detector").ToElement() == NULL){
		cerr << "TConfig -- Error reading config file\t";
		cerr << "Detector field missing in config file." << endl;
		cerr << XML_ERROR_PARSING_ELEMENT << endl;
		exit(0);
	}
	
	XMLElement* gapw = docHandle.FirstChildElement("Detector").FirstChildElement("GapWidth").ToElement();
	XMLElement* steps = docHandle.FirstChildElement("Detector").FirstChildElement("Steps").ToElement();
	XMLElement* efield = docHandle.FirstChildElement("Detector").FirstChildElement("ElectricField").ToElement();
	if ( gapw==NULL or steps==NULL or efield==NULL ){
		cerr << "TConfig -- Error reading config file\t";
		cerr << "GapWidth  or  Steps  or  ElectricField   field(s) missing in config file." << endl;
		exit(0);
	}
	gapw->QueryDoubleText( &gapWidth );
	steps->QueryIntText( &nSteps );
	efield->QueryDoubleText( &ElectricField );
	
	
	XMLElement* anodeWidthElem = docHandle.FirstChildElement( "Detector" ).FirstChildElement( "Anode" ).FirstChildElement( "width" ).ToElement();
	XMLElement* cathodeWidthElem = docHandle.FirstChildElement( "Detector" ).FirstChildElement( "Cathode" ).FirstChildElement( "width" ).ToElement();
	XMLElement* anodePermElem = docHandle.FirstChildElement( "Detector" ).FirstChildElement( "Anode" ).FirstChildElement( "permittivity" ).ToElement();
	XMLElement* cathodePermElem = docHandle.FirstChildElement( "Detector" ).FirstChildElement( "Cathode" ).FirstChildElement( "permittivity" ).ToElement();
	if ( anodeWidthElem and cathodeWidthElem and anodePermElem and cathodePermElem ) {
		anodeWidthElem->QueryDoubleText( &anodeWidth );
		anodePermElem->QueryDoubleText( &anodePermittivity );
		
		cathodeWidthElem->QueryDoubleText( &cathodeWidth );
		cathodePermElem->QueryDoubleText( &cathodePermittivity );
	} 
	else {
		cerr << "TConfig -- Error, anode and/or cathode informations missing." << endl;
		exit(0);
	}
	
	XMLElement* EbarElem = docHandle.FirstChildElement( "Detector" ).FirstChildElement( "EbarTableCalculationSteps" ).ToElement();
	if ( EbarElem )
		EbarElem->QueryIntText( &EbarTableCalculationSteps );
	else
		EbarTableCalculationSteps = 0;
	
	XMLElement* thrElem = docHandle.FirstChildElement( "Detector" ).FirstChildElement( "Threshold" ).ToElement();
	if( thrElem )
		thrElem->QueryDoubleText (&threshold);
	else {
		cerr << "TConfig -- Error, detection threshold information missing." << endl;
		exit(0);
	}
	
	/* Gas configuration */
	if(  docHandle.FirstChildElement("Gas").ToElement() == NULL){
		cerr << "TConfig -- Error reading config file\t";
		cerr << "Gas field missing in config file." << endl;
		cerr << XML_ERROR_PARSING_ELEMENT << endl;
		exit(0);
	}
	XMLElement* temp = docHandle.FirstChildElement("Gas").FirstChildElement("temperature").ToElement();
	XMLElement* pres = docHandle.FirstChildElement("Gas").FirstChildElement("pressure").ToElement();
	if (temp==NULL or pres==NULL){
		cerr << "TConfig -- Error reading config file\t";
		cerr << "Gas pressure  or  gas temperature   field(s) missing in config file." << endl;
		exit(0);
	}
	temp->QueryDoubleText( &gasTemperature );
	pres->QueryDoubleText( &gasPressure );
	
	for (XMLElement* child = docHandle.FirstChildElement("Gas").FirstChildElement("gas").ToElement(); child != NULL; child = child->NextSiblingElement("gas")) {
		string name;
		double percentage;
		
		XMLElement* elName =  child->FirstChildElement("name");
		XMLElement* elPercentage =  child->FirstChildElement("percentage");
		if (elName==NULL or elPercentage==NULL) {
			cerr << "TConfig -- Error reading config file\t";
			cerr << "Missing   gas name   or   gas percentage   information(s) in config file." << endl;
		exit(0);
		}
		name = elName->GetText();
		elPercentage->QueryDoubleText( &percentage );
		
		nGases++;
		gasNames.push_back(name);
		gasPercentage.push_back(percentage);
	}
	if (nGases < 1){
		cerr << "TConfig -- Error reading config file\t";
		cerr << "Missing gas information in config file." << endl;
		exit(0);
	}
	
	/* Simulation configuration */
	XMLElement* simElement = docHandle.FirstChildElement("Simulation").ToElement();
	if ( simElement == NULL){
		cerr << "TConfig -- Error reading config file\t";
		cerr << "Simulation field missing in config file." << endl;
		cerr << XML_ERROR_PARSING_ELEMENT << endl;
		exit(0);
	}
	
	if ( simElement->FirstChildElement( "Particle" )->FirstChildElement( "name" ) == NULL
		 or simElement->FirstChildElement( "Particle" )->FirstChildElement( "momentum" ) == NULL
		 or simElement->FirstChildElement( "Particle" )->FirstChildElement( "x0" ) == NULL
		 or simElement->FirstChildElement( "Particle" )->FirstChildElement( "theta" ) == NULL
		 or simElement->FirstChildElement( "Events" ) == NULL
		 or simElement->FirstChildElement( "Threads" ) == NULL
		 or simElement->FirstChildElement( "Generator" ) == NULL ) {
			
		cerr << "TConfig -- Error reading config file\t";
		cerr << "particle infos  or  events  or  threads  or  generator   missing in config file." << endl;
		cerr << XML_ERROR_PARSING_ELEMENT << endl;
		exit(0);
		}
	particleName = simElement->FirstChildElement( "Particle" )->FirstChildElement( "name" )->GetText();
	simElement->FirstChildElement( "Particle" )->FirstChildElement( "momentum" )->QueryDoubleText( &particleMomentum );
	simElement->FirstChildElement( "Particle" )->FirstChildElement( "x0" )->QueryDoubleText( &x0 );
	simElement->FirstChildElement( "Particle" )->FirstChildElement( "theta" )->QueryDoubleText( &theta );
	simElement->FirstChildElement( "Events" )->QueryIntText( &nEvents );
	simElement->FirstChildElement( "Threads" )->QueryIntText( &nThreads );
	generator = simElement->FirstChildElement( "Generator" )->GetText();
	
	XMLElement* outFileChild = docHandle.FirstChildElement( "Simulation" ).FirstChildElement( "OutFile" ).ToElement();
	if ( outFileChild )
		outFile += outFileChild->GetText();
	else
		outFile += "out.dat";
	
	XMLElement* seedElemnt = docHandle.FirstChildElement( "Simulation" ).FirstChildElement( "GlobalSeed" ).ToElement();
	if (seedElemnt)
		seedElemnt->QueryIntText( &globalSeed );
	else
		globalSeed = -1;
	
	XMLElement* GarSeed = docHandle.FirstChildElement( "Simulation" ).FirstChildElement( "GarfieldSeed" ).ToElement();
	if (GarSeed) {
		GarSeed->QueryIntText( &garfieldSeed );
	}
	
	if ( simElement->FirstChildElement( "Verbose" ) )
		simElement->FirstChildElement( "Verbose" )->QueryBoolText( &verbose );
		
	if ( simElement->FirstChildElement( "Snapshots" ) )
		simElement->FirstChildElement( "Snapshots" )->QueryBoolText( &snaps );
		
	if ( simElement->FirstChildElement( "VerbosityLevel" ) )
		simElement->FirstChildElement( "VerbosityLevel" )->QueryIntText( &verbosityLevel );
	
	XMLElement* SCElement = docHandle.FirstChildElement( "Simulation" ).FirstChildElement( "SingleCluster" ).ToElement();
	if ( SCElement ) {
		singleCluster = true;
		SCElement->FirstChildElement( "n0" )->QueryIntText( &n0 );
		SCElement->FirstChildElement( "x0" )->QueryDoubleText( &x0 );
	}
	else
		singleCluster = false;
		
		
	if ( simElement->FirstChildElement( "NoAvalanche" ) )
		simElement->FirstChildElement( "NoAvalanche" )->QueryBoolText( &noAvalanche );
	
	/*if ( simElement->FirstChildElement( "Efficiency" ) )
		simElement->FirstChildElement( "Efficiency" )->QueryBoolText( &computeEfficiency );*/
		
	if ( simElement->FirstChildElement( "OnlyMultiplication" ) )
		simElement->FirstChildElement( "OnlyMultiplication" )->QueryBoolText( &onlyMult );
		
	if ( simElement->FirstChildElement( "DebugOutput" ) )
		simElement->FirstChildElement( "DebugOutput" )->QueryBoolText( &debugOutput );
			
}


void TConfig::print () {
	cout << "Configuration:" << endl;
	cout << "   Detector:" << endl;
	cout << "\t gap width: " << gapWidth << endl;
	cout << "\t steps: " << nSteps << endl;
	cout << "\t electric field: " << ElectricField << endl;
	cout << "\t ebar calculation steps: " << EbarTableCalculationSteps << endl;

	cout << "\t anode: " << anodeWidth << " cm, epsilon=" << anodePermittivity << endl;
	cout << "\t cathode: " << cathodeWidth << " cm, epsilon=" << cathodePermittivity << endl;
	cout << "\t threshold: " << threshold << " pC" << endl;
	cout << "   Gas composition:" << endl;
	cout << "\t temperature: " << gasTemperature << endl;
	cout << "\t pressure: " << gasPressure << endl;
	cout << "\t gases: " << nGases << endl;
	for(int i=0; i<nGases; i++)
		cout << "\t gas " << i << " -- name: " << gasNames[i] << " - percentage: " << gasPercentage[i] << endl;
		
	cout << "   Simulation parameters:" << endl;
	cout << "\t particle: " << particleName << endl;
	cout << "\t momentum: " << particleMomentum << endl;
	cout << "\t x0: " << x0 << endl;
	cout << "\t theta: " << theta << endl;
	cout << "\t events to simulate: " << nEvents << endl;
	cout << "\t threads to use: " << nThreads << endl;
	cout << "\t outfile path: " << outFile << endl;
	cout << "\t generator: " << generator << endl;
	if (garfieldSeed != -1)
		cout << "\t garfield Seed: " << garfieldSeed << endl;
	if (globalSeed == -1)
		cout << "\t random global seed" << endl;
	else
		cout << "\t global seed: " << globalSeed << endl;
	if (singleCluster) {
		cout << "\t single cluster at step ";
		if (x0 >= 0)
			cout << x0;
		else
			cout << "random";
		cout << " containing ";
		if (n0 >=0)
			cout << n0; 
		 else
			cout << "random";
		cout << " electrons" << endl; 
	}
	if (noAvalanche) {
		cout << "\t No avalanche simulation" << endl; 
	}
	if (onlyMult) {
		cout << "\t Only multiplication is simulated. Every other processes are deactivated." << endl; 
	}
	/*if (computeEfficiency) {
		cout << "\t Efficiency runs" << endl; 
	}*/
	cout << "\t Debug outputs: " << debugOutput  << endl; 
	cout << "\t verbose: " << verbose << endl;
	cout << "\t verbosity level: " << verbosityLevel << endl;
	cout << "\t snapshots: " << snaps << endl;
	
	
}
