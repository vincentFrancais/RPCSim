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
	verbose = false;
	verbosityLevel = 0;
	snaps = false;
	
	XMLDocument doc;
	if ( doc.LoadFile( file.c_str() ) != XML_SUCCESS ){
		cerr << "TConfig -- Error opening file" << endl;
		exit(0);
	}
	
	XMLHandle docHandle( &doc );
	
	XMLElement* detectorElement = docHandle.FirstChildElement("Detector").ToElement();
	if( detectorElement == NULL){
		cerr << "TConfig -- Error reading config file" << endl;
		cerr << XML_ERROR_PARSING_ELEMENT << endl;
		exit(0);
	}
	
	detectorElement->FirstChildElement( "GapWidth" )->QueryDoubleText( &gapWidth );
	detectorElement->FirstChildElement( "Steps" )->QueryIntText( &nSteps );
	detectorElement->FirstChildElement( "ElectricField" )->QueryDoubleText( &ElectricField );
	
	for (XMLElement* child = docHandle.FirstChildElement("Detector").FirstChildElement("ResistiveLayers").FirstChildElement("layer").ToElement(); child != NULL; child = child->NextSiblingElement("layer")){
		double width, epsilon;
		child->FirstChildElement( "width" )->QueryDoubleText( &width );
		child->FirstChildElement( "resistivity" )->QueryDoubleText( &epsilon );
		
		nResisLayers++;
		resisLayersWidth.push_back(width);
		resisLayersEpsilon.push_back(epsilon);
	}
	
	docHandle.FirstChildElement("Gas").FirstChildElement("temperature").ToElement()->QueryDoubleText( &gasTemperature );
	docHandle.FirstChildElement("Gas").FirstChildElement("pressure").ToElement()->QueryDoubleText( &gasPressure );
	
	for (XMLElement* child = docHandle.FirstChildElement("Gas").FirstChildElement("gas").ToElement(); child != NULL; child = child->NextSiblingElement("gas")){
		string name;
		double percentage;
		 
		name = child->FirstChildElement("name")->GetText();
		child->FirstChildElement( "percentage" )->QueryDoubleText( &percentage );
		
		nGases++;
		gasNames.push_back(name);
		gasPercentage.push_back(percentage);
	}
	
	XMLElement* simElement = docHandle.FirstChildElement("Simulation").ToElement();
	if ( simElement == NULL){
		cerr << "TConfig -- Error reading config file" << endl;
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
	//outFile = simElement->FirstChildElement( "outFile" )->GetText();
	
	XMLElement* seedsElemnt = docHandle.FirstChildElement( "Simulation" ).FirstChildElement( "Seeds" ).FirstChildElement("seed").ToElement();
	while (seedsElemnt) {
		int s;
		seedsElemnt->QueryIntText( &s );
		seeds.push_back(s);
		seedsElemnt = seedsElemnt->NextSiblingElement();
	}
	
	XMLElement* GarSeed = docHandle.FirstChildElement( "Simulation" ).FirstChildElement( "GarfieldSeed" ).ToElement();
	if (GarSeed) {
		GarSeed->QueryIntText( &garfieldSeed );
	}
	
	if ( simElement->FirstChildElement( "Verbose" )->ToElement() )
		simElement->FirstChildElement( "Verbose" )->QueryBoolText( &verbose );
		
	if ( simElement->FirstChildElement( "Snapshots" )->ToElement() )
		simElement->FirstChildElement( "Snapshots" )->QueryBoolText( &snaps );
		
	if ( simElement->FirstChildElement( "VerbosityLevel" )->ToElement() )
		simElement->FirstChildElement( "VerbosityLevel" )->QueryIntText( &verbosityLevel );
	
	XMLElement* SCElement = docHandle.FirstChildElement( "Simulation" ).FirstChildElement( "SingleCluster" ).ToElement();
	if ( SCElement ) {
		singleCluster = true;
		SCElement->FirstChildElement( "n0" )->QueryIntText( &n0 );
		SCElement->FirstChildElement( "x0" )->QueryDoubleText( &x0 );
		if (n0 <= 0) {
			cerr << "TConfig -- n0 must be greater than 0. Aborting" << endl;
			exit(0);
		}
		else if (x0 < 0) {
			cerr << "TConfig -- n0 must be equal to or greater than 0. Aborting" << endl;
			exit(0);
		}
	}
	else
		singleCluster = false;
	
	//return config;
}

void TConfig::print () {
	cout << "Configuration:" << endl;
	cout << "   Detector:" << endl;
	cout << "\t gap width: " << gapWidth << endl;
	cout << "\t steps: " << nSteps << endl;
	cout << "\t electric field: " << ElectricField << endl;
	cout << "\t resistive layers: " << nResisLayers << endl;
	for(int i=0; i<nResisLayers; i++)
		cout << "\t layer " << i << " -- width: " << resisLayersWidth[i] << " - epsilon: " << resisLayersEpsilon[i] << endl;
		
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
	if (seeds.size() > 0) {
		cout << "\t seeds: ";
		for (uint i=0; i<seeds.size(); i++)
			cout << seeds.at(i) << "  ";
		cout << endl;  
	}
	if (singleCluster) {
		cout << "\t single cluster at " << x0 << " containing " << n0 << " electrons" << endl; 
	}
	cout << "\t verbose: " << verbose << endl;
	cout << "\t verbosity level: " << verbosityLevel << endl;
	cout << "\t snapshots: " << snaps << endl;
	
	
}
