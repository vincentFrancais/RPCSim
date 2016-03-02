#include "TConfig.hpp"

using namespace std;
using namespace tinyxml2;

TConfig readConfigFile(string fileName){
	
	TConfig config;
	
	XMLDocument doc;
	if ( doc.LoadFile( fileName.c_str() ) != XML_SUCCESS ){
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
	
	detectorElement->FirstChildElement( "GapWidth" )->QueryDoubleText( &config.gapWidth );
	detectorElement->FirstChildElement( "Steps" )->QueryIntText( &config.nSteps );
	detectorElement->FirstChildElement( "ElectricField" )->QueryDoubleText( &config.ElectricField );
	
	for (XMLElement* child = docHandle.FirstChildElement("Detector").FirstChildElement("ResistiveLayers").FirstChildElement("layer").ToElement(); child != NULL; child = child->NextSiblingElement("layer")){
		double width, epsilon;
		child->FirstChildElement( "width" )->QueryDoubleText( &width );
		child->FirstChildElement( "resistivity" )->QueryDoubleText( &epsilon );
		
		config.nResisLayers++;
		config.resisLayersWidth.push_back(width);
		config.resisLayersEpsilon.push_back(epsilon);
	}
	
	docHandle.FirstChildElement("Gas").FirstChildElement("temperature").ToElement()->QueryDoubleText( &config.gasTemperature );
	docHandle.FirstChildElement("Gas").FirstChildElement("pressure").ToElement()->QueryDoubleText( &config.gasPressure );
	
	for (XMLElement* child = docHandle.FirstChildElement("Gas").FirstChildElement("gas").ToElement(); child != NULL; child = child->NextSiblingElement("gas")){
		string name;
		double percentage;
		 
		name = child->FirstChildElement("name")->GetText();
		child->FirstChildElement( "percentage" )->QueryDoubleText( &percentage );
		
		config.nGases++;
		config.gasNames.push_back(name);
		config.gasPercentage.push_back(percentage);
	}
	
	XMLElement* simElement = docHandle.FirstChildElement("Simulation").ToElement();
	if ( simElement == NULL){
		cerr << "TConfig -- Error reading config file" << endl;
		cerr << XML_ERROR_PARSING_ELEMENT << endl;
		exit(0);
	}
	
	config.particleName = simElement->FirstChildElement( "Particle" )->FirstChildElement( "name" )->GetText();
	simElement->FirstChildElement( "Particle" )->FirstChildElement( "momentum" )->QueryDoubleText( &config.particleMomentum );
	simElement->FirstChildElement( "Particle" )->FirstChildElement( "x0" )->QueryDoubleText( &config.x0 );
	simElement->FirstChildElement( "Particle" )->FirstChildElement( "theta" )->QueryDoubleText( &config.theta );
	simElement->FirstChildElement( "Events" )->QueryIntText( &config.nEvents );
	simElement->FirstChildElement( "Threads" )->QueryIntText( &config.nThreads );
	
	XMLElement* outFileChild = docHandle.FirstChildElement( "Simulation" ).FirstChildElement( "outFile" ).ToElement();
	if ( outFileChild )
		config.outFile += outFileChild->GetText();
	else
		config.outFile += "out.dat";
	//config.outFile = simElement->FirstChildElement( "outFile" )->GetText();
	
	return config;
}

void printConfig(TConfig config){
	cout << "Configuration:" << endl;
	cout << "   Detector:" << endl;
	cout << "\t gap width: " << config.gapWidth << endl;
	cout << "\t steps: " << config.nSteps << endl;
	cout << "\t electric field: " << config.ElectricField << endl;
	cout << "\t resistive layers: " << config.nResisLayers << endl;
	for(int i=0; i<config.nResisLayers; i++)
		cout << "\t layer " << i << " -- width: " << config.resisLayersWidth[i] << " - epsilon: " << config.resisLayersEpsilon[i] << endl;
		
	cout << "   Gas composition:" << endl;
	cout << "\t temperature: " << config.gasTemperature << endl;
	cout << "\t pressure: " << config.gasPressure << endl;
	cout << "\t gases: " << config.nGases << endl;
	for(int i=0; i<config.nGases; i++)
		cout << "\t gas " << i << " -- name: " << config.gasNames[i] << " - percentage: " << config.gasPercentage[i] << endl;
		
	cout << "   Simulation parameters:" << endl;
	cout << "\t particle: " << config.particleName << endl;
	cout << "\t momentum: " << config.particleMomentum << endl;
	cout << "\t x0: " << config.x0 << endl;
	cout << "\t theta: " << config.theta << endl;
	cout << "\t events to simulate: " << config.nEvents << endl;
	cout << "\t threads to use: " << config.nThreads << endl;
	cout << "\t outfile path: " << config.outFile << endl;
}
