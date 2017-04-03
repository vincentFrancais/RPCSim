#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <utility>
#include <fcntl.h>
#include <cassert>

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"

#include "TDetector.hpp"
#include "TAvalanche1D.hpp"
#include "helper_functions.hpp"
#include "TThreadsFactory.h"
#include "TResult.hpp"
#include "TConfig.hpp"

#include "argparse.hpp"

#define SFMT_MEXP 19937

#include "SFMT/SFMT.h"
#include "SFMT/SFMT-jump.h"
#include "MT/dc.h"
#include "TRandomEngineMTDC.hpp"
#include "TRandomEngineSFMT.hpp"
#include "TRandomEngineMRG.hpp"
#include "TRandomEngineMT.hpp"

static const char * jump10_20 = "52641a85b9d278c3a590f230feb4716171f0560e4afb7e51a4582af14cdec440f60e727fa787535c0a2de7c0bd55b1b6b0d2af28fc9c6bb8a734707e90e9e2b5302c5c91fc83ac019c79badacc9e240bb826a685075106998b7734ce1b2e89d84894a930fb1ba5f03e4adf2c91bd68cf8e7c7cb40544e316eb9f504df1fe92b847b454834e44414a2ff5f1d4fc64a198d0d772fd9aba39e8bd26e5341e955c00c01f3499e6d2f7df34ccbc0e4e65a3c2f827d771862295890e6fffed9208bcbe1c7ac6fa657b7baf47192cd31c788d650ad54ac9cdea4ed6d73deeea2940e0cea10b3336fb1f2a884a827d36f296cc8a0c6131b509f13adcfd448999ccd461b476dae87ab411e7dac424b2fd44f100a7d72bd16d7f8e0574477c7687f073ffc8a039e53db66457e31e378da7787f436a55e31ba4d3aa7c46098ea3131fd0784ff7c94f1d52561024f662052f2a13d8a862c619baa89d3d5c2557c066e98ae77f152a509ff3496c14a2467bb962dd3550c87a2fe43b602292e3d534c326c3a7736400056eee5b074ed0f985d782139e3bb3ab8af26edd18300ea14da59c594c07f27bc9cc51ad777c8a659faece654693a205a9122f7c5bf42ad3f430e1599b19bd74afd43b0ad769938384250ac20176dffc80d08e9b20fb10c5ab74da5e84558164e4d36f64762a8f50c5bb10ca7f765c5bacfed61516b68b8162c6d4d9c665db748ca86e7046aa3763b3dad114049a3494b634eb902cf2e790c52f19f0cec37d8769ebbe38c341d238395687f4fe6ff1995051030575430317125a99a2ec68a82ed8ff6a22c4a196f52bbdad18ed372fa924d1e45fb033d70898d4aba98d58bca4b33dce440d857bad063b3f71e928af460c9787cb2dfdd3ebbc856d945a69256724b870c388e04f48cf277f97eae697c60152219a7a17ce6093c759a8711a4b1178528b99bbf0b99b1fdda49f513f0fe08544cb5bb4099ce4ea40c61ae27d4cacf5d6742662d08775653d4ce25c6463ad86e13ef7aa669b721f6c4481e45e630274a41e7844835c120a033314b40cb2eaa041bf6671531e569f6747ed7b70b5f68d3e189bb0a97fc7d9fccbbc02596f5201a7d63fe273ff7b0e4a0cfd52120abba0b0d264239bb873f35d37efe33343e8b1b2b36ed52d510d0859aeae76bcc31708dbcd5bcc484287ec31bd91da83e9cf1541e4224b0e26caee9a6a87439176342be07a03424fa83ae23d032c0b2ee704e5c7128cb089c1a82d80f43669f26758a58fae0af67cc19e0b23ecec463b1de4bd578fc0c28a88d2af607cdf81ae2b3010c92613f7a8f20ca566bdff57f318826183450904ac39575b1b88443e89a653cad578e68b69add392cda9011eefcf12e1e08bab39e424c9c831511ccadc2ec60931e55100476122131c58ecd088fcc8f0262884336f2320f2415f201716ae606e38742123c0d6534662024dd12ddbc7884be07c558916387c39b6275ae77d0c6c6fa0dadbcdd8d54827051a682534391b462ce6d1492af312e69b053cc817b9f0c560491f99de20748396cb459a27a56d3406bc2977f14c2b4ff66af6f73eb4aad635b1622e9e5f55f56e09dbbb39ad8b86b5110a6674b4d77b40c83d097a7802592d94e7d32cf6453b9ffb0d47f7a20cc0c94b17a4f318276e8330c43a7f4aae3e5281e924084ae06879d3597983e355bdb7a69410cd81f58921d095b6c00edcbbad9463e8e6f72a0331186ec471ba9a3f2ea3e4474e1e8260dae925afa18fc2cb0c907db96994d49118ef4e0a5d172fe2a52e06b6dea72a201428893f3baf787111ee6fcc68c9d20dabaf72bcb100a7abd20f70aa9d6a0c15084e4f8388cc9a98e154d7067ea21b0305ccd8511dd93474ccb6e5f3f2705af17cf51631fae1c587110a277c9807cadd7c468f96b76b32e3a3f34ad903fbf7a9a027593a42be182e32e6acfe71d2b9e248b2c7509c630cb3fa84aa179e38fa4a8d5927b87b25f3ef8d531ff2abca0940ff3491d82cff185aa79e877a9d833e2a821d45fdd3cca0a558509a906566be0cd455f8148fc8de21870299651c25bc01963c3286f8a9cfffe52dd38b41aa537da5677cdc93af9cfa0e4ec79ead2a0c409bf408830f6eda95db894a9f9b6fe4a360549727cd0ce2211f901819a2bef3863fb8c0e51c83bfc5a7e4a030c053f0f49afba66a8650467f39104b83b5ac2fa21f118de070eb2144ec3460ed45797e998c297015665de52de88dc6d1445488b634211ad8191e1ac20f956b0aa2dfa04a6cad82f936b78f86788362bfb68e791c1961bd167507cec7ad8fcff47897f5ace5caa09881e3bd39603c246c79400938d27b014e045e76e92e54c7b27e4f5a08e200fd3132367e727ce1f98df19f3569e07b9101d14a32f6e189f02b2c535a58eaa5686c570ee4a2560eec19855584273f600d839d045d1e00e1a7e500028bb806da420af6fd5da805267850e624f8c1e94dcef34aa2fc50c839156f80d582018670426ad091f0336bebf4c11affa38d9801e067073be34fb856198ff926162a72919c7bc8012bfaa6023e0c5784bdeecca9e1faa26913168b3624fcdc6903540e0fa1b485dbb559ca68d5c91548043cb42759894a8429e9b1ec5fa4863b21c3187ca71466198977036d0bd6d9a0430e8e49ad6590b3336bb750e5c12505304b17591be6eefbbf441585f6e4b6d98e491f69477afe094c37f49b0bb8b5d1cd1c43c0c79e6c72b2cb8d6439b4e91f22c5075928b7b72a025f81821f761056441e47f06c5b52dc26ad8b52352052aba6060faeb76962ff58ef6baf6a581fec7d8f3a0eff0e8bbc789284e86b4b1742db86822838552dcd530066ad384c33d0f435a45964b93e2873989edf4b5f54fb0dc8bb5688880fe68293887b5124b6da421e213697db122e020abeabd328a4e36b5d84bb71c6167c2cdb907912489ea21a5bd249cfdd623e4b7e1ff7d0cae118960986529c1cbb27866a19a30149cfdda9c5fb506407279ca5439353e96f20ed5696a267d8f9fc8459e49f0333c8a255d70d6a57221c3c203735dd4bb34a28f395338883f281d5d38ee7887fbdfc85186e3c57f153ab66ee2afcfa7a9368c13d2ccab127f00abababfc916394dc0c8a858fabf7ce7129c29f3a00013dc083f6a96800a8cb8a132c27410fad6b656c900276ffe50e49d60373f526907d2d60ebc4dde98e8cb7e19c9c8f6084e2da8497fc20b0a31323345fae7f4d91c19c777690dc27fa2ee7166609de5e4b151f2691c6ac8f8e165961e99ec02ad050ed240512c69d974d61451409082f83876037691ea57f8e0ee2b684d44b271e99da9774d4feedb014853ae45e1e384dcfeca72acaa9abfda85452ea3a0f0388405286116ee80988e43d80bb786dc54d56e8c49c66687f20e7fa905cd1a034d86a809d5e8b0d68fab16177d46a685bac8a44e1b053b8deaa6b13c7b3a8159a8b53519f91fc2e617f09dee476448f5c3e5e4e158c7c386f343e3eced33ac92d96de4c3b147ccacb4e9eb16";

using namespace std;
using namespace Garfield;

typedef unsigned int uint;
typedef unsigned long ulong;

double cm = 0.01;
int TAvalanche::count = 0;
int TAvalanche::countSim = 0;

pthread_mutex_t gPipeLock;
pthread_mutex_t gTrackLock;
int gPipe[2];
TResult gNullResult;


/* Simple structure for passing simulation parameters to the wrapping function */
struct ThreadData{
	TDetector * detector;
	TConfig config;
	sfmt_t sfmt;
	int id;
	
	ThreadData (TDetector * det, TConfig conf, sfmt_t status, int i) : detector(det), config(conf), sfmt(status), id(i)	{ };
};

void * wrapperFunction(void * Arg){
	
	assert(Arg != NULL);

	ThreadData* data = reinterpret_cast< ThreadData* > (Arg);
	
	TResult result;
	TAvalanche1D avalanche(data->detector, data->config, data->sfmt, data->id);
	
	sem_post(TThreadsFactory::GetInstance()->GetInitLock());
	
	pthread_mutex_lock(&gTrackLock);
	avalanche.initialiseTrackHeed();
	pthread_mutex_unlock(&gTrackLock);

	//avalanche.disableSpaceChargeEffect();
	avalanche.simulateEvent();

	result = avalanche.getResultFile();

	pthread_mutex_lock(&gPipeLock);
    write(gPipe[1], &result, sizeof(result));
    pthread_mutex_unlock(&gPipeLock);

	return 0;
}

void * WriteResults(void * Arg) {	
    TResult result;
    char * outputFile = reinterpret_cast<char *>(Arg);

    int outFD;
	
    /* Open the output file */
    string outFileBinary (outputFile);
    size_t found = outFileBinary.find(".dat");
	if (found!=string::npos)
		outFileBinary.erase(found,4);
    outFileBinary += ".bin";

    outFD = open(outFileBinary.c_str(), O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    assert(outFD != -1);

    /* Read the incoming event */
    while (read(gPipe[0], &result, sizeof(TResult)) > 0)
    {
        /* End signal => stop dealing with data */
        if (memcmp(&result, &gNullResult, sizeof(TResult)) == 0)
            break;

        /* Write the event to the output file */
        write(outFD, &result, sizeof(TResult));
    }

    close(outFD);

    return 0;
}

int main(int argc, const char** argv) {
	/* Parsing argument ... Not very useful anymore. Could be removed. */
	ArgumentParser parser;
	parser.addArgument("-c", "--config", 1, false);
	parser.addArgument("-g", "--gap",1);
	parser.addArgument("-e", "--eff", 1);
	parser.parse(argc, argv);

	/* Read config file */
	TConfig config;
	config = TConfig ( parser.retrieve<string>("config") );
	
	if (parser.count("g"))
		config.gapWidth = stod (parser.retrieve<string> ("g"));
	if ( parser.count("e") )
		config.outFile = "out/eff-"+parser.retrieve<string>("e");
		
    config.print();
    
    char outputFile[PATH_MAX];
    
    pthread_t writingThread;
    void * ret;
	
	strncpy(outputFile, config.outFile.c_str(), PATH_MAX - 1);
	outputFile[PATH_MAX - 1] = '\0';
	
    unsigned int nThreads = config.nThreads;
    unsigned long nEvents = config.nEvents;
	
	/* Init the SFMT status */
	sfmt_t SFMT;
	if ( config.globalSeed != -1 ) {		
		cout << "Init SFMT generator with seed: " << config.globalSeed << endl;
		sfmt_init_gen_rand(&SFMT, config.globalSeed);
	}
	else {
		uint seed = getUUID();
		cout << "Init SFMT generator with random seed: " << seed << endl;
		sfmt_init_gen_rand(&SFMT, seed);
		config.globalSeed = seed;
	}
    
    /* Initialize our pipe lock */
    pthread_mutex_init(&gPipeLock, 0);
    pthread_mutex_init(&gTrackLock, 0);


    /* Start our threads factory */
    TThreadsFactory::GetInstance()->SetMaxThreads(nThreads);
    
    
    /* Init our null event */
    memset(&gNullResult, 0, sizeof(TResult));
    
    
    /* Init our detector */
	TDetector* detector = new TDetector(config);
	if ( parser.count("e") ) {
		detector->setElectricField( stod(parser.retrieve<string>("e"))*1e3, 0, 0 );
		cout << "====== Efficiency simulation run, HV at " << parser.retrieve<string>("e") << " ======" << endl;
	}
	detector->initialiseDetector();
	detector->writeGasTransportParameters();
	
	/* Functions to produce data on primary inisation */
	//TAvalanche::computeClusterDensity(detector,"muon",6e7,1.5e10,600);
	//TAvalanche::computeElectronsProduction(detector,"muon",5.e9,60000);
	
	/* Here we define a Magboltz Gas in order to print its photo-absorption CS through HEED */
	//MediumMagboltz* gas = new MediumMagboltz();
	//gas->SetComposition("Ar", 100.);
	//gas->SetTemperature(293.15);
	//gas->SetPressure(760);
	//detector->setGasMixture(gas);
	//TDetector::printPACSData(gas);
	//detector->printPACSData();
	//delete gas;
	
	if (config.noAvalanche){
		delete detector;
		return 0;
	}

	
	
	/* Init struct of simulation parameters */
	ThreadData* data = new ThreadData(detector, config, SFMT, 0);
	
	
    /* Open the communication pipe */
    if (pipe(gPipe) == -1){
		goto end;
    }


    /* Start our background writing thread */
    if (pthread_create(&writingThread, 0, WriteResults, outputFile) != 0){
		goto end2;
    }


    /* Hot loop, the simulation happens here */
    for (unsigned long i = 0; i < nEvents; i++){
		/* the SFMT status is given to the thread and then jump-ahead by 10^20 numbers (ensure independant and large enough streams) */
		data->sfmt = SFMT;
		data->id = i+2; // Increment the indice by 2 because we use 2 instance of DCMT in TAvalanche
		TThreadsFactory::GetInstance()->CreateThread(wrapperFunction, data);
		SFMT_jump(&SFMT, jump10_20);
    }


    /* Wait for all the propagations to finish */
    TThreadsFactory::GetInstance()->WaitForAllThreads();

    /* Send the end signal to writer */
    write(gPipe[1], &gNullResult, sizeof(TResult));

    /* Wait for the end of the writer */
    pthread_join(writingThread, &ret);


end2:
    close(gPipe[0]);
    close(gPipe[1]);
end:
    pthread_mutex_destroy(&gPipeLock);
    pthread_mutex_destroy(&gTrackLock);
    delete detector;
    delete data;

    return 0;
}
