#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <utility>
#include <sys/stat.h>
#include <typeinfo>
#include <fcntl.h>

#include "TResult.hpp"

using namespace std;

int main(){
	
	int inFD;
	TResult result;
	
	inFD = open("out.dat", O_RDONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	if (inFD == -1){
        cerr << "Failed to open file" << endl;
        return 0;
    }
    
    while (read(inFD, &result, sizeof(result)) > 0){
		//cout << "Induced charge: " << result.fInducedSignal[0] << endl;
		cout << "Induced charge: " << result.fDiffCoeff[1] << endl;
		cout << "Induced signal: " << result.fNElectrons[1] << endl;
    }
    close(inFD);
	return 0;
}
