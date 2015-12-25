#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>
#include <vector>
#include <fstream>

#include "TResult.hpp"

using namespace std;

int main(){
	
	int inFD;
	TResult result;
	
	ofstream outFile("thrCross.dat", ios::out | ios::trunc);
	
	inFD = open("../out/out.dat", O_RDONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	if (inFD == -1){
        cerr << "Failed to open file" << endl;
        return 0;
    }
    
    while (read(inFD, &result, sizeof(result)) > 0){
		cout << "Dt: " << result.Dt << endl;
		cout << "Dx: " << result.Dx << endl;
		cout << "nSteps: " << result.iNstep << endl;
		cout << "test: " << result.thrCrossTimeStep << endl;
		cout << "status: " << result.avalStatus << endl;
		cout << "charges_size: " << result.charges_size << endl;
		cout << "chargesTot_size: " << result.chargesTot_size << endl;
		cout << "signal_size: " << result.signal_size << endl;
		cout << "========================" << endl;
		
		outFile << result.thrCrossTimeStep << endl;
    }
    close(inFD);
    outFile.close();
    
    return 0;
}
