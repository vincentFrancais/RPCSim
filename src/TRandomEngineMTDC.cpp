#include "TRandomEngineMTDC.hpp"

TRandomEngineMTDC::TRandomEngineMTDC(uint16_t id, int stateSeed, int globalSeed) : TRandomEngine() {
	fMTStruct = get_mt_parameter_id_st(32,521,id,stateSeed);
	sgenrand_mt(globalSeed, fMTStruct);
}

TRandomEngineMTDC::~TRandomEngineMTDC() {
	free_mt_struct(fMTStruct);
}

double TRandomEngineMTDC::RandU01() {
	return genrand_mt(fMTStruct)*(1.0/4294967296.0);
}

