#include "trpccalpar.h"
#include <fstream>

using namespace std;

ClassImp(TRpcCalPar);

TRpcCalPar::TRpcCalPar( ){
    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		for(Int_t l=0;l<2;l++) {
		    pars[i][j][k][l] = 0.;
		}
	    }
	}
    }
}
TRpcCalPar::TRpcCalPar(TString filename){
    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		for(Int_t l=0;l<2;l++) {
		    pars[i][j][k][l] = 0.;
		}
	    }
	}
    }
    readParams(filename);
}
Bool_t TRpcCalPar::readParams(TString filename) {

    std::ifstream file(filename.Data());
    Int_t detid;
    Int_t col, row;
    Float_t qt, toff;
    Bool_t reading = kTRUE;
    Int_t num;
    while(!file.eof()) {
	file>>detid>>row>>col>>qt>>toff;
	//printf("%i %i %i %f %f\n",detid,col,row,qt,toff);
        pars[detid][row][col][0] = qt;
        pars[detid][row][col][1] = toff;
    }
    file.close();
    return kTRUE;
}
Bool_t TRpcCalPar::writeParams(TString filename) {
    return kTRUE;
}
