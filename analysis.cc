#include "tunpacker.h"
#include <iostream>
#include "TString.h"
using namespace std;


void doStaffCalibration(char* name) {
    TString fName(name);
    Unpacker unpack =  Unpacker();
    TString day(fName(11,fName.First(".")-11));
    cout<<"DAY "<<day<<endl;
    cout<<"PATH to files "<<"../"+fName<<endl;
    unpack.setFileLookupPar("../pars/luptab.txt");
    //unpack.setFileHitFinderPar("../pars/time_calPar.txt");
    unpack.setFileHitFinderPar("../pars/time_calPar3planes.txt");
    unpack.setFileHitFinderParOut("../pars2/"+day+"_CalPars.txt");
    unpack.fillCalibration("/media/Datos2TB/tragaldabas/data/done/","../"+fName,"../qcalhistos2/",day+"_qhistos.root",1000000,0);
}

void doStaffAnalysis(char* name){
    TString fName(name);
    Unpacker unpack =  Unpacker();
    TString day(fName(11,fName.First(".")-11));
    cout<<"DAY "<<day<<endl;
    cout<<"PATH to files "<<"../"+fName<<endl;
    unpack.setFileLookupPar("../pars/luptab.txt");
    unpack.setFileHitFinderPar("../pars2/"+day+"_CalPars.txt");
    unpack.fillHistograms("/media/Datos2TB/tragaldabas/data/done/","../"+fName,"../histos2/",day+"_result_histos.root",1000000,0,0);
}


int main(int argc, char* argv[] ) {
    cout<<"Starting list "<<argv[1]<<endl;

    doStaffCalibration(argv[1]);
    //doStaffAnalysis(argv[1]);
    return 1;

}
