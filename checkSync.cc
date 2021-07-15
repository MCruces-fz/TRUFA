#include "tunpacker.h"
#include "TString.h"
using namespace std;

Int_t doSyncCheck(char* pathToData, char* name, char* lookupTable) {
    TString fName(name);
    Unpacker unpack =  Unpacker();
    //unpack.setFileLookupPar("../pars/luptab.txt");
    unpack.setFileLookupPar(lookupTable);
    return unpack.syncCheck(pathToData,fName,1000000,0);
}

int main(int argc, char* argv[] ) {
    int isSync = doSyncCheck(argv[1],argv[2],argv[3]);
    return isSync;
}
