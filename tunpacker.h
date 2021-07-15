#ifndef TUNPACKER
#define TUNPACKER

/////////////////////////////////////////////////////////////////////
//   Unpacker
/////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include <TFile.h>  
#include <TObject.h> 
#include <string>
#include "thldevent.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

using namespace std;

enum EDsState {kDsOk=0,kDsEndFile=1,kDsEndData=2,kDsError=3};

class Unpacker: public TObject
{

public:

 Unpacker(); 
 Unpacker(Int_t j,const char* name,Int_t nEvt=50000, const char* subEvtId="899|871",Int_t referenceChannel=95, const char* fpga_code = "");
 Unpacker(const char* dir,const char* name,const char* odir,Int_t nEvt=50000,
	  const char* subEvtId="",Int_t referenceChannel=95, const char* fpga_code = "",
	  Int_t min = -100000, Int_t max = -100000, Int_t quietMode = 1, Int_t fullSetup = 0, Int_t VHR = 0);
 Unpacker(const char* dir, const char* name, const char* odir,const char* ofile, Int_t nEvt);
 Unpacker(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n);
 Unpacker(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n,Int_t n2);
 Unpacker(const char* dir,const char* name, const char* odir, Int_t nEvt, TString luptab, TString calpar);
 void unpackerFast(const char* dir,const char* name, const char* odir, Int_t nEvt, TString luptab, TString calpar);


 ~Unpacker();

 void fillCalibration(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n);
 void fillHistograms(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n,Int_t n2);
 Int_t syncCheck(const char* dir, TString file, Int_t nEvt,Int_t n);


 UInt_t HexStrToInt(const char* str) {
	UInt_t t;
	std::stringstream s;
	s << std::hex << str;
	s >> t;
	return t;
 }

 // THIS PART SHOULD BE IN A DIFFERENT PLACE IS JUST A FAST SOLUTION
 // UNTIL WE HAVE A PROPERLY DEFINED RUNTIMEDATABASE!
 // initialized only once!

 TString fileLookupPar;         //! * name of the file with lookup params
 TString fileHitFinderPar;      //! * name of the file with hit finder params
 TString fileTrackFinderPar;    //! * name of the file with track finder params
 TString fileHitFinderParOut;      //! * name of the file with hit finder params

 void setFileLookupPar( TString fileName) {         fileLookupPar      = fileName;}
 void setFileHitFinderPar( TString fileName) {      fileHitFinderPar   = fileName;}
 void setFileHitFinderParOut( TString fileName) {      fileHitFinderParOut   = fileName;}
 void setFileTrackFinderPar( TString fileName) {    fileTrackFinderPar = fileName;}

 TString getFileLookupPar( void ) {        return fileLookupPar;}
 TString getFileHitFinderPar( void ) {     return fileHitFinderPar;}
 TString getFileHitFinderParOut( void ) {     return fileHitFinderParOut;}
 TString getFileTrackFinderPar( void ) {   return fileTrackFinderPar;}

 Bool_t setRootFile(const char* filename="test.root"); //to set the root output file
 string setInputFile(const char* filename); //to set the hld input file
 string setInputFile(const char* dir,const char* filename); //to set the hld input file
 Bool_t eventLoop(Int_t NbEvt=50000,Int_t startEvt=0);
 Bool_t setpEvent(Int_t subId);	 //it sets pEvent by reading hld file
 HldEvent* getpEvent(void) {return pEvent;}
 Int_t getEventNr() { return EventNr; }
 Int_t getEventLimit() { return EventLimit; }
// Bool_t eventLoopFillCal(Int_t nbEvt, Int_t startEv, TH1D** hq, TH1D** hdt);
 Bool_t eventLoopFillCal(Int_t nbEvt, Int_t startEv, TH1D** hq, TH1D** hdt, TH1D** hdt2);
 Bool_t eventLoopFillCal(Int_t nbEvt, Int_t startEv, TH1D** h1D, TH2D** h2D, TH3D** h3D);
 Int_t eventLoopSyncCheck(Int_t nbEvt,Int_t startEv);
protected:
    void setpEvent(HldEvent* evt) { pEvent=evt; }
    HldEvent* pEvent; //Current event read from file
    Int_t EventNr;   //Event Counter
    Int_t EventLimit; //Maximum event number per file
    Int_t subEvtId;
    TFile* pRootFile; // pointer to TFile with the output tree
    string inputFile; //wk 28.05
    string outputFile; //wk 28.05
    Int_t fpga_code; // address of the data source (e.g. given fpga ) decoded from hld file
    Int_t refCh;

    /*
    TH1D *hq[3*10*12];
    TH1D *hdta[3*10*12];
    TH1D* hdt[10*12*10*12+1];
    TH1D* h1D[10];
    TH2D* h2D[8];
    TH3D* h3D[6];
    */

public:
 ClassDef(Unpacker,1);
};

#endif /* !TUNPACKER */

