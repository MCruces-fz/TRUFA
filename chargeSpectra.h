#ifndef TFiller
#define TFiller

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
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

using namespace std;


class Filler: public TObject
{

public:

    Filler();
    ~Filler();
    void fillHistograms(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n);
    UInt_t HexStrToInt(const char* str) {
	UInt_t t;
	std::stringstream s;
	s << std::hex << str;
	s >> t;
	return t;
    }

    // initialized only once!
    TString fileLookupPar;            //! * name of the file with lookup params
    TString fileHitFinderPar;         //! * name of the file with hit finder params
    TString fileTrackFinderPar;       //! * name of the file with track finder params
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
    Bool_t setpEvent(Int_t subId);	 //it sets pEvent by reading hld file
    HldEvent* getpEvent(void) {return pEvent;}
    Int_t getEventNr() { return EventNr; }
    Int_t getEventLimit() { return EventLimit; }
    Bool_t eventLoopFill(Int_t nbEvt, Int_t startEv);

protected:
    void setpEvent(HldEvent* evt) { pEvent=evt; }
    HldEvent* pEvent;  //Current event read from file
    Int_t EventNr;     //Event Counter
    Int_t EventLimit;  //Maximum event number per file
    Int_t subEvtId;
    TFile* pRootFile;  // pointer to TFile with the output tree
    string inputFile;  //
    string outputFile; //
    Int_t fpga_code;   // address of the data source (e.g. given fpga ) decoded from hld file
    Int_t refCh;

    Double_t fTheta_intervals[7]; /*= {0.,
                                   1.29209663815835789e+01,
                                   1.84349488229220171e+01,
				   2.27864979995971453e+01,
                                   2.65650511770779936e+01,
                                   3.00000000000000036e+01,
				   60.};*/

   ////////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////
   // Define your histograms here!
   ////////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////
   TH1D* h_daq_active;
   TH1D* h_daq_sync;
   TH1D* h_daq_seconds;
   TH1D* h_daq_total;

   TH2F* h_q[3][10][12];





};

#endif /* !TFiller */

