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
   TH2D* h_eff[2];      // 2D efficiency along the X coordinate
   TH2D* h_flux3[6][6]; // The fifth will be "all" - 6 theta intervals!
   TH1D* h_daq_active;
   TH1D* h_daq_sync;
   TH2D* h_cell_counts[3];

   TH1D* h_eff_nn[2];
   TH2D* h_cell_counts_n1[3];
   TH2D* h_cell_counts_n2[3];
   TH2D* h_n_q[3];

   TH2D* h_n_q_clean[3];


   TH2D* h_flux2[6][6];

   TH2D* h_y_x_theta[3][12];

   TH2D* h_log10q_flux[3];
   TH2D* h_theta_time[2];





    /*
     TH1D *hq[3*10*12];
     TH1D *hdta[3*10*12];
     TH1D* hdt[10*12*10*12+1];
     TH1D* h1D[10];
     TH2D* h2D[8];
     TH3D* h3D[6];
     */







};

#endif /* !TFiller */

