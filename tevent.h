
#ifndef TEVENT
#define TEVENT

/////////////////////////////////////////////////////////////////////
// Event
//
// Event is a class used for creating the root tree.
// It contains the information unpacked from hld files.
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TBuffer.h>
#include <TClonesArray.h>
#include "teventhdr.h"
#include "TGlobal.h"
#include "TMath.h"
#include "TTimeStamp.h"


#include <iostream>
using namespace std;


typedef UInt_t UInt4;
typedef UChar_t UInt1;

class HldEvent;
class Hit;
class Event;


class Event :public TObject {
public:
    Event();
    Event(const HldEvent& HldEvt, Int_t refCh);
    virtual ~Event();
protected:
    EventHdr EvtHdr;      //contains data from event header
    const Int_t kMaxMult; //Maximum number of multiplicity
    const Int_t kMaxChannelNr;

    //data from subEvent
    UInt_t subEvtId;      //Id of subEvent according to HLD format
private:

    Int_t errors_per_event; 	// number of TDC errors per event
    Int_t referenceChannel;	// reference channel
    Int_t referenceTime;	// reference time from indicated channel if it is used

    Int_t totalNHits;		// total number of hits per one event

    TClonesArray* Hits;        		 //! * array with all hits
    TClonesArray* rpcRaw;      		 //! * array with rpcRaw
    TClonesArray* rpcHit;      		 //! * array with rpcHits
    TClonesArray* rpcHitCorr;   	 //! * array with rpcHits
    TClonesArray* RpcSaeta2Planes;       //! * array with tracks
    TClonesArray* RpcSaeta3Planes;       //! * array with tracks

    Int_t multH1;
    Int_t multH2;
    Int_t multH3;
    Float_t chargeH1;
    Float_t chargeH2;
    Float_t chargeH3;

    Int_t multT;
    Int_t multT3;

    Bool_t sync;
    // Give access from everywhere
    static TClonesArray* gHits;

    // THIS PART SHOULD BE IN A DIFFERENT PLACE IS JUST A FAST SOLUTION
    // UNTIL WE HAVE A PROPERLY DEFINED RUNTIMEDATABASE!
    // initialized only once!

    TString fileLookupPar;         //! * name of the file with lookup params
    TString fileHitFinderPar;      //! * name of the file with hit finder params
    TString fileTrackFinderPar;    //! * name of the file with track finder params

    // incident angle of the events!
    Float_t fEvAl;
    Float_t fEvBe;
    Float_t fEvGa;

    Float_t fEvAl3Planes = -100;
    Float_t fEvBe3Planes= -100;
    Float_t fEvGa3Planes= -100;

    TTimeStamp time;             //! * time stamp for calculation of the time in universal format

public:
    Int_t getEvtSize()  { return EvtHdr.getSize(); }
    Int_t getEvtDecoding()  { return EvtHdr.getDecoding(); }
    Int_t getEvtId()  { return EvtHdr.getId(); }
    Int_t getEvtSeqNr()  { return EvtHdr.getSeqNr(); }
    Int_t getEvtDate()  { return EvtHdr.getDate(); }
    Int_t getEvtTime()  { return EvtHdr.getTime(); }
    Int_t getEvtYear()  { return EvtHdr.getYear(); }
    Int_t getEvtMonth()  { return EvtHdr.getMonth(); }
    Int_t getEvtDay()  { return EvtHdr.getDay(); }
    Int_t getEvtHour()  { return EvtHdr.getHour(); }
    Int_t getEvtMinute()  { return EvtHdr.getMinute(); }
    Int_t getEvtSecond()  { return EvtHdr.getSecond(); }
    Float_t getDOY();
    Int_t getEvtPad()  { return EvtHdr.getPad(); }
    Int_t getEvtDataSize() { return EvtHdr.getDataSize();}
    Int_t getEvtPaddedSize() { return EvtHdr.getPaddedSize(); }
    Int_t getSubEvtId() { return subEvtId; }
    Int_t getReferenceTime() {return referenceTime; }
    Int_t getReferenceChannel() {return referenceChannel; }
    Int_t getErrors_per_event() {return errors_per_event; }
    Bool_t getSync() {return sync;}
    Int_t getMultH1() {return multH1;}
    Int_t getMultH2() {return multH2;}
    Int_t getMultH3() {return multH3;}
    Float_t getChargeH1() {return chargeH1;}
    Float_t getChargeH2() {return chargeH2;}
    Float_t getChargeH3() {return chargeH3;}

    Int_t getNtracks() {return multT;}
    Int_t getNtracks3() {return multT3;}

    TClonesArray* getHits(){return gHits;}
    TClonesArray* getRpcRawHits(){return rpcRaw;}
    TClonesArray* getRpcHits(){return rpcHit;}
    TClonesArray* getRpcHitsCorr(){return rpcHitCorr;}
    TClonesArray* getRpcSaeta2Planes(){return RpcSaeta2Planes;}
    TClonesArray* getRpcSaeta3Planes(){return RpcSaeta3Planes;}

    void setRpcRawHits(TClonesArray* array){rpcRaw = array;}
    void setRpcHits(TClonesArray* array)   {rpcHit = array;}
    void setRpcHitsCorr(TClonesArray* array)   {rpcHitCorr = array;}
    void setRpcSaeta2Planes(TClonesArray* array)   {RpcSaeta2Planes = array;}
    void setRpcSaeta3Planes(TClonesArray* array)   {RpcSaeta3Planes = array;}

    void setNtracks(Int_t n) {multT = n;}
    void setNtracks3(Int_t n) {multT3 = n;}

    void setErrors_per_event(Int_t e) { errors_per_event=e; }
    void setReferenceChannel(Int_t c) { referenceChannel=c; }
    void setReferenceTime(Int_t t) { referenceTime=t; }
    void setSubEvtId(Int_t sb=0) { subEvtId=sb; }
    void setNHits(Int_t n) {totalNHits = n;}
    void setMults(Int_t m1,Int_t m2,Int_t m3) {multH1=m1;multH2=m2;multH3=m3;}
    void setCharges(Int_t q1,Int_t q2,Int_t q3) {chargeH1=q1;chargeH2=q2;chargeH3=q3;}
    void setSync(Bool_t s) {sync=s;}

    void setFileLookupPar( TString fileName) {         fileLookupPar      = fileName;}
    void setFileHitFinderPar( TString fileName) {      fileHitFinderPar   = fileName;}
    void setFileTrackFinderPar( TString fileName) {    fileTrackFinderPar = fileName;}

    TString setFileLookupPar( void ) {        return fileLookupPar;}
    TString setFileHitFinderPar( void ) {     return fileHitFinderPar;}
    TString setFileTrackFinderPar( void ) {   return fileTrackFinderPar;}

    void setAngles(Float_t al,Float_t be,Float_t ga) {
        fEvAl = al;
        fEvBe = be;
        fEvGa = ga;
    }
    void setAngles3(Float_t al,Float_t be,Float_t ga) {
        fEvAl3Planes = al;
        fEvBe3Planes = be;
        fEvGa3Planes = ga;
    }

    Float_t getAl()      {return fEvAl;    }
    Float_t getBe()      {return fEvBe;    }
    Float_t getGa()      {return fEvGa;    }
    Float_t getAl3()     {return fEvAl3Planes;    }
    Float_t getBe3()     {return fEvBe3Planes;    }
    Float_t getGa3()     {return fEvGa3Planes;    }

    Float_t getPhi() {return atan2(fEvAl,fEvBe);}
    Float_t getPhiDeg() {return atan2(fEvAl,fEvBe)*TMath::RadToDeg();}
    Float_t getTheta() {return acos(fEvGa);}
    Float_t getThetaDeg() {return acos(fEvGa)*TMath::RadToDeg();}

    Float_t getPhi3() {return atan2(fEvAl3Planes,fEvBe3Planes);}
    Float_t getPhiDeg3() {return atan2(fEvAl3Planes,fEvBe3Planes)*TMath::RadToDeg();}
    Float_t getTheta3() {return acos(fEvGa3Planes);}
    Float_t getThetaDeg3() {return acos(fEvGa3Planes)*TMath::RadToDeg();}

    Hit* addHit();
    Bool_t fill(const HldEvent&);
    void clearAll(void);
public:
    ClassDef(Event,1);
};

#endif /* !TEVENT */
