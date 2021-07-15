#ifndef TTRACKF
#define TTRACKF

/////////////////////////////////////////
// Class TRpcRawF for finding and filling
// the RPC raw category
// Also fills the Syncronism signalin a special object
/////////////////////////////////////////
#include "trpchit.h"
#include "ttrack.h"
#include "TString.h"
#include "TClonesArray.h"
#include "tevent.h"
#include "TObject.h"
#include "trpccalpar.h"
#include "TMatrixF.h"

R__EXTERN Event *gEvent;

class TRpcHit;
class TClonesArray;

class TTrackF: public TObject
{
private:
    TClonesArray    *fRpcHitHits;
    TClonesArray    *fRpcHitCorr;
    TClonesArray    *fTracks;
    TClonesArray    *fTracks3Planes;    
    Int_t totalNHits;
    Int_t totalNHits3Planes;
    Int_t totalNHitsCorr;
public:
    TTrackF();
    ~TTrackF();
    Int_t execute();
    Int_t init();
    TTrack* addTrack();
    TRpcHit* addRpcHit();
    TClonesArray*  getTracks() {return fTracks;}
    TClonesArray*  getTracks3Planes(){return fTracks3Planes;}
    TClonesArray*  getRpcHitCorr() {return fRpcHitCorr;}
    TClonesArray*  getRpcHits()    {return gEvent->getRpcHits();}
    TMatrixF AVector(TMatrixF, Float_t, Float_t, Float_t, Float_t);
    TMatrixF KMatrix(TMatrixF, Float_t);
    TMatrixF Saeta2Planes(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
    ClassDef(TTrackF,1);
};

#endif /* !TTRACKF */
