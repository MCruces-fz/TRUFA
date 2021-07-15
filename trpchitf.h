#ifndef TRPCHITF
#define TRPCHITF

/////////////////////////////////////////
// Class TRpcRawF for finding and filling
// the RPC raw category
// Also fills the Syncronism signalin a special object
/////////////////////////////////////////
#include "trpchit.h"
#include "TString.h"
#include "TClonesArray.h"
#include "tevent.h"
#include "TObject.h"
#include "trpccalpar.h"
#include "tactivecells.h"

R__EXTERN Event *gEvent;

class TRpcRaw;
class TClonesArray;

class TRpcHitF: public TObject
{
private:
    TRpcCalPar      *fPar;
    TClonesArray    *fRpcRawHits;
    TClonesArray    *fRpcHitHits;
    TActiveCells    *fActiveCells;
    Int_t totalNHits;
public:
    TRpcHitF();
    ~TRpcHitF();
    Int_t execute();
    Int_t init(TString filename, TString filenameActive);
    TRpcHit* addRpcHit();
    TClonesArray*  getRpcRawHits() {return gEvent->getRpcRawHits();}
    TClonesArray*  getRpcHits()    {return fRpcHitHits;}
    ClassDef(TRpcHitF,0);
};

#endif /* !TRPHITF */
