#ifndef TRPCRAWF
#define TRPCRAWF

/////////////////////////////////////////
// Class TRpcRawF for finding and filling
// the RPC raw category
// Also fills the Syncronism signalin a special object
/////////////////////////////////////////
#include "trpcraw.h"
#include "trpclookuptable.h"
#include "TString.h"
#include "TClonesArray.h"
#include "tevent.h"
#include "TObject.h"

R__EXTERN Event *gEvent;

class TRpcLookupTable;
class TRpcRaw;
class TClonesArray;


class TRpcRawF: public TObject
{
private:
    RpcLookupTable *fLookPar;
    TClonesArray    *fTrbHits;
    TClonesArray    *fRpcRawHits;
    Int_t totalNHits;
public:
    TRpcRawF();
    ~TRpcRawF();
    Int_t execute();
    Int_t init(TString filename);
    TRpcRaw* addRpcRaw();
    TClonesArray*  getTrbRawHits() {return gEvent->getHits();}
    TClonesArray*  getRpcRawHits() {return fRpcRawHits;}
    ClassDef(TRpcRawF,0);
};

#endif /* !TRPRAWF */
