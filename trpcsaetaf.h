#ifndef TRPCSAETAF
#define TRPCSAETAF

#include "trpchit.h"
#include "trpcsaeta.h"
#include "TString.h"
#include "TClonesArray.h"
#include "tevent.h"
#include "TObject.h"
#include "trpccalpar.h"
#include "TMatrixF.h"

R__EXTERN Event *gEvent;

class TRpcSaeta;
class TClonesArray;

class TRpcSaetaF: public TObject
{
private:
	TClonesArray  *fRpcHitHits;
	TClonesArray  *fRpcHitCorr;
	TClonesArray  *fRpcSaeta2Planes;
	TClonesArray  *fRpcSaeta3Planes;
	Int_t totalNHits;
	Int_t totalNHits2Planes;
	Int_t totalNHits3Planes;
	Int_t totalNHitsCorr;
public:
	TRpcSaetaF();
	~TRpcSaetaF();
	Int_t execute();
	Int_t init();
	TRpcHit*   addRpcHit();
	TRpcSaeta* addRpcSaeta2Planes();
	TRpcSaeta* addRpcSaeta3Planes();

	TClonesArray*  getRpcHits()         {return gEvent->getRpcHits(); }
	TClonesArray*  getRpcHitCorr()      {return fRpcHitCorr;          }
	TClonesArray*  getRpcSaeta2Planes() {return fRpcSaeta2Planes;     }
	TClonesArray*  getRpcSaeta3Planes() {return fRpcSaeta3Planes;     }

    	TMatrixF InputSaeta2Planes(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
   	TMatrixF KMatrix(TMatrixF, Float_t);
    	TMatrixF AVector(TMatrixF, Float_t, Float_t, Float_t, Float_t);
	
	ClassDef(TRpcSaetaF,1);
};

#endif /* !TRPCSAETAF */
