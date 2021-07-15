#ifndef TRPCCALPAR
#define TRPCCALPAR

#include "TObject.h"
#include "TString.h"

//using namespace std;
class TRpcCalPar: public TObject
{
protected:
    Float_t pars[3][10][12][2];                //!- 4 detectors x 4 Tdc trb x 32 chan x 6
public:
    TRpcCalPar();
    TRpcCalPar(TString filename);
    ~TRpcCalPar(){}
    Float_t getTimeCal(Int_t trbnum, Int_t col,Int_t row) {return pars[trbnum][row][col][1];}
    Float_t getChargePedestal(Int_t trbnum, Int_t col,Int_t row) {return pars[trbnum][row][col][0];}
    Bool_t readParams(TString filename);
    Bool_t writeParams(TString filename);
public:
    ClassDef(TRpcCalPar,0)
};

#endif /* !TRPCCALPAR */