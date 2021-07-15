#ifndef TRPCLOOKUPTABLE
#define TRPCLOOKUPTABLE

#include "TObject.h"
#include "TString.h"

//using namespace std;
class RpcLookupTable: public TObject
{
protected:
    Float_t pars[4][4][32][6];                //!- 4 detectors x 4 Tdc trb x 32 chan x 6
    Float_t detInd[4];                        //!- cell col raw x y z
public:
    RpcLookupTable();
    RpcLookupTable(TString filename);
    ~RpcLookupTable(){}
    void getParams(Int_t num,Int_t mbch, Int_t ch, Int_t& cell, Int_t& col, Int_t& row, Float_t& x, Float_t& y, Float_t& z);
    void setParams(Int_t num,Int_t mbch, Int_t ch, Int_t cell, Int_t col, Int_t row, Float_t x, Float_t y, Float_t z);
    Bool_t readParams(TString filename);
    Bool_t writeParams(TString filename);
    Int_t getDetNum(Int_t trbId);
public:
    ClassDef(RpcLookupTable,0)
};

#endif /* !TRPCLOOKUPTABLE */
