#ifndef TRPCRAW
#define TRPCRAW

#include "TObject.h"

class TRpcRaw : public TObject {
protected:
    Int_t   fTrbnum;
    Int_t   fCell;
    Int_t   fCol;
    Int_t   fRow;
    Float_t fX;
    Float_t fY;
    Float_t fZ;
    Float_t fTime;
    Float_t fCharge;
public:
    TRpcRaw();
    ~TRpcRaw(){}
    void setHit(Int_t trbnum,Int_t  cell,Int_t  col,Int_t  row, Float_t x, Float_t y, Float_t z, Float_t time, Float_t charge);
    void getHit(Int_t& trbnum, Int_t& cell,Int_t&  col, Int_t&  row, Float_t& x, Float_t& y, Float_t& z, Float_t& time, Float_t& charge);
    ClassDef(TRpcRaw,1)
};

#endif /* !TRPRAW */
