#include "trpcraw.h"

ClassImp(TRpcRaw)

TRpcRaw::TRpcRaw( ) {
    fTrbnum  = -1;
    fCell    = -1;
    fCol     = -1;
    fRow     = -1;
    fX       = -1;
    fY       = -1;
    fZ       = -1;
    fTime    = -1;
    fCharge  = -1;
}
void TRpcRaw::setHit(Int_t trbnum,Int_t  cell,Int_t  col, Int_t  raw, Float_t x, Float_t y, Float_t z, Float_t time, Float_t charge) {
    fTrbnum  =   trbnum;
    fCell    =   cell;
    fCol     =   col;
    fRow     =   raw;
    fX       =   x;
    fY       =   y;
    fZ       =   z;
    fTime    =   time;
    fCharge  =   charge;
}

void TRpcRaw::getHit(Int_t& trbnum, Int_t&  cell,Int_t&  col,Int_t&  raw, Float_t& x, Float_t& y, Float_t& z, Float_t& time, Float_t& charge) {
    trbnum = fTrbnum;
    cell   = fCell;
    col    = fCol;
    raw    = fRow;
    x      = fX;
    y      = fY;
    z      = fZ;
    time   = fTime;
    charge = fCharge;
}


