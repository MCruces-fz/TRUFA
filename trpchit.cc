#include "trpchit.h"

ClassImp(TRpcHit)

TRpcHit::TRpcHit( ) {
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

void TRpcHit::setHit(Int_t trbnum,Int_t  cell,Int_t  col, Int_t  row, Float_t x, Float_t y, Float_t z, Float_t time, Float_t charge) {
    fTrbnum  =   trbnum;
    fCell    =   cell;
    fCol     =   col;
    fRow     =   row;
    fX       =   x;
    fY       =   y;
    fZ       =   z;
    fTime    =   time;
    fCharge  =   charge;
}

void TRpcHit::getHit(Int_t& trbnum, Int_t&  cell,Int_t&  col,Int_t&  row, Float_t& x, Float_t& y, Float_t& z, Float_t& time, Float_t& charge) {
    trbnum = fTrbnum;
    cell   = fCell;
    col    = fCol;
    row    = fRow;
    x      = fX;
    y      = fY;
    z      = fZ;
    time   = fTime;
    charge = fCharge;
}


