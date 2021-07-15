#include "ttrack.h"

ClassImp(TTrack)

TTrack::TTrack( ) {
    find0  = -1;
    find1  = -1;
//    find2  = -1;
    fAl    = -1;
    fBe    = -1;
    fGa    = -1;
    fX     = -1;
    fY     = -1;
    fZ     = -1;
    fTime  = -1;
}

void TTrack::setTrack(Float_t x0,Float_t y0,Float_t t0,Float_t al,Float_t be,Float_t ga,Int_t ind0,Int_t ind1) {
    find0  =   ind0;
    find1  =   ind1;
//    find2  =   ind2;
    fAl    =   al;
    fBe    =   be;
    fGa    =   ga;
    fX     =   x0;
    fY     =   y0;
    fZ     =   1000.;
    fTime  =   t0;
}

