#include "trpcsaeta.h"

ClassImp(TRpcSaeta)

TRpcSaeta::TRpcSaeta( ) {
	fX    = -1;
	fXP   = -1;
	fY    = -1;
	fYP   = -1;
	fZ    = -1;
	fTime = -1;
	fSlow = -1;
	fAl   = -1;
	fBe   = -1;
	fGa   = -1;
	fSaN  = -1;
	find0 = -1;
	find1 = -1;
	find2 = -1;
	fChi2 = -1;
}

/*
void TRpcSaeta::setRpcSaeta2Planes(Float_t x0,Float_t xP,Float_t y0,Float_t yP,Float_t z, Float_t t0,Float_t sl, Float_t al,Float_t be,Float_t ga, Int_t san,Int_t ind0,Int_t ind1, Float_t chi2) {
	fX    = x0;
	fXP   = xP;
	fY    = y0;
	fYP   = yP;
	fZ    = 1000.;
	fTime = t0;
	fSlow = sl;
	fAl   = al;
	fBe   = be;
	fGa   = ga;
	fSaN  = san;
	find0 = ind0;
	find1 = ind1;
	fChi2 = chi2;
}
*/

void TRpcSaeta::setRpcSaeta2Planes(Float_t x0,Float_t y0, Float_t t0,Float_t al,Float_t be,Float_t ga, Int_t ind0,Int_t ind1) {
        fX    = x0;
        //fXP   = xP;
        fY    = y0;
        //fYP   = yP;
        //fZ    = 1000.;
        fTime = t0;
        //fSlow = sl;
        fAl   = al;
        fBe   = be;
        fGa   = ga;
        //fSaN  = san;
        find0 = ind0;
        find1 = ind1;
        //fChi2 = chi2;
}


void TRpcSaeta::setRpcSaeta3Planes(Float_t x0,Float_t xP,Float_t y0,Float_t yP,Float_t z, Float_t t0,Float_t sl, Float_t al,Float_t be,Float_t ga, Int_t san,Int_t ind0,Int_t ind1,Int_t ind2, Float_t chi2) {
	fX    = x0;
	fXP   = xP;
	fY    = y0;
	fYP   = yP;
	fZ    = 1000.;
	fTime = t0;
	fSlow = sl;
	fAl   = al;
	fBe   = be;
	fGa   = ga;
	fSaN  = san;
	find0 = ind0;
	find1 = ind1;
	find2 = ind2;
	fChi2 = chi2;
}

/*
void TRpcSaeta::getRpcSaeta2Planes(Float_t& x0,Float_t& xP,Float_t& y0,Float_t& yP,Float_t& z, Float_t& t0,Float_t& sl, Float_t& al,Float_t& be,Float_t& ga, Int_t& san,Int_t& ind0,Int_t& ind1, Float_t& chi2) {
	x0   = fX;
	xP   = fXP;
	y0   = fY;
	yP   = fYP;
	z    = fZ;
	t0   = fTime;
	sl   = fSlow;
	al   = fAl;
	be   = fBe;
	ga   = fGa;
	san  = fSaN;
	ind0 = find0;
	ind1 = find1;
	chi2 = fChi2;
}
*/

void TRpcSaeta::getRpcSaeta2Planes(Float_t& x0,Float_t& y0, Float_t& t0,Float_t& al,Float_t& be,Float_t& ga, Int_t& ind0,Int_t& ind1) {
        x0   = fX;
        y0   = fY;
        t0   = fTime;
        al   = fAl;
        be   = fBe;
        ga   = fGa;
        ind0 = find0;
        ind1 = find1;
}

void TRpcSaeta::getRpcSaeta3Planes(Float_t& x0,Float_t& xP,Float_t& y0,Float_t& yP,Float_t& z, Float_t& t0,Float_t& sl, Float_t& al,Float_t& be,Float_t& ga, Int_t& san,Int_t& ind0,Int_t& ind1,Int_t& ind2, Float_t& chi2) {
	x0   = fX;
	xP   = fXP;
	y0   = fY;
	yP   = fYP;
	z    = fZ;
	t0   = fTime;
	sl   = fSlow;
	al   = fAl;
	be   = fBe;
	ga   = fGa;
	san  = fSaN;
	ind0 = find0;
	ind1 = find1;
	ind2 = find2;
	chi2 = fChi2;
}

