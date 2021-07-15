#ifndef TTMATRIX
#define TTMATRIX

#include "TObject.h"
#include "TMath.h"
#include "TMatrixF.h"

class TTMatrix : public TObject {

protected:
	Float_t fX1;
	Float_t fY1;
	Float_t fT1;
	Float_t fZ1;
	Float_t fX2;
	Float_t fY2;
	Float_t fT2;
	Float_t fZ2;

	Float_t fwx;
	Float_t fwy;
	Float_t fwz;

public:
    TTMatrix();
    ~TTMatrix(){}

    TMatrixF AVector(TMatrixF, Float_t, Float_t, Float_t, Float_t);
    TMatrixF KMatrix(TMatrixF, Float_t);
    TMatrixF InputSaeta2Planes(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
    //void InputSaeta3Planes(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);

    ClassDef(TTMatrix,1)
};

#endif /* !TTMATRIX */
