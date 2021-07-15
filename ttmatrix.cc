#include "ttmatrix.h"
#include "TMatrixF.h"

ClassImp(TTMatrix)

TTMatrix::TTMatrix( ) {
	fX1 = -1;
	fY1 = -1;
	fT1 = -1;
	fZ1 = -1;
	fX2 = -1;
	fY2 = -1;
	fT2 = -1;
	fZ2 = -1;

	fwx  = -1;
	fwy  = -1;
	fwz  = -1;
}

TMatrixF TTMatrix::InputSaeta2Planes(Float_t x1, Float_t y1, Float_t t1, Float_t z1, Float_t x2, Float_t y2, Float_t t2, Float_t z2){
	TMatrixF S2(6,1);
	fX1 = x1;
	fY1 = y1;
	fT1 = t1;
	fZ1 = z1;
	fX2 = x2;
	fY2 = y2;
	fT2 = t2;
	fZ2 = z2;
	Float_t Dz = z1-z2;
	S2[0][0] = (x2*z1-x1*z2)/Dz;
	S2[1][0] = (x1-x2)/Dz;
	S2[2][0] = (y2*z1-y1*z2)/Dz;
	S2[3][0] = (y1-y2)/Dz;
	S2[4][0] = (t2*z1-t1*z2)/Dz;
	S2[5][0] = (t1-t2)/Dz;
	return S2;
}

TMatrixF TTMatrix::AVector(TMatrixF SIn, Float_t x, Float_t y, Float_t t, Float_t z) {
	TMatrixF A(6,1);
	Float_t X0=SIn[0][0], XP=SIn[1][0], Y0=SIn[2][0], YP=SIn[3][0], T0=SIn[4][0], S=SIn[5][0], k=TMath::Sqrt(1.0+XP*XP+YP*YP);

	Float_t wx = 0;
	Float_t wy = 0;
	Float_t wz = 0;
	Float_t wt = 0;

	A[0][0] = wx*x;
	A[1][0] = z*(wx*x*k*k+S*wt*XP*(t*k+S*(XP*XP+YP*YP)*z))/(k*k);
	A[2][0] = wy*y;
	A[3][0] = z*(wy*y*k*k+S*wt*YP*(t*k+S*(XP*XP+YP*YP)*z))/(k*k);
	A[4][0] = wt*(t*k+S*(XP*XP+YP*YP)*z)/k;
	A[5][0] = wt*z*(t*k+S*(XP*XP+YP*YP)*z);
	return A;
}

TMatrixF TTMatrix::KMatrix(TMatrixF SIn, Float_t z){
        TMatrixF K(6,6);
        Float_t X0=SIn[0][0], XP=SIn[1][0], Y0=SIn[2][0], YP=SIn[3][0], T0=SIn[4][0], S=SIn[5][0], k=TMath::Sqrt(1.0+XP*XP+YP*YP);

	Float_t wx = 0;
	Float_t wy = 0;
	Float_t wz = 0;
	Float_t wt = 0;   
     
	K[0][0] = wx;
        K[0][1] = wx*z;
        K[1][0] = K[0][1];
        K[1][1] = wx*z*z+S*S*wt*XP*XP*z*z/(k*k);
        K[1][3] = S*S*wt*XP*YP*z*z/(k*k);
        K[1][4] = S*wt*XP*z/k;
        K[1][5] = S*wt*XP*z*z;
        K[2][2] = wy;
        K[2][3] = wy*z;
        K[3][1] = K[1][3];
        K[3][2] = K[2][3];
        K[3][3] = wy*z*z+S*S*wt*YP*YP*z*z/(k*k);
        K[3][4] = S*wt*YP*z/k;
        K[3][5] = S*wt*YP*z*z;
        K[4][1] = K[1][4];
        K[4][3] = K[3][4];
        K[4][4] = wt;
        K[4][5] = wt*k*z;
        K[5][1] = K[1][5];
        K[5][3] = K[3][5];
        K[5][4] = K[4][5];
        K[5][5] = wt*k*k*z*z;
        return K;
}

