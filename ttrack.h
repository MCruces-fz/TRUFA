#ifndef TTRACK
#define TTRACK

#include "TObject.h"
#include "TMath.h"

class TTrack : public TObject {
protected:
    Int_t   find0;
    Int_t   find1;
//    Int_t   find2;
    Float_t   fAl;
    Float_t   fBe;
    Float_t   fGa;
    Float_t fX;
    Float_t fY;
    Float_t fZ;
    Float_t fTime;

public:
    TTrack();
    ~TTrack(){}

    Int_t getInd(Int_t n)       {
	if(n==0) return find0;
	if(n==1) return find1;
//	if(n==2) return find2;
        return -1;
    }
    Float_t getX0()      {return fX;     }
    Float_t getY0()      {return fY;     }
    Float_t getZ0()      {return fZ;     }
    Float_t getTime()    {return fTime;  }
    Float_t getAl()      {return fAl;    }
    Float_t getBe()      {return fBe;    }
    Float_t getGa()      {return fGa;    }
    Float_t getPhi() {return atan2(fAl,fBe);}
    Float_t getPhiDeg() {return atan2(fAl,fBe)*TMath::RadToDeg();}
    Float_t getTheta() {return acos(fGa);}
    Float_t getThetaDeg() {return acos(fGa)*TMath::RadToDeg();}


    void setTime(Float_t val )   {fTime = val;  }
    void setTrack(Float_t x0,Float_t y0,Float_t t0,Float_t al,Float_t be,Float_t ga,Int_t ind0,Int_t ind1);
    void getTrack(Float_t& x0,Float_t& y0,Float_t& t0,Float_t& al,Float_t& be,Float_t& ga,Int_t& ind0,Int_t& ind1){
	x0 = fX;
	y0 = fY;
	t0 = fTime;
	al = fAl;
	be = fBe;
	ga = fGa;
        ind0 = find0;
        ind1 = find1;
    }


    ClassDef(TTrack,1)
};

#endif /* !TTRACK */
