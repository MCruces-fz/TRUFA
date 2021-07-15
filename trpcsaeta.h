#ifndef TRPCSAETA
#define TRPCSAETA

#include "TObject.h"
#include "TMath.h"

class TRpcSaeta : public TObject {

protected:
	Float_t fX;	// X coordinate
	Float_t fXP;	// X slope
	Float_t fY;	// Y coordinate
	Float_t fYP;	// Y slope
	Float_t fZ;	// Z coordinate

	Float_t fTime;	// Track time
	Float_t fSlow;	// Slowness -> 1/v

	Float_t fAl;	// Alpha angle
	Float_t fBe;	// Beta angle
	Float_t fGa;	// Gamma angle

	Int_t   fSaN;	// Saeta order
	Int_t   find0;	// Hit index
	Int_t   find1;	// Hit index
	Int_t   find2;	// Hit index

	Float_t fChi2;	// Chi-square

public:
	TRpcSaeta();
	~TRpcSaeta(){}

	Float_t getX0()	{return fX;   }
	Float_t getXP()	{return fXP;  }
	Float_t getY0()	{return fY;   }
	Float_t getYP()	{return fYP;  }
	Float_t getZ0()	{return fZ;   }

	Float_t getTime(){return fTime;}
	Float_t getSlow(){return fSlow;}

	Float_t getAl()		{return fAl;  }
	Float_t getBe()		{return fBe;  }
	Float_t getGa()		{return fGa;  }
	Float_t getPhi()	{return atan2(fAl,fBe);}
	Float_t getTheta()	{return acos(fGa);     }
	Float_t getPhiDeg()	{return atan2(fAl,fBe)*TMath::RadToDeg();}
	Float_t getThetaDeg()	{return acos(fGa)*TMath::RadToDeg();     }

	Int_t   getSaetaN(){return fSaN; }

	Int_t getInd(Int_t n)       {
	if(n==0) return find0;
	if(n==1) return find1;
	if(n==2) return find2;
	return -1;
	}

	Float_t getChi2()    {return fChi2;}

	// setRpcSaeta2Planes is not a real Saeta -> Used to compare new & previous results only. Commented the real saeta definition

	void setTime(Float_t val) {fTime = val;}
	//void setRpcSaeta2Planes(Float_t x0,Float_t xP,Float_t y0,Float_t yP,Float_t z, Float_t t0,Float_t sl, Float_t al,Float_t be,Float_t ga, Int_t san,Int_t ind0,Int_t ind1, Float_t chi2);
	void setRpcSaeta2Planes(Float_t x0,Float_t y0,Float_t t0,Float_t al,Float_t be,Float_t ga,Int_t ind0,Int_t ind1);
	void setRpcSaeta3Planes(Float_t x0,Float_t xP,Float_t y0,Float_t yP,Float_t z, Float_t t0,Float_t sl, Float_t al,Float_t be,Float_t ga, Int_t san,Int_t ind0,Int_t ind1,Int_t ind2, Float_t chi2);

	//void getRpcSaeta2Planes(Float_t& x0,Float_t& xP,Float_t& y0,Float_t& yP,Float_t& z, Float_t& t0,Float_t& sl, Float_t& al,Float_t& be,Float_t& ga, Int_t& san,Int_t& ind0,Int_t& ind1, Float_t& chi2);
	void getRpcSaeta2Planes(Float_t& x0,Float_t& y0,Float_t& t0,Float_t& al,Float_t& be,Float_t& ga,Int_t& ind0,Int_t& ind1);
	void getRpcSaeta3Planes(Float_t& x0,Float_t& xP,Float_t& y0,Float_t& yP,Float_t& z, Float_t& t0,Float_t& sl, Float_t& al,Float_t& be,Float_t& ga, Int_t& san,Int_t& ind0,Int_t& ind1,Int_t& ind2, Float_t& chi2);

	ClassDef(TRpcSaeta,1)
};

#endif /* !TRPCSAETA */
