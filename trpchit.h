#ifndef TRPCHIT
#define TRPCHIT

#include "TObject.h"

class TRpcHit : public TObject {
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
    TRpcHit();
    ~TRpcHit(){}

    Int_t getTrbnum()   {return fTrbnum;}
    Int_t getCell()     {return fCell;  }
    Int_t getCol()      {return fCol;   }
    Int_t getRow()      {return fRow;   }
    Float_t getX()      {return fX;     }
    Float_t getY()      {return fY;     }
    Float_t getZ()      {return fZ;     }
    Float_t getTime()   {return fTime;  }
    Float_t getCharge() {return fCharge;}

    void setTrbnum(Int_t num   ) {fTrbnum = num;}
    void setCell(Int_t  num  )   {fCell = num;  }
    void setCol(Int_t num   )    {fCol = num;   }
    void setRow(Int_t num   )    {fRow = num;   }
    void setX(Float_t val )      {fX = val;     }
    void setY(Float_t val )      {fY = val;     }
    void setZ(Float_t val )      {fZ = val;     }
    void setTime(Float_t val )   {fTime = val;  }
    void setCharge(Float_t val ) {fCharge = val;}

    void setHit(Int_t trbnum,Int_t  cell,Int_t  col,Int_t  row, Float_t x, Float_t y, Float_t z, Float_t time, Float_t charge);
    void getHit(Int_t& trbnum, Int_t& cell,Int_t&  col, Int_t&  row, Float_t& x, Float_t& y, Float_t& z, Float_t& time, Float_t& charge);


    ClassDef(TRpcHit,1)
};

#endif /* !TRPHIT */
