#include "trpchitf.h"
#include "trpcraw.h"
#include "trpchit.h"
#include "TRandom3.h"
#include "stdlib.h"
#include <math.h>
#include "TMath.h"

#include <iostream>
using namespace std;

ClassImp(TRpcHitF);

TRpcHitF::TRpcHitF() {
    fPar        = new TRpcCalPar();
    //fRpcRawHits    = getRpcRawHits();
    fActiveCells = new TActiveCells();

    fRpcRawHits    = gEvent->getRpcRawHits();
    fRpcHitHits = new TClonesArray("TRpcHit",1000);
}

Int_t TRpcHitF::init(TString filename, TString filenameCells) {
    fPar        = new TRpcCalPar(filename);
    //fRpcRawHits    = getRpcRawHits();
    fActiveCells = new TActiveCells(filenameCells);

    fRpcRawHits    = gEvent->getRpcRawHits();
    if(!fRpcHitHits) {
        cout<<fRpcRawHits<<endl;
        fRpcHitHits = new TClonesArray("TRpcHit",1000);
    }

    return 1;
}
TRpcHitF::~TRpcHitF() {
    delete fPar;
    delete fRpcHitHits;
    delete fActiveCells;
}
Int_t TRpcHitF::execute() {
    // Clear the TClones array for later saving it!
    // clearAll();
    totalNHits = 0;
    //cout<<"point fRpcRawHits "<<fRpcRawHits<<endl;
    fRpcHitHits->Clear("C");
    Int_t cell, col, row, trbnum;
    Float_t     x, y, z,charge,time;
    // Loop through the TrbHits
    // Int_t ntrb = fRpcRawHits->GetEntriesFast();
    fRpcRawHits    = gEvent->getRpcRawHits();
    if(!fRpcRawHits) return 1;
    Int_t nraw = fRpcRawHits->GetEntriesFast();
    Int_t m1=0;
    Int_t m2=0;
    Int_t m3=0;

    Float_t q1=0.;
    Float_t q2=0.;
    Float_t q3=0.;

    TRandom3 ran;
    ran.SetSeed(0);

    Int_t year = gEvent->getEvtYear()+1900;
    //cout<<"year "<<year<<endl;
    Float_t doy  = gEvent->getDOY();

    Float_t sortt[360] = {0.};
    Int_t sorti[360]   = {0};


    for(Int_t i=0;i<nraw;i++) {
        TRpcRaw* hitraw = (TRpcRaw*)fRpcRawHits->At(i);
        if(!hitraw) continue;
        hitraw->getHit(trbnum, cell,col,row,x,y,z,time,charge);
        time   = time*0.098   - fPar->getTimeCal(trbnum,col-1,row-1);
        charge = charge*0.098 - fPar->getChargePedestal(trbnum,col-1,row-1);

        //cout<<"trbnum  "<<trbnum<<" col "<<col-1<<" row "<<row-1<<" par "<<fPar->getTimeCal(trbnum, col-1, row-1)<<endl;


        sortt[i] = time;

        if(charge<-0.)continue;

        //cout<<" trbnum "<<trbnum<<endl;
        if(!fActiveCells ->isCellActive(year, doy, trbnum, (col-1)*10+row-1)) continue;


        //cout<<"tcal "<<fPar->getTimeCal(trbnum,col-1,row-1)<<endl;
        //cout<<"pedestal "<<fPar->getChargePedestal(trbnum,col-1,row-1)<<endl;
        // q2w  HARDCODED PARAMETERS!
        Float_t q  =  (2.9454E-04  +
                       4.7284E+01*charge  +
                       -9.3629E-1*charge*charge +
                       1.4674E-2*charge*charge*charge  +
                       -9.4780E-5*charge*charge*charge*charge +
                       2.3987E-7*charge*charge*charge*charge*charge);
        //if( charge>0 ) {

        Float_t wcorr = 0;
        if(q<150)  wcorr = exp(q*-0.0111727+0.935106)+0.628899-0.0478589;
        if(q>=150) wcorr = exp(q*-0.00373208+0.63251);

        if(trbnum==2) {
            wcorr = exp(q*-7.18517e-04+7.87721e-01);
        }


//        if(q<140) wcorr = 3.69656e+02*pow(q+2.58017e+02,-7.77774e-01)-2.43326;
//        else if(q<600) wcorr = 1.98557e+02*pow(q+2.11211e+02,-7.98902e-01)-7.42063e-01;
//        else wcorr = 1.06623e+02*pow(q+6.83683e+02,-5.28285e-01)-2.23882e+00;


        // Remove Walk Correction for this production!
        time -= wcorr;




        TRpcHit* fHit = addRpcHit();
        Float_t randx  = (.5 - ran.Rndm())*116.;
        Float_t randy  = (.5 - ran.Rndm())*111.;
        // Uncomment next 2 lines if you want a simulation-like result.
        //Float_t randx  = 0.;
        //Float_t randy  = 0.;
        fHit->setHit(trbnum, cell, col, row, x+randx, y+randy, z, time, q);
        if(trbnum==0) {
            m1++;
            q1+=charge;
        }
        if(trbnum==1) {
            m2++;
            q2+=charge;
        }
        if(trbnum==2) {
            m3++;
            q3+=charge;
        }



        //}
    }

    //cout<<"Event nhits : "<<nraw<<endl;
    //TMath::Sort(nraw,sortt,sorti,kTRUE);
    /*
    //if(m1>1 && m2>1) {
    
    for(Int_t i=0;i<nraw;i++) {
        Int_t ind = sorti[i];
        TRpcRaw* hitraw = (TRpcRaw*)fRpcRawHits->At(ind);
        if(!hitraw) continue;
        hitraw->getHit(trbnum, cell,col,row,x,y,z,time,charge);
        charge = charge*0.098 - fPar->getChargePedestal(trbnum,col-1,row-1);
        cout<<"  "<<ind<<" trb "<<trbnum<<" c "<<col<<" r "<<row<<" t "<<time*0.098-3700.<<" q "<<charge<<endl;

        time   = time*0.098   - fPar->getTimeCal(trbnum,col-1,row-1);
    }
    }
    */

    gEvent->setRpcHits(fRpcHitHits);
    gEvent->setMults(m1,m2,m3);
    gEvent->setCharges(q1,q2,q3);
    return 1;

}
TRpcHit* TRpcHitF::addRpcHit( ) {
    TClonesArray& hits = *fRpcHitHits;
    TRpcHit *rpchit = new (hits[totalNHits++]) TRpcHit();
    //cout<<"ENTRIES in addhit "<<Hits->GetEntriesFast()<<endl;
    return rpchit;
}
