#include "ttrackf.h"
#include "trpchitf.h"
#include "trpchit.h"
#include "ttrack.h"
#include "stdlib.h"
#include <math.h>

// Global constants: statistical weigth of the measured variables.
// // w=1/sigma^2
// // Distances in mm, time in ns 
Float_t wx=11.0, wy=12.0, wt=0.3;


ClassImp(TTrackF);

TTrackF::TTrackF() {
    fRpcHitHits    = gEvent->getRpcHits();
    fRpcHitCorr    = new TClonesArray("TRpcHit",1000);
    fTracks        = new TClonesArray("TTrack",1000);
}

Int_t TTrackF::init() {
    fRpcHitHits    = gEvent->getRpcHits();
    //if(!fRpcRawHits)cout<<fRpcRawHits<<endl;
    if(!fRpcHitCorr)
        fRpcHitCorr    = new TClonesArray("TRpcHit",1000);
    if(!fTracks)
        fTracks = new TClonesArray("TTrack",1000);
    return 1;
}
/*
TTrackF::~TTrackF() {
    delete fRpcHitCorr;
    delete fTracks;
}*/
Int_t TTrackF::execute() {
    // Clear the TClones array for later saving it!
    // clearAll();
    totalNHits = 0;
    fTracks     ->Clear("C");
    fRpcHitCorr ->Clear("C");

    Int_t maxIter = 100; // Maximum number of iterations until convergence
    Float_t mod, mod2;
    Float_t cutMod = 1e-5;

    Int_t   cell,  col,  row,  trbnum;
    Int_t   cell2, col2, row2, trbnum2;
    Float_t x,  y,  z,  charge,  time;
    Float_t x2, y2, z2, charge2, time2;
    Float_t kv, al, be, ga;
    Float_t evAl = 0.;
    Float_t evBe = 0.;
    Float_t evGa = 0.;
    Float_t xp, yp, x0, y0, t0, sl;

    Int_t ind0, ind1;
    Float_t thitmin = 1e6;

    TMatrixF SInit(6,1), S2Planes(6,1), S3Planes(6,1), A(6,1), K(6,6);

    fRpcHitHits    = gEvent->getRpcHits();
    if(!fRpcHitHits) return 1;
    Int_t nhits = fRpcHitHits->GetEntriesFast();


    for(Int_t i=0;i<fRpcHitHits->GetEntriesFast();i++) {
        TRpcHit* hit = (TRpcHit*)fRpcHitHits->At(i);
        if(!hit) continue;
        hit->getHit(trbnum, cell,col,row,x,y,z,time,charge);
        //cout<<time<<" "<<trbnum<<" "<<cell<<endl;
        if(trbnum!=0) continue;
        if(thitmin>time) thitmin = time;
    }

    //cout<<thitmin<<endl;

   for(Int_t i=0;i<nhits-1;i++) {
        TRpcHit* hit = (TRpcHit*)fRpcHitHits->At(i);
        if(!hit) continue;
        hit->getHit(trbnum, cell,col,row,x,y,z,time,charge);
        if(trbnum==2) continue;
        //loop in strict order in order to prevent double counting.
        for(Int_t j=i+1;j<nhits;j++) {
            TRpcHit* hit2 = (TRpcHit*)fRpcHitHits->At(j);
            if(!hit2) continue;
            hit2->getHit(trbnum2, cell2,col2,row2,x2,y2,z2,time2,charge2);
            if(trbnum2==trbnum) continue;
            if(trbnum2==2) continue;
            Float_t dist = sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2));
            //cout<<time<<" "<<time2<<" "<<dist/300.<<endl;
            if(trbnum==0 && trbnum2==1) {
                if( fabs(time2-time-dist/300.)>0.9) continue;
		    // Input Saeta
		    SInit = TTrackF::Saeta2Planes(x,y,time,z, x2,y2,time2,z2);
		    // Convergence loop
		    for (Int_t k=0;k<maxIter;k++){
			A = TTrackF::AVector(SInit, x,y,time,z) + TTrackF::AVector(SInit, x2,y2,time2,z2);
			K = TTrackF::KMatrix(SInit, z) + TTrackF::KMatrix(SInit, z2);
			mod  = TMath::Sqrt(SInit.E2Norm());
			S2Planes = K.Invert()*A;
			SInit    = S2Planes;
			mod2 = TMath::Sqrt(SInit.E2Norm());
			if (TMath::Abs(mod2-mod)<cutMod){break;} // Saeta has converged	
		    }
                ind0 = i;
                ind1 = j;
            }
            if(trbnum==1 && trbnum2==0) {
                if( fabs(time-time2-dist/300.)>0.9) continue;
		// Input Saeta
		SInit = TTrackF::Saeta2Planes(x2,y2,time2,z2, x,y,time,z);
		// Convergence loop
		for (Int_t k=0;k<maxIter;k++){
			A = TTrackF::AVector(SInit, x2,y2,time2,z2) + TTrackF::AVector(SInit, x,y,time,z);
			K = TTrackF::KMatrix(SInit, z2) + TTrackF::KMatrix(SInit, z);
			mod  = TMath::Sqrt(SInit.E2Norm());
			S2Planes = K.Invert()*A;
			SInit    = S2Planes;
			mod2 = TMath::Sqrt(SInit.E2Norm());
			if (TMath::Abs(mod-mod2)<cutMod){break;} // Saeta has converged	
		}
                ind0 = j;
                ind1 = i;

            }

	    x0 = S2Planes[0][0];
	    xp = S2Planes[1][0];
	    y0 = S2Planes[2][0];
 	    yp = S2Planes[3][0];
  	    t0 = S2Planes[4][0];
	    sl = S2Planes[5][0];

	    kv = TMath::Sqrt(1.0+xp*xp+yp*yp);

            al = xp/kv;
            be = yp/kv;
            ga = 1./kv;

            TTrack* fTrack = addTrack();
            fTrack->setTrack(x0,y0,t0, al,be,ga, ind0,ind1);

            evAl += al;
            evBe += be;
            evGa += ga;
        }
    }
    gEvent->setTracks(fTracks);
    gEvent->setRpcHitsCorr(fRpcHitCorr);
    // calculate incidence and mean vector (event direction)
    gEvent->setAngles(-100.,-100.,-100.);
    if(fTracks->GetEntriesFast()==0) {
        return 1;
}
    Float_t norm = 1./sqrt(evAl*evAl+evBe*evBe+evGa*evGa);

    evAl *= norm;
    evBe *= norm;
    evGa *= norm;

    gEvent->setAngles(evAl,evBe,evGa);

    // The formula to calculate the corrected incidence time is:
    // t0 - (evAl*x+evBe*y)/300.
    // And then the corrected hit times can be obtained as:
    Float_t trmin  = 1e5;
    Float_t trtime = 0.;
    for(Int_t i=0;i<fTracks->GetEntriesFast();i++) {
        TTrack* track = (TTrack*)fTracks->At(i);
        trtime = track->getTime()-(evAl*track->getX0()+evBe*track->getY0())/300.;
        if(trtime<trmin) trmin = trtime;
        track->setTime(trtime);
    }
    for(Int_t i=0;i<fTracks->GetEntriesFast();i++) {
        TTrack* track = (TTrack*)fTracks->At(i);
        track->setTime(track->getTime()-trmin);
    }

    totalNHitsCorr=0;

    Float_t thitmin2 = 1e6;
    for(Int_t i=0;i<fRpcHitHits->GetEntriesFast();i++) {
        TRpcHit* hit = (TRpcHit*)fRpcHitHits->At(i);
        if(!hit) continue;
        hit->getHit(trbnum, cell,col,row,x,y,z,time,charge);
        TRpcHit* fHit = addRpcHit();
        time -= (evAl*x+evBe*y)/300.;
        time -= thitmin;
        if(thitmin2>time) thitmin2 = time;
        //time -= trmin;
        fHit->setHit(trbnum, cell, col, row, x, y, z, time, charge);
    }

    for(Int_t i=0;i<fRpcHitCorr->GetEntriesFast();i++) {
        TRpcHit* fHit = (TRpcHit*)fRpcHitCorr->At(i);
        if(!fHit) continue;
        fHit->getHit(trbnum, cell,col,row,x,y,z,time,charge);
        time -= thitmin2;
        //time -= trmin;
        fHit->setHit(trbnum, cell, col, row, x, y, z, time, charge);
    }

    return 1;

}
TTrack* TTrackF::addTrack( ) {
    TClonesArray& tracks = *fTracks;
    TTrack *track = new (tracks[totalNHits++]) TTrack();
    //cout<<"ENTRIES in addhit "<<Hits->GetEntriesFast()<<endl;
    return track;
}
TRpcHit* TTrackF::addRpcHit( ) {
    TClonesArray& hits = *fRpcHitCorr;
    TRpcHit *rpchit = new (hits[totalNHitsCorr++]) TRpcHit();
    //cout<<"ENTRIES in addhit "<<Hits->GetEntriesFast()<<endl;
    return rpchit;
}

// AVector constructs the A matrix for a non linear model
TMatrixF TTrackF::AVector(TMatrixF SIn, Float_t x, Float_t y, Float_t t, Float_t z){
        TMatrixF A(6,1);
        Float_t X0=SIn[0][0], XP=SIn[1][0], Y0=SIn[2][0], YP=SIn[3][0], T0=SIn[4][0], S=SIn[5][0], k=TMath::Sqrt(1.0+XP*XP+YP*YP);
        A[0][0] = wx*x;
        A[1][0] = z*(wx*x*k*k+S*wt*XP*(t*k+S*(XP*XP+YP*YP)*z))/(k*k);
        A[2][0] = wy*y;
        A[3][0] = z*(wy*y*k*k+S*wt*YP*(t*k+S*(XP*XP+YP*YP)*z))/(k*k);
        A[4][0] = wt*(t*k+S*(XP*XP+YP*YP)*z)/k;
        A[5][0] = wt*z*(t*k+S*(XP*XP+YP*YP)*z);
        return A;
}

// KMatrix constructs the K matrix for a non linear model 
// // // The non-linear model has the form
// // // x = X0+X'*z
// // // y = Y0+Y'*z
// // // t = T0+k*S*z
TMatrixF TTrackF::KMatrix(TMatrixF SIn, Float_t z){
        TMatrixF K(6,6);
        Float_t X0=SIn[0][0], XP=SIn[1][0], Y0=SIn[2][0], YP=SIn[3][0], T0=SIn[4][0], S=SIn[5][0], k=TMath::Sqrt(1.0+XP*XP+YP*YP);
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

// Saeta2Planes calculates a saeta between 2 planes. A simple non-linear model used. 
 // Used to fill an non-zero input saeta
 // The output saeta has the form (X0,X',Y0,Y',T0,S)
TMatrixF TTrackF::Saeta2Planes(Float_t x1, Float_t y1, Float_t t1, Float_t z1, Float_t x2, Float_t y2, Float_t t2, Float_t z2){
	TMatrixF S2(6,1);
	Float_t Dz = z1-z2;
	S2[0][0] = (x2*z1-x1*z2)/Dz;
	S2[1][0] = (x1-x2)/Dz;
	S2[2][0] = (y2*z1-y1*z2)/Dz;
	S2[3][0] = (y1-y2)/Dz;
	S2[4][0] = (t2*z1-t1*z2)/Dz;
	S2[5][0] = (t1-t2)/Dz;
	return S2;
}                                                                                                                                          
