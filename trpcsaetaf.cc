#include "trpchitf.h"
#include "trpchit.h"
#include "ttmatrix.h"
#include "trpcsaeta.h"
#include "trpcsaetaf.h"
#include "stdlib.h"
#include "TMatrixF.h"
#include <math.h>

// Global constants: statistical weigth of the measured variables.
// // w=1/sigma^2
// // Distances in mm, time in ns 
// Float_t wx=11.0*sqrt(2), wy=12.0*sqrt(2), wt=0.3;
Float_t wx=116/sqrt(12), wy=111/sqrt(12), wt=0.3;

Float_t sx=116/sqrt(12), sy=111/sqrt(12), st=0.3; // Replace wx,wy,wz by sx,sy,sz
//Float_t wx = 1/pow(sx,2);
//Float_t wy = 1/pow(sy,2);
//Float_t wt = 1/pow(st,2);

Double_t vc = 299.792458; // mm/ns
Int_t nevent = 0;

//ClassImp(TTMatrix);
ClassImp(TRpcSaetaF);

TRpcSaetaF::TRpcSaetaF() {
    fRpcHitHits      = gEvent->getRpcHits();
    fRpcHitCorr      = new TClonesArray("TRpcHit",1000);
    fRpcSaeta2Planes = new TClonesArray("TRpcSaeta",1000);
    fRpcSaeta3Planes = new TClonesArray("TRpcSaeta",1000);
}

Int_t TRpcSaetaF::init() {
    fRpcHitHits      = gEvent->getRpcHits();
    if(!fRpcHitCorr)
        fRpcHitCorr      = new TClonesArray("TRpcHit",1000);
    if(!fRpcSaeta2Planes)
        fRpcSaeta2Planes = new TClonesArray("TRpcSaeta",1000);
    if(!fRpcSaeta3Planes)
        fRpcSaeta3Planes = new TClonesArray("TRpcSaeta",1000);
    return 1;
}

TRpcSaetaF::~TRpcSaetaF() {
    delete fRpcHitCorr;
    delete fRpcSaeta2Planes;
    delete fRpcSaeta3Planes;
}

Int_t TRpcSaetaF::execute() {
    // Clear the TClones array for later saving it!
    // clearAll();
    totalNHits2Planes = 0;
    totalNHits3Planes = 0;

    fRpcHitCorr 	 ->Clear("C");
    fRpcSaeta2Planes ->Clear("C");
    fRpcSaeta3Planes ->Clear("C");

    Float_t mod, mod2;
    Int_t maxIter  = 100; // Maximum number of iterations until convergence
    Float_t cutMod = 1e-5;

    Int_t   cell,  col,  row,  trbnum;
    Int_t   cell2, col2, row2, trbnum2;
    Int_t   cell3, col3, row3, trbnum3;
    Float_t x,  y,  z,  charge,  time;
    Float_t x2, y2, z2, charge2, time2;
    Float_t x3, y3, z3, charge3, time3;
    Float_t evAl = 0.;
    Float_t evBe = 0.;
    Float_t evGa = 0.;

    // RpcSaeta variables below
    Float_t xP, yP, x0, y0, t0, sl;
    Float_t xRec, yRec;
    Float_t kv, al, be, ga;
    Int_t san;
    Float_t chi2;
    Int_t ind, ind2, ind3;

    Float_t thitmin = 1e6;

    TMatrixF SInit(6,1), S2Planes(6,1), S3Planes(6,1), A(6,1), K(6,6);

    fRpcHitHits = gEvent->getRpcHits();
    if(!fRpcHitHits) return 1;

    Int_t nhits = fRpcHitHits->GetEntriesFast();
    if(nhits>10) return 1;
    //if(nhits>3) return 1;

    for(Int_t i=0;i<fRpcHitHits->GetEntriesFast();i++) {
        TRpcHit* hit = (TRpcHit*)fRpcHitHits->At(i);
        if(!hit) continue;
        hit->getHit(trbnum, cell,col,row,x,y,z,time,charge);
        //cout<<time<<" "<<trbnum<<" "<<cell<<endl;
        if(trbnum!=0) continue;
        if(thitmin>time) thitmin = time;
    }
    //-------------------------------------------------------------------------------------------------- Fit 2 planes - Old version

    Float_t sq;

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
            if(trbnum==0 && trbnum2==1 ) {
                if( fabs(time2-time-dist/300.)>0.9) continue;
                xP = (x2-x)/(z2-z);
                yP = (y2-y)/(z2-z);
                x0 = x;
                y0 = y;
                t0 = time - thitmin;
                ind  = i;
                ind2 = j;
            }
            if(trbnum==1 && trbnum2==0) {
                if( fabs(time-time2-dist/300.)>0.9) continue;
                xP = (x-x2)/(z-z2);
                yP = (y-y2)/(z-z2);
                x0 = x2;
                y0 = y2;
                t0 = time2 - thitmin;
                ind  = j;
                ind2 = i;
            }

            sq = sqrt(1.+xP*xP+yP*yP);
            al = xP/sq;
            be = yP/sq;
            ga = 1./sq;
            TRpcSaeta* fRpcSaeta2Planes = addRpcSaeta2Planes();
            fRpcSaeta2Planes->setRpcSaeta2Planes(x0,y0,t0,al,be,ga,ind,ind2);
            evAl += al;
            evBe += be;
            evGa += ga;
        }
    }

    gEvent->setRpcSaeta2Planes(fRpcSaeta2Planes);
    gEvent->setRpcHitsCorr(fRpcHitCorr);

    // calculate incidence and mean vector (event direction)

    gEvent->setAngles(-100.,-100.,-100.);
    gEvent->setAngles3(-100.,-100.,-100.);

    if(fRpcSaeta2Planes->GetEntriesFast()==0) {
        return 1;
    }

    Float_t normOld = 1./sqrt(evAl*evAl+evBe*evBe+evGa*evGa);

    evAl *= normOld;
    evBe *= normOld;
    evGa *= normOld;

    gEvent->setAngles(evAl,evBe,evGa);

    Float_t trmin  = 1e5;
    Float_t trtime = 0.;

    for(Int_t i=0;i<fRpcSaeta2Planes->GetEntriesFast();i++) {
        TRpcSaeta* RpcSaeta2 = (TRpcSaeta*)fRpcSaeta2Planes->At(i);
        trtime = RpcSaeta2->getTime()-(evAl*RpcSaeta2->getX0()+evBe*RpcSaeta2->getY0())/300.;
        if(trtime<trmin) trmin = trtime;
        RpcSaeta2->setTime(trtime);
    }
    for(Int_t i=0;i<fRpcSaeta2Planes->GetEntriesFast();i++) {
        TRpcSaeta* RpcSaeta2 = (TRpcSaeta*)fRpcSaeta2Planes->At(i);
        RpcSaeta2->setTime(RpcSaeta2->getTime()-trmin);
    }

    //---------------------------------------------------------------------------------------------------------- RpcSaeta 3 Planes
    evAl = 0.;
    evBe = 0.;
    evGa = 0.;

    //cout << "*********************************New Event with number of hits: " << nhits << endl;
    //ofstream myfile;
    //myfile.open ("example.txt");
    nevent += 1;
    //cout << nhits << " " << endl;
    for(Int_t i=0;i<nhits;i++) {
        TRpcHit* hit = (TRpcHit*)fRpcHitHits->At(i);
        if(!hit) continue;
        hit->getHit(trbnum, cell,col,row,x,y,z,time,charge);
        //cout << " nhits= " << nhits << " New event " << endl;
        // cout << " TRB= " << trbnum << endl;
        cout << nevent << " " << trbnum << " " << row << " " << col << " " << time << " " << charge << " " << endl;
        //if(trbnum!=0) continue;
        if(trbnum!=2) continue;
        //loop in strict order in order to prevent double counting.
        for(Int_t j=0;j<nhits;j++) {
            if(j==i) continue;
            TRpcHit* hit2 = (TRpcHit*)fRpcHitHits->At(j);
            if(!hit2) continue;
            hit2->getHit(trbnum2, cell2,col2,row2,x2,y2,z2,time2,charge2);
            //cout << " TRB1= " << trbnum << " TRB2= " << trbnum2 << endl;
            //if(trbnum2!=1) continue;
            if(trbnum2!=0) continue;
            for(Int_t k=0;k<nhits;k++) {
                if(k==j || k==i) continue;
                TRpcHit* hit3 = (TRpcHit*)fRpcHitHits->At(k);
                if(!hit3) continue;
                hit3->getHit(trbnum3, cell3,col3,row3,x3,y3,z3,time3,charge3);
                //if(trbnum3!=2) continue;
                //if(trbnum3!=2) continue;
                if(trbnum3==trbnum || trbnum3==trbnum2) continue;
                //cout << " TRB1= " << trbnum << " TRB2= " << trbnum2 << " TRB3= " << trbnum3 << endl;
                Float_t dist  = sqrt((x-x3)*(x-x3)+(y-y3)*(y-y3)+(z-z3)*(z-z3));
                //cout <<x<<" "<<y<<" "<<time<<" "<<charge<<" "<<x2<<" "<<y2<<" "<<time2<<" "<<charge2<<" "<<x3<<" "<<y3<<" "<<time3<<" "<<charge3<<" "<<fabs(time3-time-dist/vc) << endl;
                //cout<<"\nT1: "<<row<<" "<<col<<" "<<time<<" "<<charge<<" T2: "<<row2<<" "<<col2<<" "<<time2<<" "<<charge2<<" T3: "<<row3<<" "<<col3<<" "<<time3<<" "<<charge3<<" "<< endl;
                //if(trbnum==0 && trbnum2==1 && trbnum3==2) {
                if(trbnum==2 && trbnum2==0 && trbnum3==1) {
                    //cout << "PARECE QUE BIEN " << fabs(time-time3-dist/300.) << " "<< time3<< " "<< time<< endl;
                    //if(fabs(time3-time-dist/vc)>0.9) continue;
                    // Input Saeta2
                    // SInit = InputSaeta2Planes(x3,y3,time3,z3, x2,y2,time2,z2);
                    SInit = InputSaeta2Planes(x,y,time,z, x2,y2,time2,z2);
                    // Convergence loop
                    for (Int_t l=0;l<maxIter;l++){
                        //cout<<"Iter "<<l<<endl;
                        A = AVector(SInit, x,y,time,z) + AVector(SInit, x2,y2,time2,z2) + AVector(SInit, x3,y3,time3,z3);
                        K = KMatrix(SInit, z) + KMatrix(SInit, z2) + KMatrix(SInit, z3);
                        mod  = TMath::Sqrt(SInit.E2Norm());
                        S3Planes = K.Invert()*A;
                        SInit    = S3Planes;
                        mod2 = TMath::Sqrt(SInit.E2Norm());
                        //cout << mod  << " " << mod2 << " " << endl;
                        //cout << TMath::Abs(mod2-mod) << " " << endl;
                        if (TMath::Abs(mod2-mod)<cutMod){
                            //cout << "Saeta has converged: Iter = " << l << " "  << endl;
                            //cout<<x<<" "<<y<<" "<<time<<" "<<charge<<" "<<x2<<" "<<y2<<" "<<time2<<" "<<charge2<<" "<<x3<<" "<<y3<<" "<<time3<<" "<<charge3<<" "<<fabs(time3-time-dist/vc) << endl;
                            //cout<<"Saeta Converged :"<<row<<" "<<col<<" "<<time<<" "<<charge<<"  // "<<row2<<" "<<col2<<" "<<time2<<" "<<charge2<<" // "<<row3<<" "<<col3<<" "<<time3<<" "<<charge3<<" "<< endl;
                            break;
                        } // Saeta has converged
                        else{
                            //cout<<"Saeta didn't converge : "<<row<<" "<<col<<" "<<time<<" "<<charge<<"  //"<<row2<<" "<<col2<<" "<<time2<<" "<<charge2<<" //"<<row3<<" "<<col3<<" "<<time3<<" "<<charge3<<" "<< endl;

                        }
                    }
                    ind  = i;
                    ind2 = j;
                    ind3 = k;
                    x0 = S3Planes[0][0];
                    xP = S3Planes[1][0];
                    xRec = x0+xP*z;
                    y0 = S3Planes[2][0];
                    yP = S3Planes[3][0];
                    yRec = y0+yP*z;
                    t0 = S3Planes[4][0];
                    sl = S3Planes[5][0];

                    //cout << "Saeta : " <<x0<<" "<<xP<<" "<<y0<<" "<<yP<<" "<<t0<<" "<<sl<<" "<<endl;

                    //TMatrixF Schi = (S3Planes.Transpose())*K*S3Planes - (S3Planes.Transpose())*A - (S3Planes.Transpose())*A +

                    Double_t zpos[3] = {z,z2,z3};
                    Double_t xpos[3] = {x,x2,x3};
                    Double_t ypos[3] = {y,y2,y3};
                    Double_t tpos[3] = {time,time2,time3};
                    Double_t chi2full = 0;

                    kv = TMath::Sqrt(1.0+xP*xP+yP*yP);


                    Float_t dist1 = sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2));
                    Float_t dist2 = sqrt((x-x3)*(x-x3)+(y-y3)*(y-y3)+(z-z3)*(z-z3));

                    Float_t distZigZag = dist1 + dist2;

                    al = xP/kv;
                    be = yP/kv;
                    ga = 1./kv;

                    Float_t distSaeta = 1./ga*(z3-z2);
                    Float_t tSaeta = distSaeta*sl;
                    for(int ich=0;ich<3;ich++) {
                        Double_t xteo = x0+xP*zpos[ich];
                        Double_t yteo = y0+yP*zpos[ich];
                        Double_t tteo = t0+1./ga*zpos[ich]*sl;

                        //cout<<" plano "<<ich<<" chx "<<(xteo-xpos[ich])*(xteo-xpos[ich])/wx/wx<<" chy "<<(yteo-ypos[ich])*(yteo-ypos[ich])/wy/wy<<" cht "<<(tteo-tpos[ich])*(tteo-tpos[ich])/wt/wt<<endl;
                        //	cout<<" x "<<xteo<<" "<<xpos[ich]<<" y "<<yteo<<" "<<ypos[ich]<<" t "<<tteo<<" "<<tpos[ich]<<endl;
                        chi2full += (xteo-xpos[ich])*(xteo-xpos[ich])/wx/wx;
                        chi2full += (yteo-ypos[ich])*(yteo-ypos[ich])/wy/wy;
                        chi2full += (tteo-tpos[ich])*(tteo-tpos[ich])/wt/wt;
                    }

                    san  = 2; //temporary
                    //chi2 = (x-xRec)*(x-xRec)/wx/wx+(y-yRec)*(y-yRec)/wy/wy; //temporary
                    chi2 = chi2full;//temporary

//					if(nhits==3 && -1./sl/300.>1.2 && chi2<20 &&fRpcSaeta2Planes->GetEntriesFast()==1 ) {
                    //cout<<" dist and velocity "<<distZigZag<<" "<<distSaeta<<" vs "<< -1./sl/300.<<" vz "<<distZigZag/(time3-time2)/300.<<" t2 "<<time2 <<" t3 "<< time3 <<" "<< charge2  <<" "<<charge3 <<endl;
//} 	

                    TRpcSaeta* fRpcSaeta3Planes = addRpcSaeta3Planes();

                    fRpcSaeta3Planes->setRpcSaeta3Planes(x0,xP,y0,yP,z, t0,1./sl/vc, al,be,ga, san,ind,ind2,ind3, chi2);
                    //fRpcSaeta3Planes->setRpcSaeta3Planes(x0,xP,y0,yP,z, time-time3, distZigZag/vc , charge,charge2,charge3, san,ind,ind2,ind3, chi2);
                    //fRpcSaeta3Planes->setRpcSaeta3Planes(x0,xP,y0,yP,z, t0, (z-z2)/(time-time2)/300., al,be,ga, san,ind,ind2,ind3, chi2);

                    evAl += al;
                    evBe += be;
                    evGa += ga;
                }
            }
        }
    }

    gEvent->setRpcSaeta3Planes(fRpcSaeta3Planes); // FIXME: gEvent things
    gEvent->setRpcHitsCorr(fRpcHitCorr);
    // calculate incidence and mean vector (event direction)
    gEvent->setAngles3(-100.,-100.,-100.);


    //cout<<"entries "<< fRpcSaeta3Planes->GetEntriesFast()<<endl;

    if(fRpcSaeta3Planes->GetEntriesFast()==0) {
        return 1;
    }

    Float_t norm = 1./sqrt(evAl*evAl+evBe*evBe+evGa*evGa);

    evAl *= norm;
    evBe *= norm;
    evGa *= norm;

    gEvent->setAngles3(evAl,evBe,evGa);

    // The formula to calculate the corrected incidence time is:
    // t0 - (evAl*x+evBe*y)/300.
    // And then the corrected hit times can be obtained as:

    trmin  = 1e5;
    trtime = 0.;

    //cout<<evAl<<" "<<evBe<<" "<<evGa<<endl;

    /*
    for(Int_t i=0;i<fRpcSaeta3Planes->GetEntriesFast();i++) {
        TRpcSaeta* RpcSaeta3 = (TRpcSaeta*)fRpcSaeta3Planes->At(i);
        //cout << "Break" << RpcSaeta3 << endl;
        trtime = RpcSaeta3->getTime()-(evAl*RpcSaeta3->getX0()+evBe*RpcSaeta3->getY0())/300.;
        if(trtime<trmin) trmin = trtime;
        RpcSaeta3->setTime(trtime);
    }
    for(Int_t i=0;i<fRpcSaeta3Planes->GetEntriesFast();i++) {
        TRpcSaeta* RpcSaeta3 = (TRpcSaeta*)fRpcSaeta3Planes->At(i);
        RpcSaeta3->setTime(RpcSaeta3->getTime()-trmin);
    }
    */

    totalNHitsCorr=0;
    Float_t thitmin2 = 1e6;

    for(Int_t i=0;i<fRpcHitHits->GetEntriesFast();i++) {
        TRpcHit* hit = (TRpcHit*)fRpcHitHits->At(i);
        if(!hit) continue;
        hit->getHit(trbnum, cell,col,row,x,y,z,time,charge);
        TRpcHit* fHit = addRpcHit();
        time -= (evAl*x+evBe*y)/vc;
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


TRpcSaeta* TRpcSaetaF::addRpcSaeta2Planes( ) {
    TClonesArray& RpcSaeta2Planes = *fRpcSaeta2Planes;
    //cout << "Total Hits 2 Planes: " << totalNHits2Planes << endl;
    TRpcSaeta *RpcSaeta = new (RpcSaeta2Planes[totalNHits2Planes++]) TRpcSaeta();
    //cout<<"ENTRIES in addhit "<<Hits->GetEntriesFast()<<endl;
    return RpcSaeta;
}

TRpcSaeta* TRpcSaetaF::addRpcSaeta3Planes( ) {
    TClonesArray& RpcSaeta3Planes = *fRpcSaeta3Planes;
    //cout << "Total Hits 3 Planes: " << totalNHits3Planes << endl;
    TRpcSaeta *RpcSaeta = new (RpcSaeta3Planes[totalNHits3Planes++]) TRpcSaeta();
    //cout<<"ENTRIES in addhit "<<Hits->GetEntriesFast()<<endl;
    return RpcSaeta;
}
//---------------------------------------------------------------------------------------------------------
/*	for(Int_t i=0;i<nhits-1;i++) {
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
		    cout<<nhits<<" "<<trbnum<<" "<<trbnum2<<" "<<time<<" "<<time2<<" "<<dist/300.<<" "<<x<<" "<<x2<<" "<<endl;
		    if(trbnum==0 && trbnum2==1) {
			if( fabs(time2-time-dist/300.)>0.9) continue;
			// Input Saeta
			SInit = InputSaeta2Planes(x,y,time,z, x2,y2,time2,z2);
			// Convergence loop
			for (Int_t k=0;k<maxIter;k++){
			    A = AVector(SInit, x,y,time,z) + AVector(SInit, x2,y2,time2,z2);
			    K = KMatrix(SInit, z) + KMatrix(SInit, z2);
			    mod  = TMath::Sqrt(SInit.E2Norm());
			    S2Planes = K.Invert()*A;
			    SInit    = S2Planes;
			    mod2 = TMath::Sqrt(SInit.E2Norm());
			    if (TMath::Abs(mod2-mod)<cutMod){break;} // Saeta has converged	
			}
			ind = i;
			ind2 = j;
		    }
		    if(trbnum==1 && trbnum2==0) {
		    	if( fabs(time-time2-dist/300.)>0.9) continue;
		    	// Input Saeta
		    	SInit = InputSaeta2Planes(x2,y2,time2,z2, x,y,time,z);
		    	// Convergence loop
		    	for (Int_t k=0;k<maxIter;k++){
				A = AVector(SInit, x2,y2,time2,z2) + AVector(SInit, x,y,time,z);
				K = KMatrix(SInit, z2) + KMatrix(SInit, z);
				mod  = TMath::Sqrt(SInit.E2Norm());
				S2Planes = K.Invert()*A;
				SInit    = S2Planes;
				mod2 = TMath::Sqrt(SInit.E2Norm());
			if (TMath::Abs(mod-mod2)<cutMod){break;} // Saeta has converged	
		    	}
			ind  = j;
			ind2 = i;
		    }

		    x0 = S2Planes[0][0];
		    xP = S2Planes[1][0];
		    xRec = x0+xP*z;
		    y0 = S2Planes[2][0];
		    yP = S2Planes[3][0];
		    yRec = y0+yP*z;
		    t0 = S2Planes[4][0];
		    sl = S2Planes[5][0];

		    kv = TMath::Sqrt(1.0+xP*xP+yP*yP);

		    al = xP/kv;
		    be = yP/kv;
		    ga = 1./kv;

		    san  = 2; //temporary
		    chi2 = (x0-xRec)*(x0-xRec)*wx+(y0-yRec)*(y0-yRec)*wy; //temporary

		    TRpcSaeta* fRpcSaeta2Planes = addRpcSaeta2Planes();
		    fRpcSaeta2Planes->setRpcSaeta2Planes(x0,xP,y0,yP,z, t0,sl, al,be,ga, san,ind,ind2, chi2);

		    evAl += al;
		    evBe += be;
		    evGa += ga;

		}
	}

	gEvent->setRpcSaeta2Planes(fRpcSaeta2Planes); // FIXME: gEvent things
	gEvent->setRpcHitsCorr(fRpcHitCorr);
	// calculate incidence and mean vector (event direction)
	gEvent->setAngles(-100.,-100.,-100.);

	if(fRpcSaeta2Planes->GetEntriesFast()==0) {
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

	for(Int_t i=0;i<fRpcSaeta2Planes->GetEntriesFast();i++) {
		TRpcSaeta* RpcSaeta2 = (TRpcSaeta*)fRpcSaeta2Planes->At(i);
		trtime = RpcSaeta2->getTime()-(evAl*RpcSaeta2->getX0()+evBe*RpcSaeta2->getY0())/300.;
		if(trtime<trmin) trmin = trtime;
		RpcSaeta2->setTime(trtime);
	}
	for(Int_t i=0;i<fRpcSaeta2Planes->GetEntriesFast();i++) {
		TRpcSaeta* RpcSaeta2 = (TRpcSaeta*)fRpcSaeta2Planes->At(i);
		RpcSaeta2->setTime(RpcSaeta2->getTime()-trmin);
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

TRpcSaeta* TRpcSaetaF::addRpcSaeta2Planes( ) {
    TClonesArray& RpcSaeta2Planes = *fRpcSaeta2Planes;
    cout << " " << totalNHits << endl;
    TRpcSaeta *RpcSaeta = new (RpcSaeta2Planes[totalNHits++]) TRpcSaeta();
    //cout<<"ENTRIES in addhit "<<Hits->GetEntriesFast()<<endl;
    return RpcSaeta;
}*/
TRpcHit* TRpcSaetaF::addRpcHit( ) {
    TClonesArray& hits = *fRpcHitCorr;
    TRpcHit *rpchit = new (hits[totalNHitsCorr++]) TRpcHit();
    //cout<<"ENTRIES in addhit "<<Hits->GetEntriesFast()<<endl;
    return rpchit;
}

// AVector constructs the A matrix for a non linear model
TMatrixF TRpcSaetaF::AVector(TMatrixF SIn, Float_t x, Float_t y, Float_t t, Float_t z){
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
TMatrixF TRpcSaetaF::KMatrix(TMatrixF SIn, Float_t z){
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
TMatrixF TRpcSaetaF::InputSaeta2Planes(Float_t x1, Float_t y1, Float_t t1, Float_t z1, Float_t x2, Float_t y2, Float_t t2, Float_t z2){
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









                                                                                                                          
