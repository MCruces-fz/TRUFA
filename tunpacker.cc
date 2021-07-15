
#include "tunpacker.h"
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "thldevent.h"
#include "tevent.h"
#include "TString.h"
#include "thit.h"
#include "trpclookuptable.h"
#include "trpcrawf.h"
#include "trpchitf.h"
#include "trpcraw.h"
#include "trpchit.h"
#include "trpcsaetaf.h"
#include "trpcsaeta.h"
#include "ttmatrix.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include <fstream>
#include "TTimeStamp.h"
#include "TMath.h"
#include "TRandom3.h"
//#include <ofstream>


using namespace std;

ClassImp(Unpacker)

//______________________________________________________________________________
Unpacker::Unpacker() 
{
   EventNr=0;
   pEvent=0;
   pRootFile=0;
   subEvtId=0;
   fpga_code=0;
   refCh = -1;

   //fileLookupPar      = 0;
   //fileHitFinderPar   = 0;
   //fileTrackFinderPar = 0;
   /*
    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		hq[i*10*12+j*12+k] = new TH1D(Form("h_q_d%i_r%i_c%i",i,j,k),
				       Form("h_q_d%i_r%i_c%i;#Q [ns];Entries",i,j,k),1000,0,200);
		hdta[i*10*12+j*12+k] = new TH1D(Form("h_dta_d%i_r%i_c%i",i,j,k),
				       Form("h_dta_d%i_r%i_c%i;#Q [ns];Entries",i,j,k),250,-40,40);
	    }
	}
    }
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    hdt[i*12*10*12+j*10*12+k*12+l] = new TH1D(Form("h_dt_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       Form("h_dt_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       250,-40,40);

		}
	    }
	}
    }
    hdt[10*12*10*12] = new TH1D("h_dt_all","h_dt_all;Dt [ns];entries",500,-80,80);

    h1D[0] = new TH1D("h_rate_H1_1particle","h_rate_H1_1particle",6*24*6,90,96);
    h1D[1] = new TH1D("h_rate_H1_2particle","h_rate_H1_2particle",6*24*6,90,96);
    h1D[2] = new TH1D("h_rate_H2_1particle","h_rate_H2_1particle",6*24*6,90,96);
    h1D[3] = new TH1D("h_rate_H2_2particle","h_rate_H2_2particle",6*24*6,90,96);
    h1D[4] = new TH1D("h_phi_M1","h_phi_M1",360,-180,180);
    h1D[5] = new TH1D("h_dt_M1","h_dt_M1",1000,-10,10);
    h1D[6] = new TH1D("h_timeDistH1","h_timeDistH1",1000,0,100);
    h1D[7] = new TH1D("h_timeDistH2","h_timeDistH2",1000,0,100);
    h1D[8] = new TH1D("h_Flux1","h_Flux1",367*24*6,0,367);
    h1D[9] = new TH1D("h_DAQactive","h_DAQactive",367*24*6,0,367);

    h2D[0] = new TH2D("h_mult_H1_H2","h_mult_H1_H2",100,0,100,100,0,100);
    h2D[1] = new TH2D("h_time_H1_mult_H1","h_time_H1_mult_H1",100,0,100,100,0,100);
    h2D[2] = new TH2D("h_time_H2_mult_H2","h_time_H2_mult_H2",100,0,100,100,0,100);
    h2D[3] = new TH2D("h_dt_M1_Q0","h_dt_M1_Q0",1000,-10,10,1000,0,10000);
    h2D[4] = new TH2D("h_dt_M1_Q1","h_dt_M1_Q1",1000,-10,10,1000,0,10000);
    h2D[5] = new TH2D("h_dt_d","h_dt_d",1000,-100,100,150,1000,2500);
    h2D[6] = new TH2D("h_dt_q","h_dt_q",200,-10,10,400,0,4000);
    h2D[7] = new TH2D("h_dt_logq","h_dt_logq",200,-10,10,400,0.01,5);

    */



}

Unpacker::Unpacker(Int_t j, const char* name,Int_t nEvt ,const char* subEvtId,Int_t refChannel, const char* fpga_code)
{
    // Initialize the options
    refCh = -1;
    refCh = refChannel;
    this->subEvtId= HexStrToInt(subEvtId);
    setInputFile(name);
    EventNr=0;
    pEvent= new HldEvent(inputFile.c_str(), HexStrToInt(subEvtId), fpga_code);
    pRootFile=0;
    //  this->fpga_code=fpga_code;
    if(nEvt>0)
    {
	eventLoop(nEvt);
    }
}

Unpacker::Unpacker(const char* dir, const char* name, const char* odir,const char* ofile, Int_t nEvt)
{
    TFile* fOut = new TFile(ofile,"RECREATE");
    TH1D *hq[3*10*12];
    TH1D *hdta[3*10*12];
    TH1D* hdt[10*12*10*12];
    TH1D* hdt2[10*12*10*12];
    //TH1D**** hq[][][];
    //TH1D***** hdt[][][][];

    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		hq[i*10*12+j*12+k] = new TH1D(Form("h_q_d%i_r%i_c%i",i,j,k),
				       Form("h_q_d%i_r%i_c%i;#Q [ns];Entries",i,j,k),1000,-20,20);
		hdta[i*10*12+j*12+k] = new TH1D(Form("h_dta_d%i_r%i_c%i",i,j,k),
				       Form("h_dta_d%i_r%i_c%i;#Q [ns];Entries",i,j,k),1000,-20,20);
	    }
	}
    }
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    hdt[i*12*10*12+j*10*12+k*12+l] = new TH1D(Form("h_dt_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       Form("h_dt_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       250,-40,40);
		    hdt2[i*12*10*12+j*10*12+k*12+l] = new TH1D(Form("h_dt2_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       Form("h_dt2_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       250,-40,40);

		}
	    }
	}
    }



    refCh = -1;
    refCh = 31;
    this->subEvtId= 0;
    outputFile = Form("%s%s.root",odir,name);
    cout<<"ofile "<<outputFile<<endl;
    setInputFile(dir,name);
    EventNr=0;
    pEvent= new HldEvent(inputFile.c_str(), 999, " ", -100000,100000);
    pEvent->setQuietMode(true);
    pEvent->setFullSetup(true);
    pEvent->setVHR(false);
    pRootFile=0;
    if(nEvt>0)
    {
	eventLoopFillCal(nEvt,0,hq,hdt,hdt2);
    }

    fOut->cd();

    TF1* fg = new TF1("fg","gaus",-40,40);
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    hdt[i*12*10*12+j*10*12+k*12+l]->Fit(fg,"NQR");
		    hdta[i*12+j]->Fill(fg->GetParameter(1));
		    hdta[12*10+k*12+l]->Fill(fg->GetParameter(1));

		}
	    }
	}
    }

    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {

		hq[i*10*12+j*12+k]->Write();
                hdta[i*10*12+j*12+k]->Write();

	    }
	}
    }
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    hdt[i*12*10*12+j*10*12+k*12+l]->Write();
		    hdt2[i*12*10*12+j*10*12+k*12+l]->Write();

		}
	    }
	}
    }
    fOut->Close();




}
Unpacker::Unpacker(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n)
{
    TFile* fOut = new TFile(ofile,"RECREATE");
    TH1D *hq[3*10*12];
    TH1D *hdta[3*10*12];
    TH1D* hdt[10*12*10*12];
    TH1D* hdt2[10*12*10*12];
    //TH1D**** hq[][][];
    //TH1D***** hdt[][][][];

    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		hq[i*10*12+j*12+k] = new TH1D(Form("h_q_d%i_r%i_c%i",i,j,k),
				       Form("h_q_d%i_r%i_c%i;#Q [ns];Entries",i,j,k),1000,0,200);
		hdta[i*10*12+j*12+k] = new TH1D(Form("h_dta_d%i_r%i_c%i",i,j,k),
				       Form("h_dta_d%i_r%i_c%i;#Q [ns];Entries",i,j,k),250,-40,40);
	    }
	}
    }
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    hdt[i*12*10*12+j*10*12+k*12+l] = new TH1D(Form("h_dt_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       Form("h_dt_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       250,-40,40);
		    hdt2[i*12*10*12+j*10*12+k*12+l] = new TH1D(Form("h_dt2_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       Form("h_dt2_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       250,-40,40);

		}
	    }
	}
    }




    refCh = -1;
    refCh = 31;
    this->subEvtId= 0;
    //outputFile = Form("%s%s.root",odir,name);
    cout<<"ofile "<<outputFile<<endl;
    //TFile* file = TFile::Open("list");
    std::ifstream file(list.Data());
    TString fName;
    Int_t ntot=0;
    while(!file.eof()) {
	file>>fName;
	if(!fName||!fName.Contains("hld"))break;
	setInputFile(dir,fName.Data());
	EventNr=0;
	pEvent= new HldEvent(inputFile.c_str(), 999, " ", -100000,100000);
	pEvent->setQuietMode(true);
	pEvent->setFullSetup(true);
	pEvent->setVHR(false);
	pRootFile=0;
        
	if(nEvt>0 && ntot <nEvt )
	{
	    eventLoopFillCal(nEvt,0,hq,hdt,hdt2);
	    ntot++;
	}
        //if(ntot>10)break;
    }
    fOut->cd();

    TF1* fg = new TF1("fg","gaus",-40,40);



    Float_t val[3][10][12];
    Float_t valq[3][10][12];
    Float_t valE[3][10][12];


    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		val[i][j][k] = -999.;
		valE[i][j][k] = 2.;
	    }
	}
    }

    val[0][5][5] = 0.;

    Float_t cons[10][12][10][12];
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
                    cons[i][j][k][l] = -999;
		    //if(sqrt((i-k)*(i-k)+(j-l)*(j-l))>4.)continue;
		    if(hdt[i*12*10*12+j*10*12+k*12+l]->GetEntries()<30 || hdt[i*10*12+j*12+k]->Integral(1,250)<20) continue;
		    Int_t max = hdt[i*12*10*12+j*10*12+k*12+l]->GetMaximumBin();
		    Float_t cent = hdt[i*12*10*12+j*10*12+k*12+l]->GetBinCenter(max);
		    //Float_t cent = hdt[i*12*10*12+j*10*12+k*12+l]->GetMean();
		    //if(cent==0)continue;
		    //if(cent==0 || cent<-30 || cent>30 || hdt[i*12*10*12+j*10*12+k*12+l]->GetRMS()>14.) continue;
		    if(cent==0 || cent<-30 || cent>30 ) continue;
		    fg->SetParLimits(2,0.01,1.5);
		    fg->SetParameters(hdt[i*12*10*12+j*10*12+k*12+l]->GetMaximum(),cent,0.7);
		    hdt[i*12*10*12+j*10*12+k*12+l]->Fit(fg,"NQR","",cent-2,cent+2);
		    hdt[i*12*10*12+j*10*12+k*12+l]->Fit(fg,"QR","",cent-2,cent+2);

		    if(TMath::Abs(fg->GetParameter(1))<20 && fg->GetParameter(2)<1.5)
			cons[i][j][k][l] = fg->GetParameter(1);

		}
	    }
	}
    }

   
    for(Int_t i=5;i<6;i++) {
	for(Int_t j=5;j<6;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    val[1][k][l] = cons[i][j][k][l];
		    cout<<"val 1 "<<k<<" "<<l<<" "<<val[1][k][l]<<endl;

		}
	    }
	}
    }

    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=5;k<6;k++) {
		for(Int_t l=5;l<6;l++) {
                    if(i==5 && j==5) continue;
		    val[0][i][j] = (val[1][k][l]-cons[i][j][k][l]);
		    cout<<"val 0 "<<i<<" "<<j<<" "<<val[0][i][j]<<endl;
		}
	    }
	}
    }










    /*
    for(Int_t i=0;i<1;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		valE[i][j][k] = 2.;
		fg->SetParameters(10,0,1.);
		if(hdta[i*10*12+j*12+k]->GetEntries()<5 || hdta[i*10*12+j*12+k]->Integral(1,250)<2) {
		    hdta[i*10*12+j*12+k]->Reset();
		    continue;
		}
		Int_t max    = hdta[i*12*10+j*12+k]->GetMaximumBin();
		Float_t cent = hdta[i*12*10+j*12+k]->GetBinCenter(max);
		if(cent==0 || cent<-30 || cent>30) continue;
		//fg->SetParameters(hdta[i*12*10+j*12+k]->GetMaximum(),cent,1.);
		//hdta[i*10*12+j*12+k]->Fit(fg,"NQR","",cent-4,cent+4);
		//hdta[i*10*12+j*12+k]->Fit(fg,"NQR","",cent-4,cent+4);
		//ta[i*10*12+j*12+k]->Fit(fg,"NQR",);
		val[i][j][k] = cent;
		//valE[i][j][k] = fabs(fg->GetParameter(2));
		hdta[i*10*12+j*12+k]->Reset();
	    }
	}
    }
    */

    TH2D* h_test = new TH2D("h_test","h_test",10,0,10,12,0,12);
    //for(Int_t i=0;i<10;i++) {
    //    for(Int_t j=0;j<12;j++) {
    //        h_test->SetBinContent(i,j,(val[1][i][j]+cons[0][0][i][j]));
    //    }
    //}

    for(Int_t num = 0;num<50;num++) {

	cout<<"NUM "<<num<<endl;


	for(Int_t i=0;i<10;i++) {
	    for(Int_t j=0;j<12;j++) {
		hdta[i*12+j]->Reset();

		//h_test->Reset();

		//for(Int_t k=0;k<10;k++) {
		//    for(Int_t l=0;l<12;l++) {
			// fg->SetParameters(10,0,0.5);
			// if(sqrt((i-k)*(i-k)+(j-l)*(j-l))> 9. )continue;
		//	hdta[i*12+j]->Fill(cons[i][j][k][l]+val[1][k][l], 0.05*( (i-k)*(i-k) + (j-l)*(j-l) ) );
		//	h_test->SetBinContent(k,l,cons[i][j][k][l]+val[1][k][l]);
			// hdta[12*10+k*12+l]->Fill(-cons[i][j][k][l]+val[0][i][j]);
		//    }
		//}
		//Int_t max    = hdta[i*12+j]->GetMaximumBin();
		//Float_t cent = hdta[i*12+j]->GetBinCenter(max);
		// if(cent<-30 || cent>30) continue;
		//fg->SetParameters(hdta[i*12+j]->GetMaximum(),cent,1.);
		//hdta[i*12+j]->Fit(fg,"NQR","",cent-2,cent+2);

		//val[0][i][j] = val[0][i][j]+0.3*(val[0][i][j]+fg->GetParameter(1));
		//val[1][i][j] = -cons[i][j][i][j]+val[0][i][j];

		for(Int_t k=0;k<10;k++) {
		    for(Int_t l=0;l<12;l++) {
			//if(h_test->GetBinContent(k,l)-val[0][i][j]<5.) {
			if(i!=k &&j!=l) {
			    if(TMath::Abs(cons[i][j][k][l])<20 && TMath::Abs(cons[k][l][k][l])<20 ) {
				val[1][k][l] = cons[i][j][k][l]+val[0][i][j];// val[1][k][l]+(-h_test->GetBinContent(k,l)+val[0][i][j]);
				val[0][k][l] = -cons[k][l][k][l]+val[1][k][l];
			    }
			}
			//}
		    }
		}
	    }
	}

        /*
	for(Int_t k=0;k<10;k++) {
	    for(Int_t l=0;l<12;l++) {
		hdta[12*10+k*12+l]->Reset();
		h_test->Reset();
		for(Int_t i=0;i<10;i++) {
		    for(Int_t j=0;j<12;j++) {
			// fg->SetParameters(10,0,0.5);
			// if(sqrt((i-k)*(i-k)+(j-l)*(j-l))> 9. )continue;
			// hdta[i*12+j]->Fill(cons[i][j][k][l]+val[1][k][l]);
			hdta[12*10+k*12+l]->Fill(-cons[i][j][k][l]+val[0][i][j], 0.05*( (i-k)*(i-k) + (j-l)*(j-l) ));
		    }
		}
		Int_t max    = hdta[12*10+k*12+l]->GetMaximumBin();
		Float_t cent = hdta[12*10+k*12+l]->GetBinCenter(max);
		// if(cent<-30 || cent>30) continue;
		fg->SetParameters(hdta[12*10+k*12+l]->GetMaximum(),cent,1.);
		hdta[12*10+k*12+l]->Fit(fg,"NQR","",cent-2,cent+2);
		val[1][k][l] = val[1][k][l]+0.3*(-val[1][k][l]+fg->GetParameter(1));
		val[0][k][l] = cons[k][l][k][l]+val[1][k][l];
                */
                /*
		for(Int_t i=0;i<10;i++) {
		    for(Int_t j=0;j<12;j++) {

			if(h_test->GetBinContent(k,l)-val[0][i][j]<5.) {
			val[0][i][j] = val[0][i][j]+0.5*(-val[0][i][j]+fg->GetParameter(1));
			val[1][i][j] = -cons[i][j][i][j]+val[0][i][j];


		    }
		}
                */

	    //}
	//}

        /*
	//for(Int_t i=0;i<3;i++) {
        Int_t i=0;
	if(num%2==0){
            i=0;
	} else i=1;

	    for(Int_t j=0;j<10;j++) {
		for(Int_t k=0;k<12;k++) {
		    fg->SetParameters(10,0,1.);
		    if(hdta[i*10*12+j*12+k]->GetEntries()<4 || hdta[i*10*12+j*12+k]->Integral(1,250)<2) {
			if(num<(4))
			hdta[i*10*12+j*12+k]->Reset();
			continue;
		    }
		    Int_t max    = hdta[i*12*10+j*12+k]->GetMaximumBin();
		    Float_t cent = hdta[i*12*10+j*12+k]->GetBinCenter(max);
		    // if(cent<-30 || cent>30) continue;
		    fg->SetParameters(hdta[i*12*10+j*12+k]->GetMaximum(),cent,1.);
		    hdta[i*10*12+j*12+k]->Fit(fg,"NQR","",cent-3,cent+3);
		    hdta[i*10*12+j*12+k]->Fit(fg,"QR","",cent-2,cent+2);
		    //hdta[i*10*12+j*12+k]->Fit(fg,"NQR",);
		    //val[i][j][k] = fg->GetParameter(1);
		    //if(j==1&&k==2) cout<<"new - old val "<<i<<" "<<j<<" "<<k<<" "<<cent<<" "<<val[i][j][k]<<endl;
		    val[i][j][k] = cent;
		    //valE[i][j][k] = fabs(fg->GetParameter(2));
                    if(num<(9))
			hdta[i*10*12+j*12+k]->Reset();
		}
	    }
	    //}
	    */
    }





    /*
    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		val[i][j][k] = 0.;
		fg->SetParameters(10,0,0.5);
		if(hdta[i*10*12+j*12+k]->GetEntries()<5) {
		    hdta[i*10*12+j*12+k]->Reset();
		    continue;
		}
		Float_t cent = hdta[i*12*10+j*12+k]->GetMean();
		if(cent==0 || cent<-30 || cent>30) continue;
		//fg->SetParameters(10,cent,0.5);
                fg->SetParameters(hdta[i*12*10+j*12+k]->GetMaximum(),cent,0.5);
		hdta[i*10*12+j*12+k]->Fit(fg,"QR+","",cent-4,cent+4);
		//hdta[i*10*12+j*12+k]->Fit(fg,"NQR","",cent-4,cent+4);
		//hdta[i*10*12+j*12+k]->Fit(fg,"NQR");
		val[i][j][k] = fg->GetParameter(1);
                //hdta[i*10*12+j*12+k]->Reset();
	    }
	}
    }
    */
    h_test->Reset();
    for(Int_t k=0;k<10;k++) {
	for(Int_t l=0;l<12;l++) {
	    h_test->SetBinContent(k,l,-cons[9][10][k][l]+val[1][k][l]);
	    // hdta[12*10+k*12+l]->Fill(-cons[i][j][k][l]+val[0][i][j]);
	}
    }


    h_test->Write();


    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		valq[i][j][k] = 0;
		Int_t   maxb = hq[i*10*12+j*12+k]->GetMaximumBin();
		Float_t max  = hq[i*10*12+j*12+k]->GetBinCenter(maxb);
		fg->SetParameters(hq[i*10*12+j*12+k]->GetMaximum(),max,0.3);
                fg->SetParLimits(1,max-0.5,max+0.5);
//                fg->SetParLimit(1,max-0.5,max+0.5);
		hq[i*10*12+j*12+k]->Fit(fg,"NQR","",max-1.,max+1.);
		hq[i*10*12+j*12+k]->Fit(fg,"QR","",fg->GetParameter(1)-fg->GetParameter(2)*3.,fg->GetParameter(1)+fg->GetParameter(2)*3.);
		hq[i*10*12+j*12+k]->Fit(fg,"QR","",fg->GetParameter(1)-fg->GetParameter(2)*3.,fg->GetParameter(1)+fg->GetParameter(2)*3.);
                valq[i][j][k] = fg->GetParameter(1)+fg->GetParameter(2)*2.;
		hq[i*10*12+j*12+k]->Write();
		hdta[i*10*12+j*12+k]->Write();
	    }
	}
    }
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    hdt[i*12*10*12+j*10*12+k*12+l]->Write();

		}
	    }
	}
    }
    fOut->Close();

    ofstream fo("parhitfNewEmptyTimeChargePedestal.txt");
    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		//fo<<i<<" "<<j<<" "<<k<<" "<<valq[i][j][k]<<" "<<val[i][j][k]<<endl;
		fo<<i<<" "<<j<<" "<<k<<" "<<valq[i][j][k]<<" "<<0.0<<endl;
	    }
	}
    }
    fo.close();




}

void Unpacker::fillCalibration(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n)
{
    // n = 0, standard mode, just calculate the pedestals and exchande in the parameter file which is previously declared.
    // n = 1, special mode:  create pedestals and set time offsets to 0.

    TFile fOut(Form("%s%s",odir,ofile),"RECREATE");

    TH1D *hq[3*10*12];
    TH1D *hdta[3*10*12];
    TH1D* hdt[10*12*10*12+1];
    TH1D* hdt2[10*12*10*12+1];

    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		hq[i*10*12+j*12+k] = new TH1D(Form("h_q_d%i_r%i_c%i",i,j,k),
				       Form("h_q_d%i_r%i_c%i;#Q [ns];Entries",i,j,k),1000,0,200);
		hdta[i*10*12+j*12+k] = new TH1D(Form("h_dta_d%i_r%i_c%i",i,j,k),
				       Form("h_dta_d%i_r%i_c%i;#Q [ns];Entries",i,j,k),250,-20,20);
	    }
	}
    }
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    hdt[i*12*10*12+j*10*12+k*12+l] = new TH1D(Form("h_dt_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       Form("h_dt_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       250,-20,20);
		    hdt2[i*12*10*12+j*10*12+k*12+l] = new TH1D(Form("h_dt2_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       Form("h_dt2_r1_%i_c1_%i_r2_%i_c2_%i",i,j,k,l),
					       250,-40,40);

		}
	    }
	}
    }
    hdt[10*12*10*12] = new TH1D("h_dt_all","h_dt_all;Dt [ns];entries",500,-80,80);
    hdt2[10*12*10*12] = new TH1D("h_dt2_all","h_dt2_all;Dt [ns];entries",500,-80,80);


    //TH1D**** hq[][][];
    //TH1D***** hdt[][][][];
    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		hq[i*10*12+j*12+k]->Reset();
				  
		hdta[i*10*12+j*12+k]->Reset();
	    }
	}
    }
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    hdt[i*12*10*12+j*10*12+k*12+l]->Reset();
		    hdt2[i*12*10*12+j*10*12+k*12+l]->Reset();

		}
	    }
	}
    }
    hdt[10*12*10*12]->Reset();
    hdt2[10*12*10*12]->Reset();


    refCh = -1;
    refCh = 31;
    this->subEvtId= 0;
    //outputFile = Form("%s%s.root",odir,name);
    cout<<"ofile "<<outputFile<<endl;
    cout<<"inputlist "<<list.Data()<<endl;
    //TFile* file = TFile::Open("list");
    std::ifstream file(list.Data());
    TString fName;
    Int_t ntot=0;
    while(!file.eof()) {
	file>>fName;
	cout<<"FILENAME: "<<fName.Data()<<" Filedir "<<dir<<endl;
	if(!fName||!fName.Contains("hld"))break;
	setInputFile(dir,fName.Data());
        cout<<" File set "<<inputFile.c_str()<<endl;
	EventNr=0;
        pEvent = new HldEvent(inputFile.c_str(), 999, " ", -100000,100000);
        /*
	if(!pEvent)pEvent= new HldEvent(inputFile.c_str(), 999, " ", -100000,100000);
	else {

	    pEvent->setFile(inputFile.c_str());
	    pEvent->setSubEvtId(999);
	    pEvent->init();
	}
        */
	pEvent->setQuietMode(true);
	pEvent->setFullSetup(true);
	pEvent->setVHR(false);
        pEvent->setDebugFlag(0);
        pEvent->setDebugFlag1(0);
	pRootFile=0;
	cout<<"Starting event loop "<<endl;
        cout<<"Pointer pEvent "<<pEvent<<endl;
	if(nEvt>0 && ntot <nEvt )
	{
	    eventLoopFillCal(nEvt,0,hq,hdt,hdt2);
	    ntot++;
	}
	delete pEvent;
        pEvent = NULL;
        //if(ntot>10)break;
    }
    //file.close();
    fOut.cd();


    TF1* fg = new TF1("fg","gaus",-1000,1000);
    TF1* fgt = new TF1("fgt","gaus",-1000,1000);
    Float_t valq[3][10][12];
    hdt[10*12*10*12]->Write();
    hdt2[10*12*10*12]->Write();

    TF1* fg2 = new TF1("fg2","gaus",-80,80);
    Int_t maxb = hdt[10*12*10*12]->GetMaximumBin();
    fg2->SetParameters(hdt[10*12*10*12]->GetMaximum(),hdt[10*12*10*12]->GetBinCenter(maxb),0.5);
    hdt[10*12*10*12]->Fit(fg2,"R","",hdt[10*12*10*12]->GetBinCenter(maxb)-1.,hdt[10*12*10*12]->GetBinCenter(maxb)+1.);
    //cout<<"Sync times "<<fg2->GetParameter(1)<<endl;
    Float_t stime  = fg2->GetParameter(1);

    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		valq[i][j][k] = 0;
		Int_t   maxb = hq[i*10*12+j*12+k]->GetMaximumBin();
		Float_t max  = hq[i*10*12+j*12+k]->GetBinCenter(maxb);
		fg->SetParameters(hq[i*10*12+j*12+k]->GetMaximum(),max,0.3);
                fg->SetParLimits(1,max-0.5,max+0.5);
		hq[i*10*12+j*12+k]->Fit(fg,"NQR","",max-1.,max+1.);
		hq[i*10*12+j*12+k]->Fit(fg,"QR","",fg->GetParameter(1)-fg->GetParameter(2)*3.,fg->GetParameter(1)+fg->GetParameter(2)*3.);
		hq[i*10*12+j*12+k]->Fit(fg,"QR","",fg->GetParameter(1)-fg->GetParameter(2)*3.,fg->GetParameter(1)+fg->GetParameter(2)*3.);
		valq[i][j][k] = fg->GetParameter(1)+fg->GetParameter(2)*2.;
		if(n==0||n==1)
		    hq[i*10*12+j*12+k]->Write();
		// hdta[i*10*12+j*12+k]->Write();
	    }
	}
    }
    Double_t dt[10][12][10][12];//={0.};
    Bool_t   dtb[10][12][10][12];//={0};
    Double_t dt2[10][12][10][12];//={0.};
    Bool_t   dtb2[10][12][10][12];//={0};

    //TF1* fg = new TF1("fg","gaus",-1000,1000);
    Double_t dtc[3][10][12];//={0.};


    for(Int_t p=0;p<3;p++) {
	for(Int_t i=0;i<10;i++) {
	    for(Int_t j=0;j<12;j++) {
                dtc[p][i][j]=0.;
	    }
	}
    }




    if(n==0) {

	TRandom3 rnd;
	rnd.SetSeed(0);
	for(Int_t i=0;i<10;i++) {
	    for(Int_t j=0;j<12;j++) {
		for(Int_t k=0;k<10;k++) {
		    for(Int_t l=0;l<12;l++) {
                        dt[i][j][k][l]=0.;
                        dtb[i][j][k][l]=0;
                        dt2[i][j][k][l]=0.;
                        dtb2[i][j][k][l]=0;


			if(n==0) {
			    hdt[i*12*10*12+j*10*12+k*12+l]->Write();
			    hdt2[i*12*10*12+j*10*12+k*12+l]->Write();
			}

                        if(hdt[i*12*10*12+j*10*12+k*12+l] -> GetEntries()>20.) {
			    fgt->SetParameters(hdt[i*12*10*12+j*10*12+k*12+l]->GetMean(),hdt[i*12*10*12+j*10*12+k*12+l]->GetMean(),0.5);
			    hdt[i*12*10*12+j*10*12+k*12+l] ->Fit(fgt,"WWNQR","");
			    hdt[i*12*10*12+j*10*12+k*12+l] ->Fit(fgt,"WWNQR","",
								 fgt->GetParameter(1)-fgt->GetParameter(2)*4.,
								 fgt->GetParameter(1)+fgt->GetParameter(2)*4.);
			    dt[i][j][k][l] = fgt->GetParameter(1);
			    dtb[i][j][k][l] = 1;
			    //cout<<" i "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<dt[i][j][k][l]<<endl;

			}

			if(hdt2[i*12*10*12+j*10*12+k*12+l] -> GetEntries()>20.) {
                            fgt->SetParameters(hdt2[i*12*10*12+j*10*12+k*12+l]->GetMean(),hdt2[i*12*10*12+j*10*12+k*12+l]->GetMean(),0.5);
			    //fgt->SetParameters(1.,0.,0.1);
			    hdt2[i*12*10*12+j*10*12+k*12+l] ->Fit(fgt,"WWNQR","");
			    hdt2[i*12*10*12+j*10*12+k*12+l] ->Fit(fgt,"WWNQR","",
								  fgt->GetParameter(1)-fgt->GetParameter(2)*4.,
								  fgt->GetParameter(1)+fgt->GetParameter(2)*4.);
			    dt2[i][j][k][l] = fgt->GetParameter(1);
			    //cout<<" i2 "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<dt2[i][j][k][l]<<endl;
			    dtb2[i][j][k][l] = 1;
			}
		    }
		}
	    }
	}

	TH1F* h_dt[3][10][12];
	TH1F* h_dt2[3][10][12];
        TH1F* h_dt3[3][10][12];
        TH1F* h_dt4[3][10][12];
        TH1F* h_dt5[3][10][12];


       // cout<<"Reach this point"<<endl;

	for(Int_t p=0;p<3;p++) {
	    for(Int_t i=0;i<10;i++) {
		for(Int_t j=0;j<12;j++) {
		    //cout<<p<<" "<<i<<" "<<j<<endl;
		    h_dt[p][i][j] = new TH1F(Form("h_dt_p%i_r%i_c%i",p,i,j),Form("h_dt_p%i_r%i_c%i",p,i,j),1000,-50,50);
		}
	    }
	}

	for(Int_t i=0;i<10;i++) {
	    for(Int_t j=0;j<12;j++) {
		for(Int_t k=0;k<10;k++) {
		    for(Int_t l=0;l<12;l++) {
			if(dtb[i][j][k][l]) {
			    h_dt[0][i][j] -> Fill(0.5*dt[i][j][k][l]);
			    h_dt[1][k][l] -> Fill(-0.5*dt[i][j][k][l]);
			}
			if(dtb2[i][j][k][l]) {
			    h_dt[2][i][j] -> Fill(0.5*dt2[i][j][k][l]);
			}
		    }
		}
	    }
	}

	for(Int_t p=0;p<3;p++) {
	    for(Int_t i=0;i<10;i++) {
                for(Int_t j=0;j<12;j++) {
                    h_dt[p][i][j] ->ClearUnderflowAndOverflow();
		    Int_t maxbin = h_dt[p][i][j] -> GetMaximumBin();
		    Double_t binc = h_dt[p][i][j] -> GetBinCenter(maxbin);
		    dtc[p][i][j] = binc;
		    //cout<<"first val "<<p<<" "<<i<<" "<<j<<" "<<dtc[p][i][j]<<endl;
		}
	    }
	}

	
	for(Int_t p=0;p<3;p++) {
	    for(Int_t i=0;i<10;i++) {
		for(Int_t j=0;j<12;j++) {
		    //cout<<p<<" "<<i<<" "<<j<<endl;
		    h_dt2[p][i][j] = new TH1F(Form("h_dt2_p%i_r%i_c%i",p,i,j),Form("h_dt2_p%i_r%i_c%i",p,i,j),1000,-20,20);
		}
            }
	}

	for(Int_t it = 0;it<20;it++) {
	    //cout<<"Reached the 20 iterations "<<it<<endl;
	    for(Int_t p=0;p<3;p++) {
		for(Int_t i=0;i<10;i++) {
		    for(Int_t j=0;j<12;j++) {
			h_dt2[p][i][j] -> Reset();
		    }
		}
	    }
	    for(Int_t i=0;i<10;i++) {
		for(Int_t j=0;j<12;j++) {
		    for(Int_t k=0;k<10;k++) {
			for(Int_t l=0;l<12;l++) {
			    if(dtb[i][j][k][l]) {
				Double_t rand = rnd.Rndm();
				h_dt2[0][i][j] -> Fill((1.-rand)*(dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]));
				h_dt2[1][k][l] -> Fill(-rand    *(dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]));
			    }
			    if(dtb2[i][j][k][l]) {
				h_dt2[2][i][j] -> Fill((1.)*(dt2[i][j][k][l]+dtc[1][k][l]-dtc[2][i][j]));
			    }
			}
		    }
		}
	    }
	    for(Int_t p=0;p<3;p++) {
		for(Int_t i=0;i<10;i++) {
		    for(Int_t j=0;j<12;j++) {
			Int_t maxbin = h_dt2[p][i][j] -> GetMaximumBin();
			Double_t binc = h_dt2[p][i][j] -> GetBinCenter(maxbin);
			dtc[p][i][j] += binc;
			//if(it==19) {
			 //   cout<<"20 val "<<p<<" "<<i<<" "<<j<<" "<<dtc[p][i][j]<<endl;
			//}
		    }
		}
	    }
	}

	
	for(Int_t p=0;p<3;p++) {
	    for(Int_t i=0;i<10;i++) {
		for(Int_t j=0;j<12;j++) {
		    //cout<<p<<" "<<i<<" "<<j<<endl;
		    h_dt3[p][i][j] = new TH1F(Form("h_dt3_p%i_r%i_c%i",p,i,j),Form("h_dt3_p%i_r%i_c%i",p,i,j),1000,-10,10);
		}
	    }
        }

        for(Int_t it = 0;it<100;it++) {

            for(Int_t p=0;p<3;p++) {
                for(Int_t i=0;i<10;i++) {
                    for(Int_t j=0;j<12;j++) {
                        h_dt3[p][i][j] -> Reset();
                    }
                }
            }
            for(Int_t i=0;i<10;i++) {
                for(Int_t j=0;j<12;j++) {
                    for(Int_t k=0;k<10;k++) {
                        for(Int_t l=0;l<12;l++) {
                            if(dtb[i][j][k][l]) {
                                Double_t rand = rnd.Rndm();
                                h_dt3[0][i][j] -> Fill((1.-rand)*(dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]));
                                h_dt3[1][k][l] -> Fill(-rand    *(dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]));
                            }
                            if(dtb2[i][j][k][l]) {
                                h_dt3[2][i][j] -> Fill((1.)*(dt2[i][j][k][l]+dtc[1][k][l]-dtc[2][i][j]));
                                //cout<<dt[i][j][k][l]<<" "<<dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]<<" "<<dtc[0][i][j]<<" "<<dtc[1][k][l]<<endl;
                            }
                        }
                    }
                }
            }
            for(Int_t p=0;p<3;p++) {
                for(Int_t i=0;i<10;i++) {
                    for(Int_t j=0;j<12;j++) {
                        Int_t maxbin = h_dt3[p][i][j] -> GetMaximumBin();
                        Double_t binc = h_dt3[p][i][j] -> GetBinCenter(maxbin);
                        fg->SetParameter(1,binc);

                        /*
                         fg->SetParameters(10,binc,0.5);
                         if( (it<10 && it%3==0) || it>14 ) {
                         h_dt3[p][i][j] -> Fit(fg,"WWNQR","",-5.5,5.5);
                         h_dt3[p][i][j] -> Fit(fg,"WWNQR","",fg->GetParameter(1)-fg->GetParameter(2)*3.,fg->GetParameter(1)+fg->GetParameter(2)*3.);
                         }  else if(it<20 && it%2==0){
                         h_dt3[p][i][j] -> Fit(fg,"WWNQR","",-5.5,5.5);
                         h_dt3[p][i][j] -> Fit(fg,"NQR","",fg->GetParameter(1)-fg->GetParameter(2)*3.,fg->GetParameter(1)+fg->GetParameter(2)*3.);
                         }
                         else if(it>40 && it <60) {
                         h_dt3[p][i][j] -> Fit(fg,"WWNQR","",-5.5,5.5);
                         h_dt3[p][i][j] -> Fit(fg,"WWNQR","",fg->GetParameter(1)-fg->GetParameter(2)*3.,fg->GetParameter(1)+fg->GetParameter(2)*3.);
                         }
                         else {
                         Int_t maxbin = h_dt3[p][i][j] -> GetMaximumBin();
                         Double_t binc = h_dt3[p][i][j] -> GetBinCenter(maxbin);
                         fg->SetParameter(1,binc);
                         }
                         */
                        // Add some random!
                        //
                        Double_t rand = rnd.Rndm()-0.5;
                        rand*=(0.2-(it*0.002)) + h_dt3[p][i][j] -> GetBinWidth(maxbin)*(rnd.Rndm()-0.5);;
                        dtc[p][i][j] += fg->GetParameter(1)+rand;
                        //cout<<p<<" "<<i<<" "<<j<<" "<<dtc[p][i][j]<<endl;
                    }
                }
            }
        }



        for(Int_t p=0;p<3;p++) {
            for(Int_t i=0;i<10;i++) {
                for(Int_t j=0;j<12;j++) {
                    //cout<<p<<" "<<i<<" "<<j<<endl;
                    h_dt4[p][i][j] = new TH1F(Form("h_dt4_p%i_r%i_c%i",p,i,j),Form("h_dt4_p%i_r%i_c%i",p,i,j),1000,-5,5);
                }
            }
        }
        for(Int_t it = 0;it<100;it++) {
            for(Int_t p=0;p<3;p++) {
                for(Int_t i=0;i<10;i++) {
                    for(Int_t j=0;j<12;j++) {
                        h_dt4[p][i][j] -> Reset();
                    }
                }
            }
            for(Int_t i=0;i<10;i++) {
                for(Int_t j=0;j<12;j++) {
                    for(Int_t k=0;k<10;k++) {
                        for(Int_t l=0;l<12;l++) {
                            if(dtb[i][j][k][l]) {
                                Double_t rand = rnd.Rndm();
                                h_dt4[0][i][j] -> Fill((1.-rand)*(dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]));
                                h_dt4[1][k][l] -> Fill(-rand    *(dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]));
                            }
                            if(dtb2[i][j][k][l]) {
                                h_dt4[2][i][j] -> Fill((1.)*(dt2[i][j][k][l]+dtc[1][k][l]-dtc[2][i][j]));
                                //cout<<dt[i][j][k][l]<<" "<<dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]<<" "<<dtc[0][i][j]<<" "<<dtc[1][k][l]<<endl;
                            }
                        }
                    }
                }
            }
            for(Int_t p=0;p<3;p++) {
                for(Int_t i=0;i<10;i++) {
                    for(Int_t j=0;j<12;j++) {
                        Int_t maxbin = h_dt4[p][i][j] -> GetMaximumBin();
                        Double_t binc = h_dt4[p][i][j] -> GetBinCenter(maxbin);
                        fg->SetParameter(1,binc);
                        Double_t rand = rnd.Rndm()-0.5;
                        rand*=(0.15-(it*0.0015)) + h_dt4[p][i][j] -> GetBinWidth(maxbin)*(rnd.Rndm()-0.5);
                        dtc[p][i][j] += fg->GetParameter(1)+rand;
                        //cout<<p<<" "<<i<<" "<<j<<" "<<dtc[p][i][j]<<endl;

                    }
                }
            }
        }

        for(Int_t p=0;p<3;p++) {
            for(Int_t i=0;i<10;i++) {
                for(Int_t j=0;j<12;j++) {
                    //cout<<p<<" "<<i<<" "<<j<<endl;
                    h_dt5[p][i][j] = new TH1F(Form("h_dt5_p%i_r%i_c%i",p,i,j),Form("h_dt5_p%i_r%i_c%i",p,i,j),1000,-2,2);
                }
            }
        }
        for(Int_t it = 0;it<100;it++) {
            for(Int_t p=0;p<3;p++) {
                for(Int_t i=0;i<10;i++) {
                    for(Int_t j=0;j<12;j++) {
                        h_dt5[p][i][j] -> Reset();
                    }
                }
            }
            for(Int_t i=0;i<10;i++) {
                for(Int_t j=0;j<12;j++) {
                    for(Int_t k=0;k<10;k++) {
                        for(Int_t l=0;l<12;l++) {
                            if(dtb[i][j][k][l]) {
                                Double_t rand = rnd.Rndm();
                                h_dt5[0][i][j] -> Fill((1.-rand)*(dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]));
                                h_dt5[1][k][l] -> Fill(-rand    *(dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]));
                            }
                            if(dtb2[i][j][k][l]) {
                                h_dt5[2][i][j] -> Fill((1.)*(dt2[i][j][k][l]+dtc[1][k][l]-dtc[2][i][j]));
                                //cout<<dt[i][j][k][l]<<" "<<dt[i][j][k][l]+dtc[1][k][l]-dtc[0][i][j]<<" "<<dtc[0][i][j]<<" "<<dtc[1][k][l]<<endl;
                            }
                        }
                    }
                }
            }
            for(Int_t p=0;p<3;p++) {
                for(Int_t i=0;i<10;i++) {
                    for(Int_t j=0;j<12;j++) {
                        Int_t maxbin = h_dt5[p][i][j] -> GetMaximumBin();
                        Double_t binc = h_dt5[p][i][j] -> GetBinCenter(maxbin);
                        fg->SetParameter(1,binc);
                        if(p==2 && it>50) {
                            fg->SetParameter(1, h_dt5[p][i][j] -> GetMean());
                        }
                        Double_t rand = rnd.Rndm()-0.5;
                        rand*=(0.1-(it*0.001)) + h_dt5[p][i][j] -> GetBinWidth(maxbin)*(rnd.Rndm()-0.5);
                        dtc[p][i][j] += fg->GetParameter(1)+rand;
                        //cout<<p<<" "<<i<<" "<<j<<" "<<dtc[p][i][j]<<endl;
                    }
                }
            }
        }
    }

    for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		delete hq[i*10*12+j*12+k];
		delete hdta[i*10*12+j*12+k];
	    }
	}
    }
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    delete hdt[i*12*10*12+j*10*12+k*12+l];

		}
	    }
	}
    }
    delete hdt[10*12*10*12];






    fOut.Close();


    delete fg;
    delete fg2;














    

    //ofstream fo("parhitfNewEmptyTimeChargePedestal.txt");
    if(n==0) {
	ofstream fo(fileHitFinderParOut);
	TRpcCalPar* calPar = new TRpcCalPar(fileHitFinderPar);
	for(Int_t i=0;i<3;i++) {
	    for(Int_t j=0;j<10;j++) {
		for(Int_t k=0;k<12;k++) {
		    //fo<<i<<" "<<j<<" "<<k<<" "<<valq[i][j][k]<<" "<<val[i][j][k]<<endl;

		    Float_t calpar = calPar->getTimeCal(i, k, j) - dtc[i][j][k];
		    //cout<<i<<" "<<j<<" "<<k<<" "<<valq[i][j][k]<<" "<<calpar<<" "<<calPar->getTimeCal(i, k, j)<<" "<<dtc[i][j][k]<<endl;
		    //if(i==0)
		    fo<<i<<" "<<j<<" "<<k<<" "<<valq[i][j][k]<<" "<<calpar<<endl;
		    //else
		    //    fo<<i<<" "<<j<<" "<<k<<" "<<valq[i][j][k]<<" "<<calPar->getTimeCal(i, k, j) <<endl;
		}
	    }
	}
	fo.close();
        delete calPar;
    }

    if(n==1) {
	ofstream fo(fileHitFinderParOut);
	TRpcCalPar* calPar = new TRpcCalPar(fileHitFinderPar);
	for(Int_t i=0;i<3;i++) {
	    for(Int_t j=0;j<10;j++) {
		for(Int_t k=0;k<12;k++) {
		    //fo<<i<<" "<<j<<" "<<k<<" "<<valq[i][j][k]<<" "<<val[i][j][k]<<endl;
		    fo<<i<<" "<<j<<" "<<k<<" "<<valq[i][j][k]<<" "<< 0.0 <<endl;
		}
	    }
	}
	fo.close();
        delete calPar;
    }




}
Int_t Unpacker::syncCheck(const char* dir, TString file, Int_t nEvt,Int_t n)
{
    // n = 0, standard mode, just calculate the pedestals and exchande in the parameter file which is previously declared.
    // n = 1, special mode:  create pedestals and set time offsets to 0.

    Int_t sync = 1;
    refCh = -1;
    refCh = 31;
    this->subEvtId= 0;
    Int_t ntot=0;
    setInputFile(dir,file.Data());
    EventNr=0;
    pEvent = new HldEvent(inputFile.c_str(), 999, " ", -100000,100000);
    pEvent->setQuietMode(true);
    pEvent->setFullSetup(true);
    pEvent->setVHR(false);
    pEvent->setDebugFlag(0);
    pEvent->setDebugFlag1(0);
    pRootFile=0;
    if(nEvt>0 && ntot <nEvt )
    {
	sync = eventLoopSyncCheck(nEvt,0);
	ntot++;
    }
    delete pEvent;
    pEvent = NULL;
    return sync;

}

Unpacker::Unpacker(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n,Int_t n2)
{
    TFile* fOut = new TFile(ofile,"RECREATE");

    TH1D* h1D[8];
    TH2D* h2D[5];
    TH3D* h3D[6];

    h1D[0] = new TH1D("h_rate_H1_1particle","h_rate_H1_1particle",6*24*6,90,96);
    h1D[1] = new TH1D("h_rate_H1_2particle","h_rate_H1_2particle",6*24*6,90,96);
    h1D[2] = new TH1D("h_rate_H2_1particle","h_rate_H2_1particle",6*24*6,90,96);
    h1D[3] = new TH1D("h_rate_H2_2particle","h_rate_H2_2particle",6*24*6,90,96);
    h1D[4] = new TH1D("h_phi_M1","h_phi_M1",360,-180,180);
    h1D[5] = new TH1D("h_dt_M1","h_dt_M1",1000,-10,10);
    h1D[6] = new TH1D("h_timeDistH1","h_timeDistH1",1000,0,100);
    h1D[7] = new TH1D("h_timeDistH2","h_timeDistH2",1000,0,100);

    h2D[0] = new TH2D("h_mult_H1_H2","h_mult_H1_H2",100,0,100,100,0,100);
    h2D[1] = new TH2D("h_time_H1_mult_H1","h_time_H1_mult_H1",100,0,100,100,0,100);
    h2D[2] = new TH2D("h_time_H2_mult_H2","h_time_H2_mult_H2",100,0,100,100,0,100);
    h2D[3] = new TH2D("h_dt_M1_Q0","h_dt_M1_Q0",1000,-10,10,1000,0,10000);
    h2D[4] = new TH2D("h_dt_M1_Q1","h_dt_M1_Q1",1000,-10,10,1000,0,10000);
    h2D[5] = new TH2D("h_dt_d","h_dt_d",1000,-100,100,150,1000,2500);





    refCh = -1;
    refCh = 31;
    this->subEvtId= 0;
    //outputFile = Form("%s%s.root",odir,name);
    cout<<"ofile "<<outputFile<<endl;
    //TFile* file = TFile::Open("list");
    std::ifstream file(list.Data());
    TString fName;
    while(!file.eof()) {
	file>>fName;
	if(!fName||!fName.Contains("hld"))break;
	setInputFile(dir,fName.Data());
	EventNr=0;
	pEvent = new HldEvent(inputFile.c_str(), 999, " ", -100000,100000);
	pEvent->setQuietMode(true);
	pEvent->setFullSetup(true);
	pEvent->setVHR(false);
	pRootFile=0;
	if(nEvt>0)
	{
	    eventLoopFillCal(nEvt,0,h1D,h2D,h3D);
	}
        //delete pEvent;
    }
    fOut->cd();
    h2D[0]->Write();
    h2D[3]->Write();
    h2D[4]->Write();
    h2D[5]->Write();
    h1D[5]->Write();





    delete h1D[0];
    delete h1D[1];
    delete h1D[2];
    delete h1D[3];
    delete h1D[4];
    delete h1D[5];
    delete h1D[6];
    delete h1D[7];

    delete h2D[0];
    delete h2D[1];
    delete h2D[2];
    delete h2D[3];
    delete h2D[4];
    delete h2D[5];


    fOut->Close();
    //delete fOut;


}
void Unpacker::fillHistograms(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n,Int_t n2)
{
//    TFile* fOut = new TFile(ofile,"RECREATE");
    TFile fOut(Form("%s%s",odir,ofile),"RECREATE");

    TH1D* h1D[21];
    TH2D* h2D[23];
    TH3D* h3D[6];


    h1D[0] = new TH1D("h_rate_H1_1particle","h_rate_H1_1particle",6*24*6,90,96);
    h1D[1] = new TH1D("h_rate_H1_2particle","h_rate_H1_2particle",6*24*6,90,96);
    h1D[2] = new TH1D("h_rate_H2_1particle","h_rate_H2_1particle",6*24*6,90,96);
    h1D[3] = new TH1D("h_rate_H2_2particle","h_rate_H2_2particle",6*24*6,90,96);
    h1D[4] = new TH1D("h_phi_M1","h_phi_M1",360,-180,180);
    h1D[5] = new TH1D("h_dt_M1","h_dt_M1",1000,-10,10);
    h1D[6] = new TH1D("h_timeDistH1","h_timeDistH1",1000,0,100);
    h1D[7] = new TH1D("h_timeDistH2","h_timeDistH2",1000,0,100);
    h1D[8] = new TH1D("h_Flux1","h_Flux1",367*24*6,0,367);
    h1D[9] = new TH1D("h_DAQactive","h_DAQactive",367*24*6,0,367);


    h1D[10] = new TH1D("h_Q_M1","h_Q1_M1",367*24*6,0,367);
    h1D[11] = new TH1D("h_Q_N","h_Q_N",367*24*6,0,367);
    h1D[12] = new TH1D("h_Flux_M1","h_FluxM1",367*24*6,0,367);
    h1D[13] = new TH1D("h_Flux_M2","h_FluxM2",367*24*6,0,367);
    h1D[14] = new TH1D("h_Flux_M3","h_FluxM3",367*24*6,0,367);


    h1D[15] = new TH1D("h_Flux_s30_M1","h_Flux_s30_M1",367*24*6,0,367);
    h1D[16] = new TH1D("h_Flux_l30_M1","h_Flux_l30_M1",367*24*6,0,367);
    h1D[17] = new TH1D("h_Flux_s30_M2","h_Flux_s30_M2",367*24*6,0,367);
    h1D[18] = new TH1D("h_Flux_l30_M2","h_Flux_l30_M2",367*24*6,0,367);
    h1D[19] = new TH1D("h_Flux_s30_M3","h_Flux_s30_M3",367*24*6,0,367);
    h1D[20] = new TH1D("h_Flux_l30_M3","h_Flux_l30_M3",367*24*6,0,367);





    h2D[0] = new TH2D("h_mult_H1_H2","h_mult_H1_H2",100,0,100,100,0,100);
    h2D[1] = new TH2D("h_time_H1_mult_H1","h_time_H1_mult_H1",100,0,100,100,0,100);
    h2D[2] = new TH2D("h_time_H2_mult_H2","h_time_H2_mult_H2",100,0,100,100,0,100);
    h2D[3] = new TH2D("h_dt_M1_Q0","h_dt_M1_Q0",1000,-10,10,1000,0,10000);
    h2D[4] = new TH2D("h_dt_M1_Q1","h_dt_M1_Q1",1000,-10,10,1000,0,10000);
    h2D[5] = new TH2D("h_dt_d","h_dt_d",1000,-100,100,150,1000,2500);
    h2D[6] = new TH2D("h_dt_q","h_dt_q",200,-10,10,400,0,4000);
    h2D[7] = new TH2D("h_dt_logq","h_dt_logq",200,-10,10,400,0.01,5);

    for(Int_t i=0;i<5;i++) {
	h2D[8+i] = new TH2D(Form("h_rate2d_M1_theta%i",i),Form("h_rate2d_M1_theta%i;DOY ;#phi [Deg]",i),367*24*6,0,367,8,-180,180);
    }
    for(Int_t i=0;i<5;i++) {
	h2D[5+8+i] = new TH2D(Form("h_rate2d_M2_theta%i",i),Form("h_rate2d_M2_theta%i;DOY ;#phi [Deg]",i),367*24*6,0,367,8,-180,180);
    }
    for(Int_t i=0;i<5;i++) {
	h2D[5+5+8+i] = new TH2D(Form("h_rate2d_M3_theta%i",i),Form("h_rate2d_M3_theta%i;DOY ;#phi [Deg]",i),367*24*6,0,367,8,-180,180);
    }



    h1D[0] ->Reset();
    h1D[1] ->Reset();
    h1D[2] ->Reset();
    h1D[3] ->Reset();
    h1D[4] ->Reset();
    h1D[5] ->Reset();
    h1D[6] ->Reset();
    h1D[7] ->Reset();
    h1D[8] ->Reset();
    h1D[9] ->Reset();

    h2D[0] ->Reset();
    h2D[1] ->Reset();
    h2D[2] ->Reset();
    h2D[3] ->Reset();
    h2D[4] ->Reset();
    h2D[5] ->Reset();
    h2D[6] ->Reset();
    h2D[7] ->Reset();






    refCh = -1;
    refCh = 31;
    this->subEvtId= 0;
    //outputFile = Form("%s%s.root",odir,name);
    cout<<"ofile "<<outputFile<<endl;
    //TFile* file = TFile::Open("list");
    std::ifstream file(list.Data());
    TString fName;
    while(!file.eof()) {
	file>>fName;
	if(!fName||!fName.Contains("hld"))break;
	setInputFile(dir,fName.Data());
	EventNr=0;
	pEvent = new HldEvent(inputFile.c_str(), 999, " ", -100000,100000);
	pEvent->setQuietMode(true);
	pEvent->setFullSetup(true);
	pEvent->setVHR(false);
	pEvent->setDebugFlag(0);
	pEvent->setDebugFlag1(0);
	pRootFile=0;
	if(nEvt>0)
	{
	    eventLoopFillCal(nEvt,0,h1D,h2D,h3D);
	}
	delete pEvent;
        pEvent = NULL;
    }
    file.close();
    fOut.cd();

    /*
    h2D[0]->Write();
    h2D[3]->Write();
    h2D[4]->Write();
    h2D[5]->Write();
    h1D[5]->Write();
    h1D[8]->Write();
    h1D[9]->Write();

    h2D[6] -> Write();
    h2D[7] -> Write();
    */

    for(Int_t i=0;i<21;i++)
	h1D[i]->Write();
    for(Int_t i=0;i<23;i++)
	h2D[i]->Write();

    for(Int_t i=0;i<21;i++)
	delete h1D[i];
    for(Int_t i=0;i<23;i++)
	delete h2D[i];



    fOut.Close();


    //fOut.Close();

    //delete fOut;

}

Unpacker::Unpacker(const char* dir,const char* name, const char* odir, Int_t nEvt, const char* subEvtId,
		   Int_t refChannel, const char* fpga_code, Int_t min, Int_t max, Int_t quietMode,
		   Int_t fullSetup, Int_t VHR)
{
  refCh = -1;
  refCh = refChannel;
  this->subEvtId= HexStrToInt(subEvtId);	
  outputFile = Form("%s%s.root",odir,name);
  cout<<"ofile "<<outputFile<<endl;
  setInputFile(dir,name);
  EventNr=0;
  pEvent= new HldEvent(inputFile.c_str(), HexStrToInt(subEvtId), fpga_code, min, max); 
  
  if (quietMode == 1) {
    pEvent->setQuietMode(true);
  }
  else if(quietMode == 0) {
    pEvent->setQuietMode(false);
  }
  
  if (fullSetup == 1) {
    pEvent->setFullSetup(true);
  }
  else if(fullSetup == 0) {
    pEvent->setFullSetup(false);
  }
  
  if(VHR == 1) {
   pEvent->setVHR(true); 
  }
  else if(VHR == 0) {
   pEvent->setVHR(false); 
  }
  
  pRootFile=0;
  if(nEvt>0)
  {
      eventLoop(nEvt,0);
  }
}
Unpacker::Unpacker(const char* dir,const char* name, const char* odir, Int_t nEvt, TString luptab, TString calpar)
//const char* subEvtId,
//		   Int_t refChannel, const char* fpga_code, Int_t min, Int_t max, Int_t quietMode,
//		   Int_t fullSetup, Int_t VHR)
{
  refCh = -1;
  refCh = 31;
  this->subEvtId= HexStrToInt("999");
  outputFile = Form("%s%s.root",odir,name);
  cout<<"ofile "<<outputFile<<endl;
  setInputFile(dir,name);
  EventNr=0;
  pEvent= new HldEvent(inputFile.c_str(), HexStrToInt("999"), "999", -1000000, 1000000);

  setFileLookupPar(luptab);
  setFileHitFinderPar(calpar);

  pEvent->setQuietMode(true);
  //pEvent->setQuietMode(false);
  pEvent->setFullSetup(true);
  //pEvent->setFullSetup(false);
  //if(VHR == 1) {
  // pEvent->setVHR(true);
  //}
  //else if(VHR == 0) {
  pEvent->setVHR(false);
  //}

  pEvent->setDebugFlag(0);
  pEvent->setDebugFlag1(0);


  pRootFile=0;
  if(nEvt>0)
  {
      eventLoop(nEvt,0);
  }

  delete pEvent;
  //delete pRootFile;
}
//______________________________________________________________________________
void Unpacker::unpackerFast(const char* dir,const char* name, const char* odir, Int_t nEvt, TString luptab, TString calpar)
//const char* subEvtId,
//		   Int_t refChannel, const char* fpga_code, Int_t min, Int_t max, Int_t quietMode,
//		   Int_t fullSetup, Int_t VHR)
{
  refCh = -1;
  refCh = 31;
  this->subEvtId= HexStrToInt("999");
  outputFile = Form("%s%s.root",odir,name);
  cout<<"ofile "<<outputFile<<endl;
  setInputFile(dir,name);
  EventNr=0;
  pEvent= new HldEvent(inputFile.c_str(), HexStrToInt("999"), "999", -1000000, 1000000);

  setFileLookupPar(luptab);
  setFileHitFinderPar(calpar);

  pEvent->setQuietMode(true);
  //pEvent->setQuietMode(false);
  pEvent->setFullSetup(true);
  //pEvent->setFullSetup(false);
  //if(VHR == 1) {
  // pEvent->setVHR(true);
  //}
  //else if(VHR == 0) {
  pEvent->setVHR(false);
  //}

  pEvent->setDebugFlag(0);
  pEvent->setDebugFlag1(0);


  pRootFile=0;
  if(nEvt>0)
  {
      eventLoop(nEvt,0);
  }

  //delete pEvent;
  //delete pRootFile;
}
//______________________________________________________________________________
Unpacker::~Unpacker() 
{
  if(pEvent)	delete pEvent;
  if(pRootFile)	delete pRootFile;

  /*
  for(Int_t i=0;i<3;i++) {
	for(Int_t j=0;j<10;j++) {
	    for(Int_t k=0;k<12;k++) {
		delete hq[i*10*12+j*12+k];
		delete hdta[i*10*12+j*12+k];
	    }
	}
    }
    for(Int_t i=0;i<10;i++) {
	for(Int_t j=0;j<12;j++) {
            for(Int_t k=0;k<10;k++) {
		for(Int_t l=0;l<12;l++) {
		    delete hdt[i*12*10*12+j*10*12+k*12+l];

		}
	    }
	}
    }
    delete hdt[10*12*10*12];

    delete h1D[0];
    delete h1D[1];
    delete h1D[2];
    delete h1D[3];
    delete h1D[4];
    delete h1D[5];
    delete h1D[6];
    delete h1D[7];
    delete h1D[8];
    delete h1D[9];

    delete h2D[0];
    delete h2D[1];
    delete h2D[2];
    delete h2D[3];
    delete h2D[4];
    delete h2D[5];
    delete h2D[6];
    delete h2D[7];
    */























}

//______________________________________________________________________________
Bool_t Unpacker::setRootFile(const char* filename/*="test.root" */)
{
  if(pRootFile)
  {
	  delete pRootFile;
	  pRootFile=new TFile(filename,"RECREATE");	  
  }
  else	pRootFile=new TFile(filename,"RECREATE");	  

return kTRUE;
}
 
//______________________________________________________________________________
string Unpacker::setInputFile(const char* dir,const char* filename)
{
    stringstream Name(Form("%s%s",dir,filename));
   string Tmp,NoWhite;
   //it strips all whitespaces from name
   while(!Name.eof())
   {
 	Tmp.clear();
 	Name>>skipws >>Tmp;
 	NoWhite+=Tmp;
   }
  inputFile=NoWhite;  
   return inputFile;
}

string Unpacker::setInputFile(const char* filename)
{
   stringstream Name(filename);
   string Tmp,NoWhite;
   //it strips all whitespaces from name
   while(!Name.eof())
   {
 	Tmp.clear();
 	Name>>skipws >>Tmp;
 	NoWhite+=Tmp;
   }
  inputFile=NoWhite;  
   return inputFile;
}
	
//______________________________________________________________________________
Bool_t Unpacker::setpEvent(Int_t Id)
//Id is the subevent id 
{
   if(inputFile.empty())
   {
      return kFALSE;
   }
   subEvtId=Id;
   return kTRUE;
}

//______________________________________________________________________________
Bool_t Unpacker::eventLoop(Int_t nbEvt,Int_t startEv)
//Loop over all events, data written to the root tree
{
    if(pEvent==0)
    {
	cout<<"Error: no pEvent set"<<endl;
	return kFALSE;
    }
    else
    {
	//if(!pRootFile)
	if(1)
	{
	    char* t = new char[(outputFile + ".root").length() + 1];
	    strcpy(t, (outputFile + ".root").c_str());
	    setRootFile(t);
	}


	// TClonesArray* arr = new TClonesArray("Hit");
	// En el caso de monitoring NO NOS INTERESA CREAR EL TREE;
	// ES MAS SENCILLO HACER OTRO EVENT LOOP eventLoopMonitor(start, stop)....
        //RpcLookupTable* look= new RpcLookupTable("luptab.txt");


	Int_t trbnum, cell, col, raw;
	Float_t     x, y, z;

	//look->getParams(837, 0, 1, cell, col, raw, x, y, z);

	TTree *tree = new TTree("T","Tree");
	tree->SetAutoSave(1000000000); //autosave when 1 GB written
	Event* event = new Event();

	Int_t split = 2;
	Int_t bsize = 64000000;
	cout<< " ********************************************************* "<<endl;
	cout<< " Tree T is being created  "<<endl;
	cout<< " All the available branches are defined in the event class "<<endl;
	cout<< " ********************************************************* "<<endl;


	TRpcRawF* rawFinder = new TRpcRawF();
	rawFinder -> init(fileLookupPar);

	TRpcHitF* calFinder = new TRpcHitF();
	calFinder -> init(fileHitFinderPar,"/media/Datos2TB/damian/tragaldabas/soft_TT/GoodActiveCells.root");

	TRpcSaetaF* trackFinder = new TRpcSaetaF();
        trackFinder->init();


	TClonesArray*  rpchits          = rawFinder  ->getRpcRawHits();
	TClonesArray*  rpccalhits       = calFinder  ->getRpcHits();
	TClonesArray*  rpccalhitscorr   = trackFinder->getRpcHitCorr();
	TClonesArray*  RpcSaeta2Planes  = trackFinder->getRpcSaeta2Planes();
	TClonesArray*  RpcSaeta3Planes  = trackFinder->getRpcSaeta3Planes();

	//cout<<"rpchits1 "<<rpchits<<endl;
        
        TClonesArray* hits = event->getHits();
	//event = new Event();

        // Full output!
	tree->Branch("event","Event", &event, bsize,split);
	tree->Branch("hits","TClonesArray",&hits,bsize,split);
	tree->Branch("rpcraw","TClonesArray", &rpchits,bsize,split);
	tree->Branch("rpchit","TClonesArray", &rpccalhits,bsize,split);
	tree->Branch("RpcSaeta2Planes","TClonesArray", &RpcSaeta2Planes,bsize,split);
	tree->Branch("RpcSaeta3Planes","TClonesArray", &RpcSaeta3Planes,bsize,split);
	tree->Branch("rpccalhitscorr","TClonesArray", &rpccalhitscorr,bsize,split);
	//tree->DropBranchFromCache("Hits");
	//rpchits->BypassStreamer();
	Bool_t sync = kTRUE;

	for(Int_t i=0; i< nbEvt; i++)
	{
	    if(!(pEvent->execute())) // The important business happens here!
	    {
		cout<<"END OF FILE"<<endl;
		cout<<"Number of Events: "<<EventNr<<endl;
		break;
	    }
            //ANADIR AQUI EL CONTINUE PARA SALTARSE EVENTOS, especialmente para monitoring!
	    if(!event)
		event = new Event(*pEvent, refCh);
            event->clearAll();
	    event->fill(*pEvent);
            // This is for streaming it to the TREE
            event->setSync(sync);
	    rawFinder -> execute();
            sync = event->getSync();
            //rpchits = rawFinder->getRpcRawHits();
	    //cout<<" rpcraw "<<rpchits->GetEntriesFast()<<endl;
	    calFinder -> execute();
            trackFinder -> execute();

            //if(rh)
	    tree->Fill();
	    //delete event;
	    EventNr++;
	}

	pRootFile->Write();
        pRootFile->Clone();
	delete tree;
        delete event;
	delete rawFinder;
        delete trackFinder;
	delete calFinder;
        delete pRootFile;
	return kTRUE;
    }
}
Int_t Unpacker::eventLoopSyncCheck(Int_t nbEvt,Int_t startEv)
//Loop over all events, data written to the root tree
{
    if(pEvent==0)
    {
	cout<<"Error: no pEvent set"<<endl;
	return 1;
    }
    else
    {
	Int_t trbnum, cell, col, raw;
	Float_t     x, y, z;

	Event* event = new Event();

	TRpcRawF* rawFinder = new TRpcRawF();
	rawFinder -> init(fileLookupPar);

        TClonesArray* hits = event->getHits();
	Bool_t sync = kTRUE;
	for(Int_t i=0; i< nbEvt; i++)
	{
	    if(!(pEvent->execute())) // The important business happens here!
	    {
		break;
	    }
            //ANADIR AQUI EL CONTINUE PARA SALTARSE EVENTOS, especialmente para monitoring!
	    if(!event)
		event = new Event(*pEvent, refCh);
            event->clearAll();
	    event->fill(*pEvent);
            event->setSync(sync);
	    rawFinder -> execute();
            sync = event->getSync();
            if(!sync) return 1;
	    EventNr++;
	}

        delete event;
	delete rawFinder;
	if(sync) return 0;
	else return 1;
    }
    return 1;
}



Bool_t Unpacker::eventLoopFillCal(Int_t nbEvt, Int_t startEv, TH1D** hq, TH1D** hdt,TH1D** hdt2)
//Loop over all events, data written to the root tree
{
    if(pEvent==0)
    {
	cout<<"Error: no pEvent set"<<endl;
	return kFALSE;
    }
    else
    {

	//if(!pRootFile)
	//{
	//    char* t = new char[(outputFile + ".root").length() + 1];
	//    strcpy(t, (outputFile + ".root").c_str());
	//    setRootFile(t);
	//}

	// TClonesArray* arr = new TClonesArray("Hit");
	// En el caso de monitoring NO NOS INTERESA CREAR EL TREE;
	// ES MAS SENCILLO HACER OTRO EVENT LOOP eventLoopMonitor(start, stop)....
        //RpcLookupTable* look= new RpcLookupTable(fileLookupPar);


	Int_t trbnum, cell, col, row;
	Float_t     x, y, z;
        Float_t time, charge;

	Int_t trbnum2, cell2, col2, row2;
	Float_t     x2, y2, z2;
        Float_t time2, charge2;

	Int_t trbnum3, cell3, col3, row3;
	Float_t     x3, y3, z3;
        Float_t time3, charge3;

	//look->getParams(837, 0, 1, cell, col, raw, x, y, z);

	//TTree *tree = new TTree("T","Tree");
	//tree->SetAutoSave(1000000000); //autosave when 1 GB written
//	Event* event = new Event(*pEvent, refCh);
	Event* event = new Event();

	//Int_t split = 2;
	//Int_t bsize = 64000;
	//cout<< " ********************************************************* "<<endl;
	//cout<< " Tree T is being created  "<<endl;
	//cout<< " All the available branches are defined in the event class "<<endl;
	//cout<< " ********************************************************* "<<endl;

	cout << "fileLookupPar " <<fileLookupPar<<endl;
	TRpcRawF* rawFinder = new TRpcRawF();
	rawFinder -> init(fileLookupPar);

        /*
	cout << "fileHitPar "<<fileHitFinderPar <<endl;
	TRpcHitF* calFinder = new TRpcHitF();
	calFinder -> init(fileHitFinderPar,"/media/Datos2TB/korna/tragaldabas/pars/active_cells_2015_2016.root");
	//calFinder -> init("parhitfNew4.txt");
        */


        cout<<"params loaded "<<endl;
	TClonesArray*  rpchits    = rawFinder->getRpcRawHits();

        //TClonesArray*  rpccalhits = calFinder->getRpcHits();

	//cout<<"rpchits1 "<<rpchits<<endl;
	TClonesArray* hits = event->getHits();
	//event = new Event();

	//tree->Branch("hits","TClonesArray",&hits,bsize,split);
        // Full output!
	//tree->Branch("event","Event", &event, bsize,split);
	//tree->Branch("rpcraw","TClonesArray", &rpchits);
	//tree->DropBranchFromCache("Hits");
	//rpchits->BypassStreamer();


        TRandom3 rand;
        rand.SetSeed(0);


	for(Int_t i=0; i< nbEvt; i++)
	{
	    if(!(pEvent->execute())) // The important business happens here!
	    {
		cout<<"END OF FILE"<<endl;
		cout<<"Number of Events: "<<EventNr<<endl;
		break;
	    }
            //ANADIR AQUI EL CONTINUE PARA SALTARSE EVENTOS, especialmente para monitoring!
	    if(!event)
		event = new Event(*pEvent, refCh);

            
	    event->clearAll();
            
	    event->fill(*pEvent);
            // This is for streaming it to the TREE

	    rawFinder -> execute();
            //calFinder -> execute();
	    rpchits = rawFinder->getRpcRawHits();
	    Int_t nrpchits = rpchits->GetEntriesFast();
	    //cout<<"rpchits2 "<<rpchits<<endl;
            Int_t dnum0=0;
            Int_t dnum1=0;
            Int_t dnum2=0;

	    TRpcRaw* rh = 0;
	    TRpcRaw* rh2 = 0;
	    for(Int_t n=0;n<nrpchits;n++) {
		rh = (TRpcRaw*)rpchits->At(n);
		rh->getHit(trbnum, cell,col,row,x,y,z,time,charge);
		row-=1;
		col-=1;
		//if(nrpchits>3)
		for(Int_t n2=0;n2<nrpchits;n2++) {
                    if(n==n2)continue;
		    rh2 = (TRpcRaw*)rpchits->At(n2);
		    rh2->getHit(trbnum2, cell2,col2,row2,x2,y2,z2,time2,charge2);
		    row2-=1;
		    col2-=1;
		    if(charge2>1000&& TMath::Sqrt((row2-row)*(row2-row) + (col2-col)*(col2-col))<1.1)
			hq[trbnum*10*12+row*12+col] -> Fill(charge*0.098);
		}
		if(charge>600) {
		    if(trbnum==0) dnum0++;
		    if(trbnum==1) dnum1++;
		    if(trbnum==2) dnum2++;
		}
	    }

            if(dnum0==1 && dnum1==1) {
                Float_t x0,y0,z0,time0,q0,x1,y1,z1,time1,q1,q2;
                Int_t col0,row0,col1,row1;
                for(Int_t n=0;n<nrpchits;n++) {
                    TRpcRaw* rh = 0;
                    //TRpcRaw* rh2 = 0;
                    //for(Int_t n=0;n<nrpchits;n++) {
                    rh = (TRpcRaw*)rpchits->At(n);
                    rh->getHit(trbnum, cell,col,row,x,y,z,time,charge);
                    if(charge>600 && trbnum==0) {
                        x0 = x;
                        y0 = y;
                        z0 = z;
                        time0 = time*0.098 + 0.098*(rand.Rndm()-0.5);
                        col0=col-1;
                        row0=row-1;
                        q0=charge;
                    }
                    if(charge>600 && trbnum==1) {
                        x1 = x;
                        y1 = y;
                        z1 = z;
                        time1 = time*0.098 + 0.098*(rand.Rndm()-0.5);
                        col1=col-1;
                        row1=row-1;
                        q1=charge;
                    }
                }
                hdt[row0*12*10*12+col0*10*12+row1*12+col1]->Fill(time1 - time0 - TMath::Sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))/299.792458);
                if(q0>100 &&q1>100)
                    hdt[10*12*10*12]->Fill(time1 - time0 - TMath::Sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))/299.792458);
            }


            if(dnum2==1 && dnum1==1) {
                Float_t x0,y0,z0,time0,q0,x1,y1,z1,time1,q1,q2;
                Int_t col0,row0,col1,row1;
                for(Int_t n=0;n<nrpchits;n++) {
                    TRpcRaw* rh = 0;
                    rh = (TRpcRaw*)rpchits->At(n);
                    rh->getHit(trbnum, cell,col,row,x,y,z,time,charge);
                    if(charge>600 && trbnum==2) {
                        x0 = x;
                        y0 = y;
                        z0 = z;
                        time0 = time*0.098 + 0.098*(rand.Rndm()-0.5);
                        col0=col-1;
                        row0=row-1;
                        q0=charge;
                    }
                    if(charge>600 && trbnum==1) {
                        x1 = x;
                        y1 = y;
                        z1 = z;
                        time1 = time*0.098 + 0.098*(rand.Rndm()-0.5);
                        col1=col-1;
                        row1=row-1;
                        q1=charge;
                    }
                }
                hdt2[row0*12*10*12+col0*10*12+row1*12+col1]->Fill(time1 - time0 - TMath::Sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))/299.792458);
                if(q0>100 &&q1>100)
                    hdt2[10*12*10*12]->Fill(time1 - time0 - TMath::Sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))/299.792458);
            }







            /*



            TRpcHit* rhit;
	    rpchits = calFinder->getRpcHits();
            nrpchits = rpchits->GetEntriesFast();
	    
	    if(event->getMultH1()==1 && event->getMultH2() == 1) {
		for(Int_t n=0;n<nrpchits;n++) {
		    rhit = (TRpcHit*)rpchits->At(n);
		    rhit->getHit(trbnum, cell,col,row,x,y,z,time,charge);
		    if(trbnum==0) {
			//rh0 = rh;
			x0 = x;
			y0 = y;
			z0 = z;
			time0 = time;
			col0=col-1;
			row0=row-1;
                        q0=charge;
		    }
		    if(trbnum==1) {
			//rh1 = rh;
			x1 = x;
			y1 = y;
			z1 = z;
			time1 = time;
			col1=col-1;
			row1=row-1;
                        q1= charge;
		    }
		}
		//if(q1>620&&q0>620)
		hdt[row0*12*10*12+col0*10*12+row1*12+col1]->Fill(time1 - time0 - TMath::Sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))/299.792458);
                if(q0>100 &&q1>100)
		    hdt[10*12*10*12]->Fill(time1 - time0 - TMath::Sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))/299.792458);

	    }

	    if(event->getMultH3()==1 && event->getMultH2() == 1) {
		for(Int_t n=0;n<nrpchits;n++) {
		    rhit = (TRpcHit*)rpchits->At(n);
		    rhit->getHit(trbnum, cell,col,row,x,y,z,time,charge);
		    if(trbnum==2) {
			//rh0 = rh;
			x2 = x;
			y2 = y;
			z2 = z;
			time2 = time;
			col2=col-1;
			row2=row-1;
                        q2=charge;
		    }
		    if(trbnum==1) {
			//rh1 = rh;
			x1 = x;
			y1 = y;
			z1 = z;
			time1 = time;
			col1=col-1;
			row1=row-1;
                        q1= charge;
		    }
		}
		//if(q1>620&&q0>620)
		hdt2[row2*12*10*12+col2*10*12+row1*12+col1]->Fill(time1 - time2 - TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))/299.792458);
                if(q2>100 &&q1>100)
		    hdt2[10*12*10*12]->Fill(time1 - time2 - TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))/299.792458);

	    }
            */
	    //cout<<"pointer rh "<<rh<<endl;
	    //Float_t time,charge;
	    //if(rh) {
	    //rh->getHit(trbnum, cell,col,raw,x,y,z,time,charge);
	    //cout<<"t q "<<time<<" "<<charge<<endl;
	    //}
	    //TClonesArray* arr = event->getHits();
            //cout<<"POINTER UNPACKER "<<arr<<endl;
	    //Int_t arrentr = arr->GetEntriesFast();
            //cout<<" ENTRIES FAST "<<arrentr<<endl;
	    //tree->Fill();
	    //delete event;
	    EventNr++;
	}

	//pRootFile->Write();

	//delete tree;
	//delete tree;
        delete event;
        delete rawFinder;
	//delete calFinder;
        //delete pRootFile;





	return kTRUE;
    }
}
Bool_t Unpacker::eventLoopFillCal(Int_t nbEvt, Int_t startEv, TH1D** h1D, TH2D** h2D, TH3D** h3D)
//Loop over all events, data written to the root tree
{
    if(pEvent==0)
    {
	cout<<"Error: no pEvent set"<<endl;
	return kFALSE;
    }
    else
    {

	//if(!pRootFile)
	//{
	//    char* t = new char[(outputFile + ".root").length() + 1];
	//    strcpy(t, (outputFile + ".root").c_str());
	//    setRootFile(t);
	//}

	// TClonesArray* arr = new TClonesArray("Hit");
	// En el caso de monitoring NO NOS INTERESA CREAR EL TREE;
	// ES MAS SENCILLO HACER OTRO EVENT LOOP eventLoopMonitor(start, stop)....

	cout<<"loading lookup "<<fileLookupPar<<endl;
	//RpcLookupTable* look= new RpcLookupTable(fileLookupPar);


	Int_t trbnum, cell, col, row;
	Float_t     x, y, z;
        Float_t time, charge;

	Int_t trbnum2, cell2, col2, row2;
	Float_t     x2, y2, z2;
        Float_t time2, charge2;

	//look->getParams(837, 0, 1, cell, col, raw, x, y, z);

	//TTree *tree = new TTree("T","Tree");
	//tree->SetAutoSave(1000000000); //autosave when 1 GB written
//	Event* event = new Event(*pEvent, refCh);
	Event* event = new Event();

	//Int_t split = 2;
	//Int_t bsize = 64000;
	//cout<< " ********************************************************* "<<endl;
	//cout<< " Tree T is being created  "<<endl;
	//cout<< " All the available branches are defined in the event class "<<endl;
	//cout<< " ********************************************************* "<<endl;


	TRpcRawF* rawFinder = new TRpcRawF();
	rawFinder -> init(fileLookupPar);

	TRpcHitF* calFinder = new TRpcHitF();
	calFinder -> init(fileHitFinderPar,"/media/Datos2TB/damian/tragaldabas/soft_TT/active_cells_2015_2016.root");

	TRpcSaetaF* trackFinder = new TRpcSaetaF();
	trackFinder->init();

	TClonesArray*  rpchits          = rawFinder->getRpcRawHits();
        TClonesArray*  rpccalhits       = calFinder->getRpcHits();
	TClonesArray*  RpcSaeta2Planes  = trackFinder->getRpcSaeta2Planes();
	TClonesArray*  RpcSaeta3Planes  = trackFinder->getRpcSaeta3Planes();
	TClonesArray*  rpccalhitscorr   = trackFinder->getRpcHitCorr();

	//cout<<"rpchits1 "<<rpchits<<endl;
        
        TClonesArray* hits = event->getHits();
	//event = new Event();

	//tree->Branch("hits","TClonesArray",&hits,bsize,split);
        // Full output!
	//tree->Branch("event","Event", &event, bsize,split);
	//tree->Branch("rpcraw","TClonesArray", &rpchits);
	//tree->DropBranchFromCache("Hits");
	//rpchits->BypassStreamer();
	Bool_t sync = kTRUE;
	Int_t prevsec = -100;
	for(Int_t i=0; i< nbEvt; i++)
	{
	    if(!(pEvent->execute())) // The important business happens here!
	    {
		cout<<"END OF FILE"<<endl;
		cout<<"Number of Events: "<<EventNr<<endl;
		break;
	    }
            //ANADIR AQUI EL CONTINUE PARA SALTARSE EVENTOS, especialmente para monitoring!
	    if(!event)
		event = new Event(*pEvent, refCh);

            event->clearAll();
	    event->fill(*pEvent);
            // This is for streaming it to the TREE



	    event->setSync(sync);
	    rawFinder -> execute();
            sync = event->getSync();
	    calFinder -> execute();

            trackFinder -> execute();





            h2D[0]->Fill(event->getMultH1(),event->getMultH2());

	    Float_t timeEvt = TTimeStamp::GetDayOfYear(event->getEvtDay(),event->getEvtMonth(),event->getEvtYear());
            timeEvt+= event->getEvtHour()/24.;
            timeEvt+= event->getEvtMinute()/24./60.;

	    if(event->getMultH1()==1 && event->getMultH2()==1) {
		
	    }
	    if(prevsec != event->getEvtSecond() ) {
		prevsec = event->getEvtSecond();
		h1D[9]->Fill(timeEvt);
	    }

	    rpccalhits = calFinder->getRpcHits();
	    Int_t nrpchits = rpccalhits->GetEntriesFast();
	    TRpcHit* rh = 0;
	    TRpcHit* rh2 = 0;


	    Float_t x0,y0,z0,time0,q0,x1,y1,z1,time1,q1;
	    Int_t col0,row0,col1,row1;


            //h1D[10] = new TH1D("h_Q_M1","h
            //h1D[11] = new TH1D("h_Q_N","h_
            //h1D[12] = new TH1D("h_Flux_M1"
            //h1D[13] = new TH1D("h_Flux_M2"
            //h1D[14] = new TH1D("h_Flux_M3"

	    Int_t ntracks = RpcSaeta2Planes->GetEntriesFast();
	    //Int_t ntracks = RpcSaeta3Planes->GetEntriesFast();

            Int_t theta = event->getThetaDeg();
	    Int_t phi   = event->getPhiDeg();


	    Int_t thetaBin = floor(theta/10.);
            if(thetaBin>4) thetaBin = 4;
            if(thetaBin<0) continue;


	    if(ntracks==1 && event->getMultH1()==1 &&event->getMultH2()==1) {

		h2D[thetaBin+8]->Fill(timeEvt,phi);

		h1D[12]->Fill(timeEvt);
		if(theta>0.) {
		    if(theta<30.) {
			h1D[15]->Fill(timeEvt);
		    }
		    else {
			h1D[16]->Fill(timeEvt);
		    }
		}
	    }
	    if(ntracks>1 && event->getMultH1()==2 && event->getMultH2()==2) {
		h1D[13]->Fill(timeEvt);

		h2D[thetaBin+8+5]->Fill(timeEvt,phi);

		if(theta>0.) {
		    if(theta<30.) {
			h1D[17]->Fill(timeEvt);
		    }
		    else {
			h1D[18]->Fill(timeEvt);
		    }
		}
		//if(event->getPhiDeg()>0. && event->getPhiDeg()<60.)

	    }
	    if(ntracks>2 && event->getMultH1()>2 && event->getMultH2()>2) {
		h1D[14]->Fill(timeEvt);

		h2D[thetaBin+8+5+5]->Fill(timeEvt,phi);

		if(theta>0.) {
		    if(theta<30.) {
			h1D[19]->Fill(timeEvt);
		    }
		    else {
			h1D[20]->Fill(timeEvt);
		    }
		}
	    }

	    if(event->getMultH1()==1 &&event->getMultH2()==1) {
		for(Int_t n=0;n<nrpchits;n++) {
		    rh = (TRpcHit*)rpccalhits->At(n);
		    rh->getHit(trbnum, cell,col,row,x,y,z,time,charge);

		    if(trbnum==0) {
			if(charge>0 && charge<1000) {
			    h1D[10] -> Fill(timeEvt,charge);
			    h1D[11] -> Fill(timeEvt);
			}
			//rh0 = rh;
			x0 = x;
			y0 = y;
			z0 = z;
			time0 = time;
			col0=col-1;
			row0=row-1;
                        q0=charge;
		    }
		    if(trbnum==1) {
			//rh1 = rh;
			x1 = x;
			y1 = y;
			z1 = z;
			time1 = time;
			col1=col-1;
			row1=row-1;
                        q1= charge;
		    }
		}
		//if(col1==col0)
		Float_t tdif = time1 - time0 - TMath::Sqrt( (x1-x0)*(x1-x0)
						    +(y1-y0)*(y1-y0)
						    +(z1-z0)*(z1-z0))/299.792458;

		h2D[5]->Fill(time1 - time0 - TMath::Sqrt( (x1-x0)*(x1-x0)
					      +(y1-y0)*(y1-y0)
						  +(z1-z0)*(z1-z0))/299.792458,
			     TMath::Sqrt( (x1-x0)*(x1-x0)
				  +(y1-y0)*(y1-y0)
				  +(z1-z0)*(z1-z0)));





		if(TMath::Abs(tdif)<2.) {
		    h1D[8]->Fill(timeEvt);
		}
		/*
		h2D[5]->Fill(time0-time1,
			     TMath::Sqrt( (x1-x0)*(x1-x0)
				  +(y1-y0)*(y1-y0)
				  +(z1-z0)*(z1-z0)));
                */
		if(q1>0 && q0>0){
		    h1D[5]->Fill(time1 - time0-TMath::Sqrt( (x1-x0)*(x1-x0)
						  +(y1-y0)*(y1-y0)
						    +(z1-z0)*(z1-z0))/299.792458 );



                    h2D[6]  -> Fill(time1 - time0-TMath::Sqrt( (x1-x0)*(x1-x0)
						  +(y1-y0)*(y1-y0)
						    +(z1-z0)*(z1-z0))/299.792458,q0);
                    h2D[7]  -> Fill(time1 - time0-TMath::Sqrt( (x1-x0)*(x1-x0)
						  +(y1-y0)*(y1-y0)
						    +(z1-z0)*(z1-z0))/299.792458,log10(q0));


		    //if(col1==col0&&row1==row0){
		}

                Float_t wcorr = 0.;

		if(q1>400&&q1<900) {
		    //if(q0<150)  wcorr = exp(q0*-0.0111727+0.935106)+0.628899-0.0478589;
		    //if(q0>=150) wcorr = exp(q0*-0.00373208+0.63251);

		    h2D[3]->Fill(time1 - (time0) -TMath::Sqrt( (x1-x0)*(x1-x0)
						    +(y1-y0)*(y1-y0)
						    +(z1-z0)*(z1-z0))/299.792458,q0);
		}
		if(q0>400&&q0<900) {

		    //if(q1<150)  wcorr = exp(q1*-0.0111727+0.935106)+0.628899-0.0478589;
		    //if(q1>=150) wcorr = exp(q1*-0.00373208+0.63251);


		    h2D[4]->Fill(time1 - time0-TMath::Sqrt( (x1-x0)*(x1-x0)
						    +(y1-y0)*(y1-y0)
						    +(z1-z0)*(z1-z0))/299.792458,q1);
		}
		//}

	    }



            /*
	    Int_t dnum0=0;
            Int_t dnum1=0;
            Int_t dnum2=0;

            
	    TRpcHit* rh = 0;
	    TRpcHit* rh2 = 0;
	    for(Int_t n=0;n<nrpchits;n++) {
		rh = (TRpcHit*)rpccalhits->At(n);
		rh->getHit(trbnum, cell,col,row,x,y,z,time,charge);
		row-=1;
		col-=1;
		//if(nrpchits>3)
		for(Int_t n2=0;n2<nrpchits;n2++) {
                    if(n==n2)continue;
		    rh2 = (TRpcRaw*)rpccalhits->At(n2);
		    rh2->getHit(trbnum2, cell2,col2,row2,x2,y2,z2,time2,charge2);
		    row2-=1;
		    col2-=1;
		    if(charge2>1000&& TMath::Sqrt((row2-row)*(row2-row) + (col2-col)*(col2-col))<2)
			hq[trbnum*10*12+row*12+col] -> Fill(charge*0.098);
		}
		if(charge>620) {
		    if(trbnum==0) dnum0++;
		    if(trbnum==1) dnum1++;
		    if(trbnum==2) dnum2++;
		}
	    }

	    Float_t x0,y0,z0,time0,q0,x1,y1,z1,time1,q1;
            Int_t col0,row0,col1,row1;
	    if(dnum0==1 && dnum1 == 1) {
		for(Int_t n=0;n<nrpchits;n++) {
		    rh = (TRpcRaw*)rpchits->At(n);
		    rh->getHit(trbnum, cell,col,row,x,y,z,time,charge);
		    if(trbnum==0) {
			//rh0 = rh;
			x0 = x;
			y0 = y;
			z0 = z;
			time0 = time;
			col0=col-1;
			row0=row-1;
                        q0=charge;
		    }
		    if(trbnum==1) {
			//rh1 = rh;
			x1 = x;
			y1 = y;
			z1 = z;
			time1 = time;
			col1=col-1;
			row1=row-1;
                        q1= charge;
		    }
		}
		if(q1>620&&q0>620)
		hdt[row0*12*10*12+col0*10*12+row1*12+col1]->Fill(time0*0.098
						  - time1*0.098
						  -TMath::Sqrt( (x1-x0)*(x1-x0)
							+(y1-y0)*(y1-y0)
							+(z1-z0)*(z1-z0))/299.792458);

	    }
            */
	    //cout<<"pointer rh "<<rh<<endl;
	    //Float_t time,charge;
	    //if(rh) {
	    //rh->getHit(trbnum, cell,col,raw,x,y,z,time,charge);
	    //cout<<"t q "<<time<<" "<<charge<<endl;
	    //}
	    //TClonesArray* arr = event->getHits();
            //cout<<"POINTER UNPACKER "<<arr<<endl;
	    //Int_t arrentr = arr->GetEntriesFast();
            //cout<<" ENTRIES FAST "<<arrentr<<endl;
	    //tree->Fill();
	    //delete event;
	    EventNr++;
	}

	//pRootFile->Write();
        delete event;
        delete rawFinder;
	delete calFinder;
        //delete pRootFile;

	//delete tree;
	return kTRUE;
    }
}
