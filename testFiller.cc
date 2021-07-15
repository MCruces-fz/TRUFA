#include "tevent.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "trpchit.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include "TF1.h"

void Analysis(Char_t* files){
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Wed Jun 15 18:13:29 2016 by ROOT version5.18/00b)
//   from TTree T/Tree
//   found on file: tr16160092743.hld.root.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tr16160092743.hld.root.root");
   //if (!f) {
   //   f = new TFile("tr16160092743.hld.root.root");
   //}
   TChain * T= new TChain("T","T");
   T->Add(files);
   //T->Print();
   //TTree *T = (TTree*)gDirectory->Get("T");

   //Declaration of leaves types
   Event* event = new Event();
   TClonesArray* hits = new TClonesArray("TRpcHit");
   /*
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UInt_t          EvtHdr_fUniqueID;
   UInt_t          EvtHdr_fBits;
   Int_t           EvtHdr_size;
   Int_t           EvtHdr_decoding;
   Int_t           EvtHdr_id;
   Int_t           EvtHdr_seqNr;
   Int_t           EvtHdr_date;
   Int_t           EvtHdr_time;
   Int_t           EvtHdr_year;
   Int_t           EvtHdr_month;
   Int_t           EvtHdr_day;
   Int_t           EvtHdr_hour;
   Int_t           EvtHdr_minute;
   Int_t           EvtHdr_second;
   Int_t           EvtHdr_pad;
   Int_t           EvtHdr_dataSize;
   Int_t           EvtHdr_paddedSize;
   Int_t           kMaxMult;
   Int_t           kMaxChannelNr;
   UInt_t          subEvtId;
   Int_t           errors_per_event;
   Int_t           referenceChannel;
   Int_t           referenceTime;
   Int_t           totalNHits;
   */
   Int_t           multH1;
   Int_t           multH2;
   Int_t           multH3;
   /*
   Bool_t          sync;
   Float_t         fEvAl;
   Float_t         fEvBe;
   Float_t         fEvGa;
   UInt_t          hits_fUniqueID[365];
   UInt_t          hits_fBits[365];
   Int_t           hits_channel[365];
   Int_t           hits_TDC[365];
   Int_t           hits_nHits[365];
   Int_t           hits_trbNum[365];
   Int_t           hits_width[365];
   Int_t           hits_leadTime1[365];
   Int_t           hits_trailTime1[365];
   Int_t           hits_leadTime2[365];
   Int_t           hits_trailTime2[365];
   Int_t           hits_leadTime3[365];
   Int_t           hits_trailTime3[365];
   Int_t           hits_leadTime4[365];
   Int_t           hits_trailTime4[365];
   Int_t           hits_spikesCtr[365];
   UInt_t          rpcraw_fUniqueID[352];
   UInt_t          rpcraw_fBits[352];
   Int_t           rpcraw_fTrbnum[352];
   Int_t           rpcraw_fCell[352];
   Int_t           rpcraw_fCol[352];
   Int_t           rpcraw_fRow[352];
   Float_t         rpcraw_fX[352];
   Float_t         rpcraw_fY[352];
   Float_t         rpcraw_fZ[352];
   Float_t         rpcraw_fTime[352];
   Float_t         rpcraw_fCharge[352];
   UInt_t          rpchit_fUniqueID[226];
   UInt_t          rpchit_fBits[226];
   Int_t           rpchit_fTrbnum[226];
   Int_t           rpchit_fCell[226];
   Int_t           rpchit_fCol[226];
   Int_t           rpchit_fRow[226];
   Float_t         rpchit_fX[226];
   Float_t         rpchit_fY[226];
   Float_t         rpchit_fZ[226];
   Float_t         rpchit_fTime[226];
   Float_t         rpchit_fCharge[226];
   UInt_t          tracks_fUniqueID[1294];
   UInt_t          tracks_fBits[1294];
   Int_t           tracks_find0[1294];
   Int_t           tracks_find1[1294];
   Float_t         tracks_fAl[1294];
   Float_t         tracks_fBe[1294];
   Float_t         tracks_fGa[1294];
   Float_t         tracks_fX[1294];
   Float_t         tracks_fY[1294];
   Float_t         tracks_fZ[1294];
   Float_t         tracks_fTime[1294];
   UInt_t          rpccalhitscorr_fUniqueID[226];
   UInt_t          rpccalhitscorr_fBits[226];
   Int_t           rpccalhitscorr_fTrbnum[226];
   Int_t           rpccalhitscorr_fCell[226];
   Int_t           rpccalhitscorr_fCol[226];
   Int_t           rpccalhitscorr_fRow[226];
   Float_t         rpccalhitscorr_fX[226];
   Float_t         rpccalhitscorr_fY[226];
   Float_t         rpccalhitscorr_fZ[226];
   Float_t         rpccalhitscorr_fTime[226];
   Float_t         rpccalhitscorr_fCharge[226];
   */
   // Set branch addresses.
   //T->SetBranchStatus("*",0);
   T->SetBranchAddress("event",&event);
   T->SetBranchAddress("rpchit",&hits);
   /*
   T->SetBranchAddress("fUniqueID",&fUniqueID);
   T->SetBranchAddress("fBits",&fBits);
   T->SetBranchAddress("EvtHdr.fUniqueID",&EvtHdr_fUniqueID);
   T->SetBranchAddress("EvtHdr.fBits",&EvtHdr_fBits);
   T->SetBranchAddress("EvtHdr.size",&EvtHdr_size);
   T->SetBranchAddress("EvtHdr.decoding",&EvtHdr_decoding);
   T->SetBranchAddress("EvtHdr.id",&EvtHdr_id);
   T->SetBranchAddress("EvtHdr.seqNr",&EvtHdr_seqNr);
   T->SetBranchAddress("EvtHdr.date",&EvtHdr_date);
   T->SetBranchAddress("EvtHdr.time",&EvtHdr_time);
   T->SetBranchAddress("EvtHdr.year",&EvtHdr_year);
   T->SetBranchAddress("EvtHdr.month",&EvtHdr_month);
   T->SetBranchAddress("EvtHdr.day",&EvtHdr_day);
   T->SetBranchAddress("EvtHdr.hour",&EvtHdr_hour);
   T->SetBranchAddress("EvtHdr.minute",&EvtHdr_minute);
   T->SetBranchAddress("EvtHdr.second",&EvtHdr_second);
   T->SetBranchAddress("EvtHdr.pad",&EvtHdr_pad);
   T->SetBranchAddress("EvtHdr.dataSize",&EvtHdr_dataSize);
   T->SetBranchAddress("EvtHdr.paddedSize",&EvtHdr_paddedSize);
   T->SetBranchAddress("kMaxMult",&kMaxMult);
   T->SetBranchAddress("kMaxChannelNr",&kMaxChannelNr);
   T->SetBranchAddress("subEvtId",&subEvtId);
   T->SetBranchAddress("errors_per_event",&errors_per_event);
   T->SetBranchAddress("referenceChannel",&referenceChannel);
   T->SetBranchAddress("referenceTime",&referenceTime);
   T->SetBranchAddress("totalNHits",&totalNHits);
   T->SetBranchAddress("multH1",&multH1);
   T->SetBranchAddress("multH2",&multH2);
   T->SetBranchAddress("multH3",&multH3);
   T->SetBranchAddress("sync",&sync);
   T->SetBranchAddress("fEvAl",&fEvAl);
   T->SetBranchAddress("fEvBe",&fEvBe);
   T->SetBranchAddress("fEvGa",&fEvGa);
   T->SetBranchAddress("hits",&hits);
   T->SetBranchAddress("hits.fUniqueID",hits_fUniqueID);
   T->SetBranchAddress("hits.fBits",hits_fBits);
   T->SetBranchAddress("hits.channel",hits_channel);
   T->SetBranchAddress("hits.TDC",hits_TDC);
   T->SetBranchAddress("hits.nHits",hits_nHits);
   T->SetBranchAddress("hits.trbNum",hits_trbNum);
   T->SetBranchAddress("hits.width",hits_width);
   T->SetBranchAddress("hits.leadTime1",hits_leadTime1);
   T->SetBranchAddress("hits.trailTime1",hits_trailTime1);
   T->SetBranchAddress("hits.leadTime2",hits_leadTime2);
   T->SetBranchAddress("hits.trailTime2",hits_trailTime2);
   T->SetBranchAddress("hits.leadTime3",hits_leadTime3);
   T->SetBranchAddress("hits.trailTime3",hits_trailTime3);
   T->SetBranchAddress("hits.leadTime4",hits_leadTime4);
   T->SetBranchAddress("hits.trailTime4",hits_trailTime4);
   T->SetBranchAddress("hits.spikesCtr",hits_spikesCtr);
   //T->SetBranchAddress("rpcraw",&rpcraw_);
   T->SetBranchAddress("rpcraw.fUniqueID",rpcraw_fUniqueID);
   T->SetBranchAddress("rpcraw.fBits",rpcraw_fBits);
   T->SetBranchAddress("rpcraw.fTrbnum",rpcraw_fTrbnum);
   T->SetBranchAddress("rpcraw.fCell",rpcraw_fCell);
   T->SetBranchAddress("rpcraw.fCol",rpcraw_fCol);
   T->SetBranchAddress("rpcraw.fRow",rpcraw_fRow);
   T->SetBranchAddress("rpcraw.fX",rpcraw_fX);
   T->SetBranchAddress("rpcraw.fY",rpcraw_fY);
   T->SetBranchAddress("rpcraw.fZ",rpcraw_fZ);
   T->SetBranchAddress("rpcraw.fTime",rpcraw_fTime);
   T->SetBranchAddress("rpcraw.fCharge",rpcraw_fCharge);
   
   T->SetBranchAddress("rpchit.fUniqueID",rpchit_fUniqueID);
   T->SetBranchAddress("rpchit.fBits",rpchit_fBits);
   T->SetBranchAddress("rpchit.fTrbnum",rpchit_fTrbnum);
   T->SetBranchAddress("rpchit.fCell",rpchit_fCell);
   T->SetBranchAddress("rpchit.fCol",rpchit_fCol);
   T->SetBranchAddress("rpchit.fRow",rpchit_fRow);
   T->SetBranchAddress("rpchit.fX",rpchit_fX);
   T->SetBranchAddress("rpchit.fY",rpchit_fY);
   T->SetBranchAddress("rpchit.fZ",rpchit_fZ);
   T->SetBranchAddress("rpchit.fTime",rpchit_fTime);
   T->SetBranchAddress("rpchit.fCharge",rpchit_fCharge);
   //T->SetBranchAddress("tracks",&tracks_);
   T->SetBranchAddress("tracks.fUniqueID",tracks_fUniqueID);
   T->SetBranchAddress("tracks.fBits",tracks_fBits);
   T->SetBranchAddress("tracks.find0",tracks_find0);
   T->SetBranchAddress("tracks.find1",tracks_find1);
   T->SetBranchAddress("tracks.fAl",tracks_fAl);
   T->SetBranchAddress("tracks.fBe",tracks_fBe);
   T->SetBranchAddress("tracks.fGa",tracks_fGa);
   T->SetBranchAddress("tracks.fX",tracks_fX);
   T->SetBranchAddress("tracks.fY",tracks_fY);
   T->SetBranchAddress("tracks.fZ",tracks_fZ);
   T->SetBranchAddress("tracks.fTime",tracks_fTime);
   //T->SetBranchAddress("rpccalhitscorr",&rpccalhitscorr_);
   T->SetBranchAddress("rpccalhitscorr.fUniqueID",rpccalhitscorr_fUniqueID);
   T->SetBranchAddress("rpccalhitscorr.fBits",rpccalhitscorr_fBits);
   T->SetBranchAddress("rpccalhitscorr.fTrbnum",rpccalhitscorr_fTrbnum);
   T->SetBranchAddress("rpccalhitscorr.fCell",rpccalhitscorr_fCell);
   T->SetBranchAddress("rpccalhitscorr.fCol",rpccalhitscorr_fCol);
   T->SetBranchAddress("rpccalhitscorr.fRow",rpccalhitscorr_fRow);
   T->SetBranchAddress("rpccalhitscorr.fX",rpccalhitscorr_fX);
   T->SetBranchAddress("rpccalhitscorr.fY",rpccalhitscorr_fY);
   T->SetBranchAddress("rpccalhitscorr.fZ",rpccalhitscorr_fZ);
   T->SetBranchAddress("rpccalhitscorr.fTime",rpccalhitscorr_fTime);
   T->SetBranchAddress("rpccalhitscorr.fCharge",rpccalhitscorr_fCharge);
   */
//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
     // disable all branches
   T->SetBranchStatus("event",1);  // activate branchname
   T->SetBranchStatus("rpchit",1);  // activate branchname

   Int_t nentries = T->GetEntries();
   //nentries = 20000;
   //Long64_t nbytes = 0;
   cout<<"Events in the tree: "<<nentries<<endl;
   TH1D* hdt[10][12];
   for(Int_t i=0;i<10;i++) {
       for(Int_t j=0;j<12;j++) {
	   hdt[i][j] = new TH1D(Form("hdt_%i_%i",i,j),Form("hdt_%i_%i",i,j),1000,-20,20);
       }
   }

   TH2D* h_its  = new TH2D("hits","hits",100,-1200,1200,100,-1200,1200);
   TH2D* h_itse = new TH2D("hitse","hitse",100,-1200,1200,100,-1200,1200);


   TH2D* h_ith1  = new TH2D("hitsh1","hitsh1",100,-1200,1200,100,-1200,1200);
   TH2D* h_ith2  = new TH2D("hitsh2","hitsh2",100,-1200,1200,100,-1200,1200);
   TH2D* h_ith3  = new TH2D("hitsh3","hitsh3",100,-1200,1200,100,-1200,1200);



   for (Int_t i=0; i<nentries;i++) {
       if(i%100000==0)cout<<i<<endl;
       T->GetEntry(i);
       //cout<<event->getMultH1()<<" "<<event->getMultH2()<<" "<<event->getMultH3()<<endl;
       multH1 = event->getMultH1();
       multH2 = event->getMultH2();
       multH3 = event->getMultH3();
       //cout<<event->getRpcHits()<<endl;
       //
       //hits = (TClonesArray*)event->getRpcHits();
       //cout<<hits->GetEntriesFast()<<endl;
       Int_t nh =hits->GetEntriesFast();
       //cout<<multH1<<" "<<multH2<<endl;
       if(multH1 == multH2 && multH1==1 && multH3>0) {
	   //for(
	   //hdt[rpchit_fRow[1]][rpchit_fCol[1]]->Fill(rpchit_fTime[1]-rpchit_fTime[0]);
	   Int_t col, row;
	   Float_t tref, tm;
	   tref = -10000;
	   Bool_t valid = kFALSE;
           Float_t x,y,z,xref,yref,zref;
	   Float_t xref2,yref2,zref2;
	   Float_t xref3,yref3,zref3;
           Float_t xteo, yteo;
	   //cout<<"such event exists "<<endl;
           Int_t counter=0;
	   for(Int_t j=0;j<nh;j++) {
	       TRpcHit* rpchit = (TRpcHit*)hits->At(j);
	       if(rpchit->getTrbnum()==0) {
		   tref = rpchit->getTime();
                   zref = rpchit->getZ();
                   xref = rpchit->getX();
		   yref = rpchit->getY();

		   h_ith1->Fill(xref,yref);

	       }

	       if(rpchit->getTrbnum()==1) {
		   zref2 = rpchit->getZ();
                   xref2 = rpchit->getX();
                   yref2 = rpchit->getY();
                   h_ith2->Fill(xref2,yref2);
	       }



	       if(rpchit->getTrbnum()==2) {
                   counter++;
                   xref3 = rpchit->getX();
                   yref3 = rpchit->getY();
		   h_ith3->Fill(xref3,yref3);
	       }
	   }

	   xteo = (xref-xref2)/(zref-zref2) * (1826-zref) + xref;
	   yteo = (yref-yref2)/(zref-zref2) * (1826-zref) + yref;

	   h_its  -> Fill(xteo,yteo);

           Bool_t once = kTRUE;
	   if(counter>1)continue;
	   for(Int_t j=0;j<nh;j++) {
	       TRpcHit* rpchit = (TRpcHit*)hits->At(j);
	       if(rpchit->getTrbnum()==2) {
		   // Valid
                   // valid = kTRUE;
		   tm  = rpchit->getTime();
		   col = rpchit->getCol()-1;
		   row = rpchit->getRow()-1;

		   z = rpchit->getZ();
                   x = rpchit->getX();
                   y = rpchit->getY();

		   //cout<<"xteo "<<xteo<<" x "<<x<<" yteo "<<yteo<<" y "<<y<<endl;
		   if(once && fabs(xteo-x)<250 && fabs(yteo-y)<250) {
		       h_itse -> Fill(xteo,yteo);
                       once = kFALSE;
		   }
		   if(rpchit->getCharge()>40)
		       //hdt[row][col]->Fill(tm-tref+sqrt(pow(x-xref,2)+pow(y-yref,2)+pow(z-zref,2))/(2.99792458000000011e+02)-999.);
		       hdt[row][col]->Fill(tm-tref+sqrt(pow(x-xref,2)+pow(y-yref,2)+pow(z-zref,2))/(2.99792458000000011e+02));
	       }
	   }
           //if(valid)
       }


   }

   /*
   ifstream fin;
   fin.open("../pars/time_calPar.txt");
   Float_t pars[3][10][12][2];
   Int_t plane,row,col;
   Float_t val1,val2;
   while(!fin.eof()) {
       fin>>plane>>row>>col>>val1>>val2;
       pars[plane][row][col][0] = val1;
       pars[plane][row][col][1] = val2;
   }

   TF1* fg  = new TF1("fg","gaus",-50,50);

   for(Int_t i=0;i<10;i++) {
       for(Int_t j=0;j<12;j++) {
	   hdt[i][j]->Fit(fg,"NQ");
	   hdt[i][j]->Fit(fg,"NQR","",fg->GetParameter(1)-1.5,fg->GetParameter(2)+1.5);
           pars[2][i][j][1] = fg->GetParameter(1);

       }
   }
   */
   /*
   ofstream fout;
   fout.open("test.txt");
   for(Int_t d=0;d<3;d++) {
       for(Int_t i=0;i<10;i++) {
	   for(Int_t j=0;j<12;j++) {
	       fout<<d<<" "<<i<<" "<<j<<" "<<pars[d][i][j][0]<<" "<<pars[d][i][j][1]<<endl;
	       //for(Int_t p=0;p<2;p++)
	       //    fout<<pars[d][i][j][p];
               //fout<<endl;
	   }
       }
   }
   fout.close();
   */

   TFile* fOut = new TFile("fouttest2.root","recreate");
   /*
   for(Int_t i=0;i<10;i++) {
       for(Int_t j=0;j<12;j++) {
	   hdt[i][j]->Fit(fg,"NQ");
	   hdt[i][j]->Fit(fg,"NQR","",fg->GetParameter(1)-1.5,fg->GetParameter(2)+1.5);
       }
   }
   */
   h_its->Write();
   h_itse->Write();

   h_ith1 -> Write();
   h_ith2 -> Write();
   h_ith3 -> Write();



   for(Int_t i=0;i<10;i++) {
       for(Int_t j=0;j<12;j++) {
	   hdt[i][j]->Write();
       }
   }
   fOut->Close();

}
int main(int argc, char* argv[] ) {
    cout<<"Starting list "<<argv[1]<<endl;
    Analysis(argv[1]);
    return 1;

}