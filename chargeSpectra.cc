#include "chargeSpectra.h"
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
#include "ttrackf.h"
#include "ttrack.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include <fstream>
#include "TTimeStamp.h"
#include "TMath.h"
#include <algorithm>
#include <vector>
//#include <ofstream>


using namespace std;


//______________________________________________________________________________
Filler::Filler()
{
    EventNr=0;
    pEvent=0;
    pRootFile=0;
    subEvtId=0;
    fpga_code=0;
    refCh = -1;

    //TH2F* h_q[3][10][12];
    //TH2F* h_q_pedestal[3][10][12];

    for(Int_t i=0;i<3;i++) {
        for(Int_t j=0;j<10;j++) {
            for(Int_t k=0;k<12;k++) {
                h_q[i][j][k] = new TH2F(Form("h_q_det%i_row%i_column%i",i,j,k),Form("h_q_det%i_row%i_column%i",i,j,k),140*12,60,200,502,-1,250);
            }
        }
    }

    h_daq_active     = new TH1D("h_daq","h_daq;D.O.Y.;Frequency",366*24*6,0.,366.);
    h_daq_sync       = new TH1D("h_daq_sync","h_daq_sync;D.O.Y.;Frequency",366*24*6,0.,366.);
    h_daq_seconds    = new TH1D("h_daq_seconds","h_daq_seconds;D.O.Y.;Frequency",366*24*6,0.,366.);
    h_daq_total      = new TH1D("h_daq_total","h_daq_total;D.O.Y.;Frequency",366*24*6,0.,366.);


    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

}

void Filler::fillHistograms(const char* dir, TString list, const char* odir,const char* ofile, Int_t nEvt,Int_t n)
{
//    TFile* fOut = new TFile(ofile,"RECREATE");
    TFile fOut(Form("%s%s",odir,ofile),"RECREATE");
    refCh = -1;
    refCh = 31;
    this->subEvtId= 0;
    cout<<"ofile "<<outputFile<<endl;
    std::ifstream file(list.Data());
    TString fName;
    while(!file.eof()) {
        file>>fName;
        cout<<" file  "<<fName<<" directory "<<dir<<endl;
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
            eventLoopFill(nEvt,0);
        }
        delete pEvent;
        pEvent = NULL;
    }

    file.close();
    fOut.cd();

    // Do your staff to the histograms




    h_daq_active  -> Write();
    h_daq_sync    -> Write();
    h_daq_seconds -> Write();
    h_daq_total   -> Write();


    for(Int_t i=0;i<3;i++) {
        for(Int_t j=0;j<10;j++) {
            for(Int_t k=0;k<12;k++) {
                h_q[i][j][k] -> Write();
            }
        }
    }



    fOut.Close();




}

//______________________________________________________________________________
Filler::~Filler()
{
    if(pEvent)	delete pEvent;
    if(pRootFile)	delete pRootFile;

    // Delete all your histograms!

}

//______________________________________________________________________________
Bool_t Filler::setRootFile(const char* filename/*="test.root" */)
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
string Filler::setInputFile(const char* dir,const char* filename)
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

string Filler::setInputFile(const char* filename)
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
Bool_t Filler::setpEvent(Int_t Id)
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
Bool_t Filler::eventLoopFill(Int_t nbEvt, Int_t startEv)
//Loop over all events, data written to the root tree
{
    if(pEvent==0)
    {
        cout<<"Error: no pEvent set"<<endl;
        return kFALSE;
    }
    else
    {
        cout<<"loading lookup "<<fileLookupPar<<endl;




        Int_t trbnum, cell, col, row;
        Float_t     x, y, z;
        Float_t time, charge;

        Int_t trbnum2, cell2, col2, row2;
        Float_t     x2, y2, z2;
        Float_t time2, charge2;

        // Will execute the unpacker
        Event* event = new Event();

        // Will execute the calibrater
        TRpcRawF* rawFinder = new TRpcRawF();
        rawFinder -> init(fileLookupPar);

        RpcLookupTable* fLookPar = new RpcLookupTable(fileLookupPar);

        /*
        // Will execute the hit finder
        TRpcHitF* calFinder = new TRpcHitF();
        //calFinder -> init(fileHitFinderPar,"/media/Datos2TB/korna/tragaldabas/pars/active_cells_2015_2016.root");
        calFinder -> init(fileHitFinderPar,"/media/Datos2TB/korna/tragaldabas/soft2/GoodActiveCells_2015_2016.root");

        // Will execute the track finder and
        // the three planes track finder!
        TTrackF* trackFinder = new TTrackF();
        trackFinder->init();
        */
        // Pointers to the clonesarray of every category of data
        TClonesArray*  rpchits          = rawFinder->getRpcRawHits();
        //TClonesArray*  rpccalhits       = calFinder->getRpcHits();
        //TClonesArray*  tracks           = trackFinder->getTracks();
        //TClonesArray*  tracks3planes    = trackFinder->getTracks3Planes();
        //TClonesArray*  rpccalhitscorr   = trackFinder->getRpcHitCorr();

        TClonesArray* hits = event->getHits();

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
            if(!event)
                event = new Event(*pEvent, refCh);
            event->clearAll();
            event->fill(*pEvent);







            event        -> setSync(sync);
            rawFinder    -> execute();
            sync = event -> getSync();
            //calFinder    -> execute();
            //trackFinder  -> execute();

            Double_t timeEvt = TTimeStamp::GetDayOfYear(event->getEvtDay(),event->getEvtMonth(),event->getEvtYear());
            timeEvt+= event->getEvtHour()/24.;
            timeEvt+= event->getEvtMinute()/24./60.;
            timeEvt+= event->getEvtSecond()/24./60./60.;



            h_daq_total    -> Fill(timeEvt);

            // histogram with all events.
            // Then the unsyncronized will be obtained by subtraction

            // Do not include event out of sync!
            if(!sync) {
                h_daq_sync->Fill(timeEvt);
                continue;
            }

            TClonesArray*  fTrbHits = gEvent->getHits();
            Int_t ntrb = fTrbHits->GetEntriesFast();
            Bool_t sinc1=0;
            Bool_t sinc2=0;
            Bool_t sinc3=0;
            for(Int_t i=0;i<ntrb;i++) {
                Hit* hittrb = (Hit*)fTrbHits->At(i);
                if(!hittrb) continue;
                Int_t trb    = hittrb->getTRBNum();
                Int_t tdc    = hittrb->getTDC();
                Int_t chan   = hittrb->getChannel();
                Int_t trbnum = fLookPar->getDetNum(trb);
                if(trbnum==-1)continue;
                //Float_t time     = hittrb->getLeadTime1();
                //Float_t charge     = hittrb->getWidth();
                //fLookPar -> getParams(trbnum, tdc, chan, cell, col, raw, x, y, z);
                if(trbnum==0&&chan==30&&tdc==1)sinc1=1;
                if(trbnum==1&&chan==30&&tdc==1)sinc2=1;
                if(trbnum==2&&chan==30&&tdc==1)sinc3=1;
            }
            //  Get the DOY with time in decimal units!
            // this is the "number of seconds active"

            //if(prevsec==-100) {
            //    prevsec = event->getEvtSecond();
            //    h_daq_active->Fill(timeEvt);
            //}

            if(prevsec != event->getEvtSecond() ) {
                prevsec = event->getEvtSecond();
                h_daq_seconds  -> Fill(timeEvt);
                //h_daq_active->Fill(timeEvt);
            }
            //if(sinc1 &&sinc2 && sinc3) {
            //    h_daq_active->Fill(timeEvt);
            //}
            if(event->getEvtYear() == 2017 ) {
                if(sinc1 &&sinc2) {
                    h_daq_active->Fill(timeEvt);
                }
            } else {
                if(sinc1 &&sinc2 && sinc3) {
                    h_daq_active->Fill(timeEvt);
                }
            }

            // For loop over TRB hits
            //
            //


            // For loop over cal hits
            //
            //

            // For loop over hits!
            rpchits = rawFinder->getRpcRawHits();
            Int_t nrpchits = rpchits->GetEntriesFast();

            for(Int_t j=0;j<nrpchits;j++) {
                TRpcRaw* rpchit = (TRpcRaw*)rpchits->At(j);
                if(rpchit) {
                    rpchit->getHit(trbnum,cell,col,row,x,y,z,time,charge);
                    //cout<<"trbnum "<<trbnum<<" "<<row<<" "<<col<<endl;
                    h_q[trbnum][row-1][col-1] -> Fill(timeEvt,charge*0.098);
                }
            }
        }
        delete event;
        delete rawFinder;
        return kTRUE;
    }
}

void doStaffAnalysis(char* path,char* name){
    TString fName(name);
    TString fPath(path);
    Filler fill =  Filler();
    TString day(fName(11,fName.First(".")-11));
    cout<<"DAY "<<day<<endl;
    cout<<"PATH to files "<<"../"+fName<<endl;
    fill.setFileLookupPar("../pars/luptable_corr_20180423.txt");
    fill.setFileHitFinderPar("../2018DST/pars/"+day+"_CalPars.txt");
    fill.fillHistograms(fPath,"../"+fName,"../2018DSTqstudy/results/",day+"_result_histos.root",1000000,0);
}


int main(int argc, char* argv[] ) {
    cout<<"Starting list "<<argv[2]<<" data path: "<<argv[1]<<endl;
    doStaffAnalysis(argv[1],argv[2]);
    return 1;

}










