#include "thit.h"
#include "trpcrawf.h"


ClassImp(TRpcRawF);

TRpcRawF::TRpcRawF() {
    fLookPar    = new RpcLookupTable();
    fTrbHits    = getTrbRawHits();
    //if(!fRpcRawHits)
    fRpcRawHits = new TClonesArray("TRpcRaw",1000);
}

Int_t TRpcRawF::init(TString filename) {
    fLookPar    = new RpcLookupTable(filename);
    fTrbHits    = getTrbRawHits();
    //if(!fRpcRawHits)
    fRpcRawHits = new TClonesArray("TRpcRaw",1000);
    return 1;
}
TRpcRawF::~TRpcRawF() {
    delete fLookPar;
    delete fRpcRawHits;

}
Int_t TRpcRawF::execute() {
    // Clear the TClones array for later saving it!
    // clearAll();
    totalNHits = 0;
    //cout<<"point fRpcRawHits "<<fRpcRawHits<<endl;
    fRpcRawHits->Clear("C");
    Int_t cell, col, raw;
    Float_t     x, y, z;
    // Loop through the TrbHits
    // Int_t ntrb = fRpcRawHits->GetEntriesFast();
    Int_t ntrb = fTrbHits->GetEntriesFast();

    Int_t sinc1=0;
    Int_t sinc2=0;
    Int_t sinc3=0;

    // We have to read the number of TRBs in the setup!

    Int_t numTrbActive[3] = {0,0,0};

    //cout<<"new event "<<endl;
    for(Int_t i=0;i<ntrb;i++) {
	Hit* hittrb = (Hit*)fTrbHits->At(i);
	if(!hittrb) continue;
	Int_t trb    = hittrb->getTRBNum();
	Int_t tdc    = hittrb->getTDC();
	Int_t chan   = hittrb->getChannel();
	Int_t trbnum = fLookPar->getDetNum(trb);
	if(trbnum==-1)continue;
	Float_t time     = hittrb->getLeadTime1();
	Float_t charge     = hittrb->getWidth();
	fLookPar -> getParams(trbnum, tdc, chan, cell, col, raw, x, y, z);
        //if(chan==30)
	//    cout<<"cell "<<cell<<" trb "<<trb<<" tdc "<<tdc<<" chan "<<chan<<" trbn "<<trbnum<<endl;
	if(cell!=-1) {
	    TRpcRaw* fRaw = addRpcRaw();
            fRaw->setHit( trbnum, cell, col, raw, x, y, z, time, charge);
   	    //cout <<trbnum<< " " <<raw<< " " <<col<< " " <<time<<" "<<charge<<" " << endl;
        }
        if(chan!=31 && cell==-1) {
            if(chan==30 && tdc==1) {
            }
            else {
                //cout<<"cell "<<cell<<" trb "<<trb<<" tdc "<<tdc<<" chan "<<chan<<" trbn "<<trbnum<<" time "<<time<<" ch "<<charge<<endl;
            }
        }
        if(trbnum==0)numTrbActive[0]=1;
        if(trbnum==1)numTrbActive[1]=1;
        if(trbnum==2)numTrbActive[2]=1;

	if(trbnum==0&&chan==30&&tdc==1)sinc1=1;
        if(trbnum==1&&chan==30&&tdc==1)sinc2=1;
        if(trbnum==2&&chan==30&&tdc==1)sinc3=1;


    }








    Int_t trbActive = numTrbActive[0]+numTrbActive[1]+numTrbActive[2];

    if(sinc1==0&&gEvent->getSync())          gEvent->setSync(1);
    if(sinc1==0&&!gEvent->getSync())         gEvent->setSync(0);
    if((sinc1==1||sinc2==1 || sinc3==1)&& (trbActive==(sinc1+sinc2+sinc3)) ) gEvent->setSync(1);
    if(sinc1!=sinc2 || sinc1!=sinc3 || sinc2!=sinc3) gEvent->setSync(0);

    //cout<<gEvent->getSync()<<endl;
    gEvent->setRpcRawHits(fRpcRawHits);
    //setRpcHits();

    return 1;

}
TRpcRaw* TRpcRawF::addRpcRaw( ) {
    TClonesArray& hits = *fRpcRawHits;
    TRpcRaw *rpcraw = new (hits[totalNHits++]) TRpcRaw();
    //cout<<"ENTRIES in addhit "<<Hits->GetEntriesFast()<<endl;
    return rpcraw;
}

