#include "thit.h"
#include "tevent.h"
#include "thldevent.h"

using namespace std;

ClassImp(Event)

TClonesArray *Event::gHits = 0;
R__EXTERN Event *gEvent;


Event::Event() : kMaxMult(10), kMaxChannelNr(256)
{
    // Create an Event object.
    // When the constructor is invoked for the first time, the class static
    // variable gHits is 0 and the TClonesArray gHits is created.
    gEvent = this;
    if (!gHits)
	gHits = new TClonesArray("Hit", 1000);         // Me da igual que se haga la comprobacion
    Hits = gHits;

    //HitsRaw = new TClonesArray("TObject",1000);
    clearAll();

}
Event::Event(const HldEvent& HldEvt, Int_t refCh) : kMaxMult(HldEvt.getkMaxMult()), kMaxChannelNr(HldEvt.getkMaxChannelNr())
{
    // NADA DE ESTO DEBERIA ESTAR AQUI...
    // TODOS ESTOS ELEMENTOS DEBEN INICIALIZARSE DESDE LAS CLASES CORRESPONDIENTES
    // COMO RPCUNPACK
    // COMO RPCCALIBRATER
    // TRACK FINDER
    // Y ASI.
    // esto se inicializa solo una vez.
    // Despues para cada evento empezamos.
    // El clear si que hace falta con cada evento pero nada mas!

    // setReferenceChannel(refCh);
    if (!gHits)
	gHits = new TClonesArray("Hit", 1000);
    Hits = gHits;
    clearAll();
    //HitsRaw = new TClonesArray("TObject",1000);
    fill(HldEvt);
}

Event::~Event() {
    //delete [] Hits;
    //delete [] gHits;
    //delete  Hits;     //! * array with all hits
    //delete  rpcRaw;   //! * array with rpcRaw
    //delete  rpcHit;   //! * array with rpcHits
    clearAll();
}

void Event::clearAll() {
    //Hits->Clear();
    totalNHits = 0;
    errors_per_event = 0;

    multT = 0.;
    multT3 = 0.;


    //setReferenceChannel(127);  -- gk 18.05.2011
    setReferenceTime(-100000);
    EvtHdr.clearAll();
}

Bool_t Event::fill(const HldEvent& HldEvt) {
    //Hit* pCurrentHit = new Hit();
    EvtHdr.fill(HldEvt);       // ESTO SI QUE LO USARE
    setSubEvtId(HldEvt.getSubEvtId());
    /*
    Int_t leadTime;
    Int_t widTime;
    Int_t trailTime;
    Int_t leadMult;
    
    bool multihitReference = false;

    for (Int_t i = 0; i < kMaxChannelNr; i++) {
	leadMult = HldEvt.getLeadingMult(i);

	if (leadMult < 1)
	    continue; //Leading Time data for this channel does not exist
	pCurrentHit =  addHit();
	pCurrentHit -> setChannel(i);
	pCurrentHit -> setTDC(i);

	pCurrentHit->setNHits(leadMult);
	// gk check if there was more hits on reference channel
	if (leadMult > 1 && i == getReferenceChannel()) {
	    multihitReference = true;
	}
	// fill time info for channel: TDC, channel
	for (Int_t chmult = 0; (chmult < leadMult && chmult < 4); chmult++)
	{
	    leadTime  = HldEvt.getLeadingTime(i, chmult);
	    trailTime = HldEvt.getTrailingTime(i, chmult);
	    widTime   = HldEvt.getWidthTime(i, chmult);
	    pCurrentHit->fillTimeAndWidth(leadTime, widTime);
	    if (chmult == 0)
	    {
		pCurrentHit->setLeadTime1(leadTime);
		pCurrentHit->setTrailTime1(trailTime);
	    }
	    else if (chmult == 1)
	    {
		pCurrentHit->setLeadTime2(leadTime);
		pCurrentHit->setTrailTime2(trailTime);
	    }
	    else if (chmult == 2)
	    {
		pCurrentHit->setLeadTime3(leadTime);
		pCurrentHit->setTrailTime3(trailTime);
	    }
	    else if (chmult == 3)
	    {
		pCurrentHit->setLeadTime4(leadTime);
		pCurrentHit->setTrailTime4(trailTime);
	    }
	    pCurrentHit->setWidth(leadTime, trailTime);
	} //for
    } // for(Int_t i=0; i<128; i++)

    if(multihitReference == true) {
	int t = 128 * (getReferenceChannel() / 128) - 1;  // get the initial channel of trb with selected reference channel
	Int_t refTimesMean = 0;
	int mults[4] = {-1, -1, -1, -1};
	Int_t vals[4] = {0, 0, 0, 0};
	Int_t temp = 0;

	for(Int_t i = 0; i < 4; i++) {  // iterate through multiplicity of first ref channel
	    mults[0] = i;
	    mults[1] = -1;
	    mults[2] = -1;
	    mults[3] = -1;
	    vals[0] = HldEvt.getLeadingTime(t + 32, i);

	    for(Int_t j = 0; j < 4; j++) {
		temp = HldEvt.getLeadingTime(t + 64, j);
		if (temp < vals[0] + 200 && temp > vals[0] - 200) {
		    mults[1] = j;
		    vals[1] = temp;
		}
	    }

	    for(Int_t j = 0; j < 4; j++) {
		temp = HldEvt.getLeadingTime(t + 96, j);
		if (temp < vals[0] + 200 && temp > vals[0] - 200) {
		    mults[2] = j;
		    vals[2] = temp;
		}
	    }

	    for(Int_t j = 0; j < 4; j++) {
		temp = HldEvt.getLeadingTime(t + 128, j);
		if (temp < vals[0] + 200 && temp > vals[0] - 200) {
		    mults[3] = j;
		    vals[3] = temp;
		}
	    }

	    if (mults[1] != -1 && mults[2] != -1 && mults[3] != -1) {
		break;
	    }
	}
	//     cerr<<"Found multihit on reference channel, selecting hit "<<mults[getReferenceChannel() / 128]<<endl;
	setReferenceTime(HldEvt.getLeadingTime(getReferenceChannel(), mults[getReferenceChannel() / 128]));
    }
    else {
	setReferenceTime(HldEvt.getLeadingTime(getReferenceChannel(), 0));
    }
    */ // TODO ESTO LO MOVI MIENTRAS QUE SE LEE EL EVENTO


    setErrors_per_event(HldEvt.errors_per_event);
    return kTRUE;
}
Float_t Event::getDOY() {
    time.Set(EvtHdr.getYear()+1900,EvtHdr.getMonth(),EvtHdr.getDay(),0.,0.,1.,0,kTRUE,0);
    Float_t doy = time.GetDayOfYear() + EvtHdr.getHour()/24. + EvtHdr.getMinute()/24./60.;
    return doy;
}

//______________________________________________________________________________
Hit *Event::addHit()
{
    TClonesArray& hits = *Hits;
    Hit *hit = new (hits[totalNHits++]) Hit();
    return hit;
}
Event *gEvent;
