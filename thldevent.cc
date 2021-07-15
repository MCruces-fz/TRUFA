
#include "thldevent.h"
#include "TError.h"
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <math.h>



R__EXTERN Event* gEvent;

using namespace std;

ClassImp(HldEvent)


Bool_t HldEvent::setFile(const Text_t *name)
{
    if (access(name, F_OK) != 0)
    {
	Error("HldEvent::setFile()", "Error during access of file %s !", name);
	exit(EXIT_FAILURE);
    }

    stat(name, &status);
    if (file)
    {
	((ifstream*) file)->open(name);
	if (!((ifstream*) file)->is_open())
	{
	    Error("HldEvent::setFile()", "Error during opening file %s !", name);
	    exit(EXIT_FAILURE);
	}
	if (!((ifstream*) file)->good())
	{
	    Error("Hl    dEvent::setFile()",
		  "Error flags discovered in stream of file %s !", name);
	    exit(EXIT_FAILURE);
	}

    }
    else
    {
	file = new std::ifstream(name);
	if (!((ifstream*) file)->is_open())
	{
	    Error("HldEvent::setFile()", "Error during opening file %s !", name);
	    exit(EXIT_FAILURE);
	}
	if (!((ifstream*) file)->good())
	{
	    Error("HldEvent::setFile()",
		  "Error flags discovered in stream of file %s !", name);
	    exit(EXIT_FAILURE);
	}
    }
    return kTRUE;
}

Bool_t HldEvent::read()
{
    if (file)
    {
	if (file->eof())
	{
	    return kFALSE;
	}
	if (!((ifstream*) file)->good())
	{
	    Error("HldEvent::read()",
		  "Error flags discovered in input stream before read!");
	    exit(1);
	}

	if (pData)
	{
	    delete[] pData;
	    pData = 0;
	}
	if (file->read((char *) (pHdr), getHdrSize()))
	{
	    if (isSwapped())
		swapHdr();
	    if (getSize() > getHdrSize())
	    {

#warning "Hardcoded maximum event size is 1000000"
		if (getDataLen() > 1000000)
		    return kFALSE;
		pData = new UInt4[getDataLen()];
		if (file->read((char*) (pData), getDataSize()))
		    file->ignore(getPaddedSize() - getSize());
	    }
	}

	if (!((ifstream*) file)->good())
	{
	    if (!file->eof())
	    {
		Error("HldEvent::read()",
		      "Error flags discovered in input stream after read!");
	    }
	    return kFALSE;
	}
	else
	{
	    return kTRUE;
	}

    }
    return kFALSE;
}

Bool_t HldEvent::readSubEvt(size_t i)
{
    UInt4* p;
    if (i)
	p = subEvt[i - 1].getPaddedEnd();
    else
	p = pData;
    if (p < getPaddedEnd())
	subEvt[i] = HldSubEvent(p);
    else
	return kFALSE;
    return kTRUE;
}

Bool_t HldEvent::execute()
{
    //cout<<" total hits previous "<<totalNHits<<endl;
    totalNHits = 0;
    TClonesArray* hits = gEvent->getHits();
    hits->Clear();
    //cout<<"ENTRIES in execute "<<hits->GetEntriesFast()<<endl;
    clearAll();
    if (read())
    {
	for (size_t idx = 0; idx < lastSubEvtIdx; idx++)
	    *subEvtTable[idx].p = 0;

	for (size_t i = 0; i < maxSubEvts && readSubEvt(i); i++)
	{ // loop subevts
	    Bool_t unpacked = kFALSE;

	    for (size_t idx = 0; idx < lastSubEvtIdx; idx++)
	    { // loop unpackers
		if (kTRUE)
		{
		    subEvt[i].swapData();
		    *subEvtTable[idx].p = subEvt + i;
		    //        cout << "pass a pointer to subevt to its unpacker" << endl;
		    unpacked = kTRUE;
		}
	    }
	    if (isWritable && !unpacked)
	    {
		// to assure that the swapData() can work
		if (((subEvt[i].getWordSize()) == 1)
		    || ((subEvt[i].getWordSize()) == 2)
		    || ((subEvt[i].getWordSize()) == 4))
		{
		    subEvt[i].swapData();
		}
		else
		    Warning("execute", "Corrupted SubEvent, SubId %x",
			    subEvt[i].getId());
	    }
            // Decode the subevent with the corresponding id of the TRB!
	    if (pSubEvt)
	    {
		// Call RPC unpacker!
		// In the RPC unpacker we should get access to the "rpccal" objects!
                // the HLD event should have all the data already unpacked!
		Int_t res = decode();
                //cout<<"decode out "<<res<<" correcction reftime "<<endl;
		correctRefTimeCh31();
                //cout<<" correcction hecha! "<<endl;
                // store them to the  trbraw data category
                fillTrbRaw(subEvt[i].getId());
		gEvent->setNHits(totalNHits);


	    }
	}  
	return kTRUE;
    }
    else
	return kFALSE;

}
void HldEvent::fillTrbRaw(Int_t trbnum)  {

    TClonesArray* hits = gEvent->getHits();
    //cout<<"GETTING THE POINTER TO HITS "<<hits<<endl;
    //if(!hits)
    //    hits = new TClonesArray("Hit",1000);
    Bool_t flag      = kFALSE;
    Int_t refCh      = -1;
    Int_t jo         = -1;
    Float_t corrTime = 0.;
    for(Int_t jj=0; jj<4; jj++) {

	//// First check if there is ANY data on a TDC, if not, skip it
	//flag  = kFALSE;
	jo    = jj*32;
	//for(Int_t ll=0; ll<32; ll++) {
	//    if( trbLeadingMult[jo+ll]>0){
	//	flag = kTRUE;
	//	break;
	//    }
	//}
	//cout<<" entre 0 "<<endl;

	//if(!flag) continue;
	// ELSE do the correction
	//else{
	    for(Int_t ll=0; ll<32; ll++) {  // For all TDC channels
		Int_t ii;
		ii=jo+ll;
		if (trbLeadingMult[ii]<1 && trbTrailingTotalMult[ii]<1) continue;
		Hit* hit = addHit(hits);
		Int_t max = trbTrailingTotalMult[ii];
		if(max>=10) {
		    //cout<<"TOO MANY TRAILINGS!! trbnum "<<trbnum<<" mbo "<<jj<<" chan "<<ll<<endl;
		    continue;
		}
		// Store first leading and last traling
		hit->fillTimeAndWidthAndTrb(trbLeadingTime[ii][0],trbADC[ii][0],trbnum);

		//cout<<" chan "<<ll
		//    <<" TDC " <<jj
		//    <<" nhits "<<trbLeadingMult[ii]
		//    <<" trbLeadingTime[ii][0] "<<trbLeadingTime[ii][0]
		//    <<" trbTrilTime[ii][0] "<<trbTrailingTime[ii][max-1]
		//    <<endl;

		hit->setChannel(ll);
		hit->setTDC(jj);
		hit->setNHits(trbLeadingMult[ii]);

		if(trbTrailingTime[ii][0]>-100)  {
		    hit->setTrailTime1(trbTrailingTime[ii][0]);
		}
		if(trbTrailingTime[ii][1]>-100)  {
		    hit->setTrailTime2(trbTrailingTime[ii][1]);
		    hit->setLeadTime2(trbLeadingTime[ii][1]);

		}
		if(trbTrailingTime[ii][2]>-100)  {
		    hit->setTrailTime3(trbTrailingTime[ii][2]);
		    hit->setLeadTime3(trbLeadingTime[ii][2]);
		}
		if(max>0 && trbTrailingTime[ii][max-1]>-100)  {
		    hit->setTrailTime4(trbTrailingTime[ii][max-1]);
		    hit->setLeadTime4(trbLeadingTime[ii][max-1]);
		}

                //cout<<" entre "<<endl;

	    }
	//}
    }
    //cout<<"ENTRIES in fillTRB "<<hits->GetEntriesFast()<<endl;
    return;
}

Hit *HldEvent::addHit(TClonesArray* Hits)
{
    TClonesArray& hits = *Hits;
    Hit *hit = new (hits[totalNHits++]) Hit();
    //cout<<"ENTRIES in addhit "<<Hits->GetEntriesFast()<<endl;
    return hit;
}

//______________________________________________________________________________
Bool_t HldEvent::swap()
{
    //only swapping correctly the header
    if (read())
    {
	return kTRUE;
    }
    else
    {
	return kFALSE;
    }
}

//from HldBase
//______________________________________________________________________________
void HldEvent::swap2(UShort_t* ptr, const size_t len) const
{
    for (size_t i = 0; i < len; ++i)
    {
	UShort_t tmp = ptr[i];
	UChar_t* t = (UChar_t*) &tmp;
	UChar_t* p = (UChar_t*) &ptr[i];
	p[0] = t[1];
	p[1] = t[0];
    }
}

//______________________________________________________________________________
void HldEvent::swap4(UInt4* p, const size_t len) const
{
    UInt4 tmp;
    UInt1* d;
    UInt1* const s = (UInt1*) &tmp;
    for (size_t i = 0; i < len; i++)
    {
	d = (UInt1*) (p + i);
	tmp = p[i];
	d[0] = s[3];
	d[1] = s[2];
	d[2] = s[1];
	d[3] = s[0];
    }
}

//end HldBase

//from HldEvt
//______________________________________________________________________________
Int_t HldEvent::appendSubEvtIdx()
{
    subEvtTable[lastSubEvtIdx].id = getSubEvtId();
    //cout<<"subevtID "<< getSubEvtId()   <<endl;
    subEvtTable[lastSubEvtIdx].p = getpSubEvt();
    if (lastSubEvtIdx == maxSubEvtIdx - 1)
    {
	printf("\nMax. nb of unpackers (%d) exceeded!\n\n", maxSubEvtIdx);
	return 0;
    }
    else
	return ++lastSubEvtIdx;
}

//end HldEvt

void HldEvent::clearAll(void)
{

    //totalNHits = 0;  // ONLY ONCE PER EVENT
    for (Int_t i = 0; i <  128; i++)
    {
	for (Int_t k = 0; k < kMaxMult; k++)
	{
	    LeadingTime[i][k] = -1000000;
	    WidthTime[i][k] = -1000000;
	    TrailingTime[i][k] = -1000000;
	    trbADC[i][k] = -1;
	}
	LeadingMult[i] = 0;
	WidthMult[i] = 0;
	TrailingMult[i] = 0;
	SpikesCtr[i] = 0;
    }
    errors_per_event = 0;

    for(Int_t i=0; i<128; i++){
	for(Int_t k=0; k<10; k++){

	    trbLeadingTime[i][k]     = -1000000;
	    trbTrailingTime[i][k]    = -1000000;
	    trbADC[i][k]           = -1;

	}

	trbDataExtension[i]     = 0;
	trbLeadingMult[i]       = 0;
	trbTrailingMult[i]      = 0;
	trbTrailingTotalMult[i] = 0;
    }

    highResModeOn    = kFALSE;
    correctINLboard  = kFALSE;




}

Int_t HldEvent::correctRefTimeCh31(void) {

  ///////////////////////////////////////////
  // Reference signal from channel 31 of TDC
  // is used for correction of corresponding TDC chns
  // this is not the final version of
  // hardware design
  // call this function only if channel 31 contains
  // reference time
  ///////////////////////////////////////////

  Bool_t flag      = kFALSE;
  Int_t refCh      = -1; 
  Int_t jo         = -1;
  Float_t corrTime = 0.;


  // scan all arrays which contain time data and do correction
  for(Int_t jj=0; jj<4; jj++) {      

    // First check if there is ANY data on a TDC, if not, skip it
    flag  = kFALSE;
    jo    = jj*32;
    for(Int_t ll=0; ll<32; ll++) {  
      if( trbLeadingMult[jo+ll]>0){
	flag = kTRUE;
	break;
      }
    }
    if(!flag) continue;

    //If there is data, check reference channel
    refCh    = jo + 31;
    corrTime = trbLeadingTime[refCh][0];


    // If NO reference time -> set all to default values
    if(corrTime <= -1000000 ){
      //if(!quietMode ) {
	//if( !gHades->isCalibration() )Info("correctRefTimeCh","No Ref Time! SubEvtID 0x%x, TRBid 0x%x",subEvtId,uTrbNetAdress);
      //}
      clearAll();
    }
    
    // ELSE do the correction 
    else{
      for(Int_t ll=0; ll<32; ll++) {  // For all TDC channels
	Int_t ii;
	ii=jo+ll;
	for(Int_t kk=0; kk<10; kk++) {
	  trbLeadingTime[ii][kk]  = trbLeadingTime[ii][kk]  - corrTime + 40000;// We want poitive times here; 20000 ~ 2us
	  trbTrailingTime[ii][kk] = trbTrailingTime[ii][kk] - corrTime + 40000;// thats bigger than max time window of TDC
	}
      }
    }
  }
  
  return 0;
}

//______________________________________________________________________________
Bool_t HldEvent::fill_lead(Int_t ch, Int_t time)
{
    Int_t leadMult = LeadingMult[ch];
    Int_t trailMult = TrailingMult[ch];
    if(leadMult - trailMult == 0) {
	LeadingMult[ch]++;
	if (leadMult < kMaxMult)  {
	    LeadingTime[ch][leadMult] = time;
	}
	return ((leadMult + 1) <= kMaxMult);
    }
    else {
	LeadingTime[ch][leadMult - 1] = time;
	SpikesCtr[ch]++;
    }
}

//______________________________________________________________________________
Bool_t HldEvent::fill_width(Int_t ch, Int_t time)
{

    //hk added
    Int_t widMult = WidthMult[ch]; //width Multiplicity
    Int_t leadMult = LeadingMult[ch]; //Leading Multiplicity

    WidthMult[ch]++;

    if (widMult < kMaxMult)
    {
	if (widMult <= leadMult + 1)
	{
	    WidthTime[ch][widMult] = time;
	}
	else
	{
	    return kFALSE;
	}
    }

    return ((widMult + 1) <= kMaxMult);
}

//______________________________________________________________________________
Bool_t HldEvent::fill_trail(Int_t ch, Int_t time)
{

    Int_t trailMult = TrailingMult[ch]; //Trailing Multiplicity
    Int_t leadMult = LeadingMult[ch]; //Leading Multiplicity

    if (trailMult < kMaxMult)
    {
	if (leadMult - trailMult == 1)
	{
	    TrailingMult[ch]++;
	    TrailingTime[ch][trailMult] = time;
	}
	else
	{
	    SpikesCtr[ch]++;

	    return kFALSE;
	}
    }

    return ((trailMult + 1) <= kMaxMult);

}

//______________________________________________________________________________
void HldEvent::PrintTdcError(UInt_t e)
{
    Char_t *e_str[15] =
    {
	"Hit lost in group 0 from read-out FIFO overflow",
	"Hit lost in group 0 from L1 buffer overflow",
	"Hit error have been detected in group 0",
	"Hit lost in group 1 from read-out FIFO overflow",
	"Hit lost in group 1 from L1 buffer overflow",
	"Hit error have been detected in group 1",
	"Hit lost in group 2 from read-out FIFO overflow",
	"Hit lost in group 2 from L1 buffer overflow",
	"Hit error have been detected in group 2",
	"Hit lost in group 3 from read-out FIFO overflow",
	"Hit lost in group 3 from L1 buffer overflow",
	"Hit error have been detected in group 3",
	"Hits rejected because of programmed event size limit",
	"Event lost (trigger FIFO overflow)",
	"Internal fatal chip error has been detected"
    };

    if (e == 0)
	return;// No Error

    cout << "=== TRB/TDC Error analysis:" << endl;
    for (Int_t i = 0; i < 15; i++)
    {
	if (e & 0x1)
	{
	    cout << e_str[i] << endl;
	}
	e >>= 1;
    }
    cout << "===" << endl;
}

//______________________________________________________________________________
Int_t HldEvent::decode(Int_t num)
{
    clearAll();

    Int_t TdcId;
    Int_t nCountTDC = 0;
    UInt_t nEvtNr;

    Int_t minWindow = -100000;
    Int_t maxWindow = 100000;

    UInt_t nSizeCounter = 0;
    UInt_t nTdcEvtId = 0;
    UInt_t TdcDataLen = 0;
    UInt_t uBlockSize = 0;
    UInt_t subeventFound = false;
    UInt_t* data = pSubEvt->getData();
    UInt_t* end = pSubEvt->getEnd();

    Int_t currentFpga = -1;


    //data++;
    //test to jum to last TrbHeader dataworld
    data+=2;

    if((*data & 0xFFFF0000) == 0xBE010000)
    {
	return (kFALSE);
    }

    //  uBlockSize = *data & 0xFF;
    uBlockSize = (*data>>16) & 0xFFFF;

    //cout<<" blocksize "<< uBlockSize;
    //uBlockSize = (*(data+1)) & 0xFFFF;
    //cout<<" blocksize2 "<< uBlockSize<<endl;
    if(uBlockSize > (UInt_t)(end-data)) cout<<" block size different from end-data "<<endl;

    nEvtNr = (*data >> 8) & 0xFF;

    nSizeCounter++;// First one already processed


    bool foundLeadingEdge = false;

    bool printDebug = false;
    //bool printDebug = true;


    //data+=2;
    data++;



    while (data < end)
    {
	UInt_t dataword;
	dataword = *data;//[ii];
	nSizeCounter++;
	UInt_t* endOfSubevent;

        /*
	if (fullSetup == true) {
	    // gk 09.12.11
	    //gk find if the subevent is the one from endpoints list, skip others
	    if (subeventFound == false) {
		currentFpga = -1;
		// loop over registered endpoint addresses
		for(int i = 0; i < 4; i++) {
		    //cerr<<"FPGA: "<<i<<" "<<fpgasAddr[i]<<endl;
		    // found a matching one
		    //cout<<" *data " << ((UInt_t)(*data&0xFFFF))<<endl;
		    if(*data&0xFFFF==0x0367)//cout<<"HEY ESTO VA BIEN"<<endl;
			if (((*data) & 0xffff) == fpgasAddr[i]) {
			    endOfSubevent = data + (((*data) & 0xffff0000) >> 16) + 1;
			    data++;
			    subeventFound = true;
			    currentFpga = i;

			    if(!quietMode) printf("Subevent found on fpga %d\n", currentFpga);
			    //if(quietMode) printf("Subevent found on fpga %d\n", currentFpga);
			    //cout<<"Subevent found on fpga  "<<currentFpga<<endl;

			    break;
			}
		}

		// if none matches, skip to the next subevent
		if (currentFpga == -1) {
		    data += (((*data) & 0xffff0000) >> 16) + 1;
		}
		continue;
	    }

	    }
            */
	// in case of single trb setup
	//else {


	    // gk 20.12.10
	    // fpgaAddr is a number used to select the source of data (given fpga) - set in constructor
	    //65535

	    /*      if((((*data) & 0xffff) != fpgaAddr) && (subeventFound == false)) {
	     data += ((*data) & 0xffff0000) >> 16;
	     data++;
	     continue;
	     }
	     else if((((*data) & 0xffff) == 65535) && (subeventFound == false)) {
	     end = data + (((*data) & 0xffff0000) >> 16) + 1;
	     subeventFound = true;
	     data++;
	     continue;
	     }  */
	    /*
	     cout<<"what i have here>?   "<<((*data) & 0xffff)<<endl;
	     if (subeventFound == false) {
	     currentFpga = -1;
	     // loop over registered endpoint addresses
	     for(int i = 0; i < fpgasNum; i++) {
	     cerr<<"FPGA: "<<i<<" "<<fpgasAddr[i]<<endl;
	     cout<<"FpgaAddress for comparing : "<< (UInt_t)((*data) & 0xffff)<<endl;
	     // found a matching one
	     if (((*data) & 0xffff) == fpgasAddr[i]) {
	     endOfSubevent = data + (((*data) & 0xffff0000) >> 16) + 1;
	     data++;
	     subeventFound = true;
	     currentFpga = i;

	     if(!quietMode) printf("Subevent found on fpga %d\n", currentFpga);
	     if(quietMode) printf("Subevent found on fpga %d\n", currentFpga);
	     cout<<"Subevent found on fpga  "<<currentFpga<<endl;

	     break;
	     }
	     }

	     // if none matches, skip to the next subevent
	     if (currentFpga == -1) {
	     data += (((*data) & 0xffff0000) >> 16) + 1;
	     }
	     continue;
	     }
	     */

	    subeventFound = true;
	    currentFpga = 0;
	//}

	//if( (((Int_t)((*data) & 0xffff)-99)%100==0 || ((Int_t)((*data) & 0xffff)-71)%100==0 )&& ((Int_t)((*data) & 0xffff))>100  )
	//cout<<"FpgaAddress: "<< (UInt_t)((*data) & 0xffff)<<endl;

	// process data in case subevent is found
	if (subeventFound == true) {

	    if(dataword == 0xDEADFACE)
	    {
		break;
	    }

	    TdcId = (dataword >> 24) & 0xF;// might be wrong for TRB board
	    TdcId = nCountTDC;

	    if (TdcDataLen > 0)
		TdcDataLen++;
	    switch (dataword >> 28)
	    {// Raw TDC Data
	    case 0:
		{// Group Header
		    // gk reset tdc  counter for each TRB board
		    //cout<<" new board found "<<endl;
		    //cout<<"address "<<(UInt_t)(*data&0xffff)<<endl;
		    nCountTDC = 0;
		}
	    case 2:
		{// TDC Header


		    if(!quietMode)
			printf("TRB unpack: Found TDC %d Header $%04X $%04X\n",TdcId,(dataword>>12)&0xFFF,dataword&0xFFF);

		    if (nCountTDC > 0 && nTdcEvtId != ((dataword >> 12) & 0xFFF))
		    {
			//wk    if(!quietMode)Error("TRB unpack","TDCs have different EventIds ******* Event Mixing *******");
			if (!quietMode)
			    printf("nTdcEvtId: %06X   dataword:  %06X  nEvtNr: %02X\n",
				   nTdcEvtId, ((dataword >> 12) & 0xFFF), nEvtNr);
			//               exit();
			//               return(kFALSE);
		    }
		    if (nEvtNr != ((dataword >> 12) & 0xFF))
		    {
			//wk    if(!quietMode)Error("TRB unpack","TDC EventIds != Main EventId ******* Event Mixing *******");
			if (!quietMode)
			    printf("nTdcEvtId: %06X   dataword:  %06X  nEvtNr: %02X\n",
				   nTdcEvtId, ((dataword >> 12) & 0xFFF), nEvtNr);
			//               exit();

			//               return(kFALSE);
		    }
		    if(!quietMode)
			printf("nTdcEvtId: %06X   dataword:  %06X  nEvtNr: %02X\n" , nTdcEvtId ,((dataword>>12)&0xFFF),nEvtNr);
		    nTdcEvtId = (dataword >> 12) & 0xFFF;

		    TdcDataLen = 1;
		    break;
		}
	    case 3:
		{// TDC Trailer

		    if(!quietMode)
			printf("TRB unpack: Found TDC %d Trailer $%04X $%04X\n",TdcId,(dataword>>12)&0xFFF,dataword&0xFFF);

		    if (TdcDataLen != (dataword & 0xFFF))
		    {
			if (!quietMode)
			    printf("TRB unpack: TdcDataLen %d != %d ", TdcDataLen, dataword & 0xFFF);
		    }
		    TdcDataLen = 0;
		    nCountTDC++;
		    break;
		}

	    case 4:
		{// TDC DATA Leading and width

		    if(!quietMode)
			printf("TRB unpack: Found TDC %d Lead Data $%08X\n",TdcId,dataword);

		    // gk 13.02.12 added support for VHR mode
		    if(VHR == false) {
			Int_t nData, nChannel, nWidth;
			nChannel = (dataword >> 19) & 0x1f; // decode channel

			nChannel += (TdcId * 32);
			nData = dataword & 0x7ffff; // decode 19bit data

			if(!quietMode)
			    printf("(Chan,Data) %3d, %d\n",nChannel,nData);

			// gk in case the search window is defined
			if (nData >= minWindow && nData <= maxWindow && minWindow != -100000) {
			    fill_lead(nChannel, nData);
			    //foundLeadingEdge = true;
			}
			// operate without search window
			else if(minWindow == -100000) {
			    fill_lead(nChannel, nData);
			    //foundLeadingEdge = true;
			}
		    }
		    else {  // VHR decoding
			Int_t nData, nChannel, nTDC;
			nChannel  = (dataword >> 21) & 0x7;
			nTDC      = (dataword >> 24) & 0xf;
			nChannel += (nTDC * 32);
			nData     = ((dataword & 0x7ffff) << 2) | ((dataword >> 19) & 0x3);

			if (!quietMode)
			    printf("VHR(Chan,Data) %3d, %d\n",nChannel,nData);

			fill_lead(nChannel, nData);
			//foundLeadingEdge = true;
		    }
		    //}
		    break;
		}
	    case 5:
		{// TDC DATA Trailing

		    if(!quietMode)
			printf("TRB unpack: Found TDC %d Trail Data $%08X\n",TdcId,dataword);
		    //if(!quietMode) printf("FPGA code: %s amount: %d\n", fpgaAddr, fpgasNum);

			// gk 13.02.12 added support for VHR mode
			if (VHR == false) {
			    Int_t nData, nChannel;
			    nChannel = (dataword >> 19) & 0x1f; // decode channel

			    //shift by tdc number and endpoint number
			    nChannel += (TdcId * 32);

			    nData = dataword & 0x7ffff; // decode 19bit data

			    if(!quietMode)
				printf("(Chan,Data) %3d, %d\n",nChannel,nData);

			    // this is for SINGLE LEADING/TRAILING EDGE measurements only!!!
			    //if (foundLeadingEdge == true) {
			    fill_trail(nChannel, nData);
			    //  foundLeadingEdge = false;
			    //} else {
			    // if(!quietMode) printf("!!!Rejecting that trail edge \n");
			    //}
			}
			else {
			    Int_t nData, nChannel, nTDC;
			    nChannel  = (dataword >> 21) & 0x7;
			    nTDC      = (dataword >> 24) & 0xf;
			    nChannel += (nTDC * 32);
			    nData     = ((dataword & 0x7ffff) << 2) | ((dataword >> 19) & 0x3);

			    if (!quietMode)
				printf("VHR(Chan,Data) %3d, %d\n",nChannel,nData);

			    fill_trail(nChannel, nData);
			    //foundLeadingEdge = false;
			}
			break;
		}
	    case 6:
		{// TDC ERROR


		    printDebug = true;

		    incErrorNr();
		    if ((dataword & 0x7FFF) == 0x1000)
		    {// special case for non fatal errors
			//wk     if(!quietMode)Info("TRB unpack","TDC Event Size Limit exceeded!\n");
			if (!quietMode)
			    printf("(TDC %d Error Event Size Limit: $%08X)\n", TdcId, dataword);
		    }
		    else
		    {
			//wk     if(!quietMode)Info("TRB unpack","Found TDC Error(s)!\n");
			if (!quietMode)
			    printf("TDC %d Error $%04X ($%08X)\n", TdcId, dataword & 0x7FFF, dataword);
			if (!quietMode)
			    PrintTdcError(dataword & 0x7FFF);
		    }
		    break;
		}
	    case 7:
		{// Debug Info
		    //wk      if(!quietMode)Error("TRB unpack","Found DEBUG Info");
		    if (!quietMode)
			printf("TRB unpack: TDC %d: Found Debug Info $%08X", TdcId, dataword);
		    break;
		}
	    default:
		{// not defined!
		    //wk       if(!quietMode)Error("TRB unpack","Found UNDEFINED data");
		    if (!quietMode)
			printf("TRB unpack: TDC %d: Found undefined $%08X", TdcId, dataword);
		    break;
		}
	    }

	}  //  end of if(subeventFound == true)


	++data;

	// reset flag when the end of subevent is reached
	if (data == endOfSubevent) {
	    subeventFound = false;
	}

    }

    if(!quietMode)printf("==== Unpacker end (%d)\n",subEvtId);


    //if (LeadingMult[31] == 0 || LeadingMult[63] == 0 || LeadingMult[95] == 0 || LeadingMult[127] == 0) {
    if(LeadingMult[31] != TrailingMult[31] || LeadingMult[63] != TrailingMult[63] || LeadingMult[95] != TrailingMult[95]) {
	if(!quietMode) printf("%d, %d, %d, %d, %d, %d \n", LeadingMult[31], TrailingMult[31], LeadingMult[63], TrailingMult[63], LeadingMult[95],  TrailingMult[95]);
	printDebug = true;
    }


    if(!quietMode && printDebug == true) {
	//if(true) {
	data = pSubEvt->getData();
	UInt_t *start = data;
	while (data < end)
	{
	    UInt_t dataword;
	    dataword = *data;

	    printf("%08x\t", dataword);

	    if (((data - start) % 4 == 3) && (data != start)) printf("\n");

	    data++;
	}

	printf("\n============================ END OF SUB EVENT ======================================\n\n");
    }

    return (kTRUE);
}

//wk from HldUnpack
//______________________________________________________________________________
HldEvent::HPP const HldEvent::getpSubEvt(void)
{
    //Return a pointer to the subevent read by this unpacker
    return &pSubEvt;
}

Int_t HldEvent::getLeadingTime(Int_t channel, Int_t mult) const
//______________________________________________________________________________
{
    if ((channel < kMaxChannelNr) && (mult < kMaxMult)) {
	return LeadingTime[channel][mult];

    }
    return -1; //channel or multiplicity out of range
}

Int_t HldEvent::getWidthTime(Int_t channel, Int_t mult) const
//______________________________________________________________________________
{

    if ((channel < kMaxChannelNr) && (mult < kMaxMult))
	return WidthTime[channel][mult];
    return -1; //channel or multiplicity out of range
}

Int_t HldEvent::getTrailingTime(Int_t channel, Int_t mult) const
//______________________________________________________________________________
{

    if ((channel < kMaxChannelNr) && (mult < kMaxMult))
	return TrailingTime[channel][mult];
    return -1; //channel or multiplicity out of range
}

Int_t HldEvent::getLeadingMult(Int_t channel) const
//______________________________________________________________________________
{

    if (channel < kMaxChannelNr)
	return LeadingMult[channel];
    return -1; //channel out of range
}

Int_t HldEvent::getTrailingMult(Int_t channel) const
//______________________________________________________________________________
{

    if (channel < kMaxChannelNr)
	return TrailingMult[channel];
    return -1; //channel out of range
}

Int_t HldEvent::getSpikesCtr(Int_t channel) const
//______________________________________________________________________________
{

    if (channel < kMaxChannelNr)
	return SpikesCtr[channel];
    return -1; //channel out of range
}
//end
Int_t HldEvent::decode(void) {

  clearAll();

  UInt_t nEvt         = 0;  // Evt SeqNumber
  UInt_t nSizeCounter = 0;  // should go to the class member !!!!!
  UInt_t nTdcEvtId    = 0;
  UInt_t TdcDataLen   = 0;
  Int_t TdcId         = 0;
  uSubBlockSize    = 0;
  uTrbNetAdress    = 0;
  trbDataVer       = 0;
  trbExtensionSize = 0;
  trbDataPairFlag  = kFALSE;
  Int_t nCountTDC = 0;
  Bool_t  printDebugSpecial = kFALSE;

  Bool_t tryRecover_1 = kFALSE;
  Bool_t tryRecover_2 = kFALSE;
  Int_t DEBUG_LEVEL = 0;
  //Bool_t doINLCorrection = kFALSE;
  Bool_t kDsSkip = kTRUE;
  //Bool_t printTdcError  = kFALSE;


  Bool_t wasRecovered=kFALSE;

  //uStartPosition = 0;

  UInt_t* data = pSubEvt->getData();
  UInt_t* end  = pSubEvt->getEnd();
  /*
  UInt_t* data = pSubEvt->getData();
  UInt_t* end = pSubEvt->getEnd();
  */
  //nEvt = gHades->getCurrentEvent()->getHeader()->getEventSeqNumber();


  if(debugFlag > 0) {
    cout<<"-0-EvNb "<<nEvt<<" -nSizeCounter: "<<nSizeCounter<<", "<<hex<<"data word: "<<*data<<dec<<"  --NetSubevent(TRB) Header"<<endl;
  }

  // data+=1;
  // jump to current SubSubEvent
  //data = data + uStartPosition;
  UInt_t* data_c1 = data; // copy for TRB extension loop
  UInt_t* data_c2 = data; // copy for TDC loop


  if(debugFlag > 0){
    cout<<"-0a-EvNb "<<nEvt<<" -uStartPosition: "<<uStartPosition<<", "<<hex<<"data word: "<<*data<<dec<<"  --"<<endl;
  }

  data+=2;
  uSubBlockSize=(*data>>16)&0xFFFF;   // No. of the data words from one TRB board
  uTrbNetAdress=*data&0xFFFF;         // TRB net Adddress (Old TRB subevent id)
  if(uSubBlockSize > (UInt_t)(end-data)){
      Error("decode()","Event %d --> SubBlkSize out of subevt. Blksize=%d.",nEvt,uSubBlockSize);
      cout<<"SIZE "<<(UInt_t)(end-data)<<endl;
  }
    //data+=2;

  /*
  if(uSubBlockSize > (UInt_t)(end-data)){
      if(tryRecover_1)Error("decode()","Event %d --> SubBlkSize out of subevt. Blksize=%d. Try to recover.",nEvt,uSubBlockSize);
      else            Error("decode()","Event %d --> SubBlkSize out of subevt. Blksize=%d.",nEvt,uSubBlockSize);

      if(!tryRecover_1) {
	  uStartPosition=0;
	  return 0;
      }

      //-----------------------------------------------------------------------------
      // try to recover the data by looping forward until next
      // valid address is found. (apr12 beam, day 116)
      Int_t        ct= 0;
      Bool_t reCover = kFALSE;

      while(data<end){ // still inside subevent
	  ct++;
          data +=1;
	  uSubBlockSize=(*data>>16)&0xFFFF;   // No. of the data words from one TRB board
          uTrbNetAdress=*data&0xFFFF;         // TRB net Adddress (Old TRB subevent id)
	  uStartPosition+=1;

	  if( ((   uTrbNetAdress>=0x4000
	        && uTrbNetAdress<=0x5000 )
                || uTrbNetAdress==0x5555)
	     && uSubBlockSize>0
	     && uSubBlockSize<200){

	      data_c1 = data;
              data_c2 = data;
	      reCover = kTRUE;
              wasRecovered=kTRUE;
	      break;
	  }
      }
      if(!reCover) { // bad luck
	  cout<<"Recover failed."<<endl;
	  uStartPosition = 0;
	  return 0;
      }
      //-----------------------------------------------------------------------------

    }
    */

  // check TRBNetAdress
  // IF TRBNet debug information
  /*
  if ( uTrbNetAdress == 0x5555 ) {
    //decodeTrbNet(data,subEvtId);
    data+=(uSubBlockSize)+1;

    uStartPosition = 0; // reset the Start position for the next Event
    if( debugFlag > 0){
      cout << "++++> stopping decode(), dataword: " << hex << *data << dec << endl;
    }
    return 100;
  }
  */
  // IF CTS -> skip & 
  // skip if its a hub in the 0x8000 range (added by JanM 20140430
  /*
  if(uTrbNetAdress <= 0x00FF || ((uTrbNetAdress & 0xF000) == 0x8000)){
    return 10;
  }
  */
  /*
  // ELSE unpack TRB - TDC data
  // ------------------------- 
  // 1. Should INL correction be done ?
  // -------------------------
  if (correctINL == kTRUE){

    trbinlcorr = trbaddressmap->getBoard(uTrbNetAdress);
    
    // if there is no mapping table it is not possible to decide if the INL correction should be done
    if(trbinlcorr == NULL ){
      Error("decode()","Event %d -->  No mapping table data for TRB net adress 0x%x",nEvt,uTrbNetAdress);
      return kFALSE;
    }
    
    else{
      
      // decide if INL correction should be done for this board
      if (trbinlcorr->getNChannels() > 0  ){
	correctINLboard = kTRUE;
      }
      else{
	if( debugFlag > 0) Warning("decode()","Event %d -->  No INL corr data for TRB net adress 0x%x",nEvt,uTrbNetAdress);
	correctINLboard = kFALSE;
      }
    }
  }
  */
  // ------------------------- 
  // 2. control the SubSubEvent size
  // -------------------------
  /*
  UInt_t uSubBlockSizeTest = (*(data+1)&0xFFFF);
  if(uSubBlockSize != uSubBlockSizeTest){
    if(tryRecover_2)Error("decode","Event %d -->  SubBlkSize from TRBnet ( %04x ) != SubBlkSize from subsubevt ( %04x ). Try to recover.",nEvt,uSubBlockSize,*(data+1)&0xFFFF);
    else            Error("decode","Event %d -->  SubBlkSize from TRBnet ( %04x ) != SubBlkSize from subsubevt ( %04x )",nEvt,uSubBlockSize,*(data+1)&0xFFFF);

    if(!tryRecover_2) {
	uStartPosition=0;
	return 0;
    }


    //-----------------------------------------------------------------------------
    // try to recover the data by looping forward until next
    // valid address is found and correct the uSubBlockSize  (apr12 beam, day 116)
    UInt_t        ct= 0;
    Bool_t reCover = kFALSE;

    UInt_t* data2=data;
    UInt_t  uSubBlockSize2;
    UInt_t  uTrbNetAdress2;
    UInt_t  uStartPosition2=uStartPosition;
    while(data2<end){ // still inside subevent
	ct++;
	data2 +=1;
	uSubBlockSize2=(*data2>>16)&0xFFFF;   // No. of the data words from one TRB board
	uTrbNetAdress2=*data2&0xFFFF;         // TRB net Adddress (Old TRB subevent id)
	uStartPosition2+=1;

	if( ((   uTrbNetAdress2>=0x4000
	      && uTrbNetAdress2<=0x5000 )
	     || uTrbNetAdress2==0x5555)
	   && uSubBlockSize2>0
	   && uSubBlockSize2<200){

	    if(ct == uSubBlockSizeTest+1) {

		uSubBlockSize =uSubBlockSizeTest;
		reCover = kTRUE;
		wasRecovered=kTRUE;
                break;
	    }
	    if(ct == uSubBlockSize+1) {

		reCover = kTRUE;
		wasRecovered=kTRUE;
                break;
	    }
	}
    }
    if(!reCover) {
	cout<<"Recover failed."<<endl;
	uStartPosition = 0;
	return 0;
    }
    //-----------------------------------------------------------------------------
  }
  
  if( debugFlag > 0){
    cout<<"-1a-EvNb "<<nEvt<<" -nSizeCounter: "<<nSizeCounter<<", "<<hex<<"data word: "<<*data<<dec
	<<" SubSize: "<<uSubBlockSize<< hex << " TrbNetAdress "<<uTrbNetAdress<<dec<<endl;
  }
  */
  // ------------------------- 
  // 3. Loop over the TRB extension
  // ------------------------- 
  data_c1+=2;// go to the dataword where TRB extension is defined
  trbDataVer = (*data_c1>>24)&0xFF;
  /*
  trbDataPairFlag = (*data_c1>>16)&0x1;

  trbExtensionSize = *data_c1&0x0000FFFF;

  if(!trbExtensionSize)cout<<" WTF "<<endl;
  if((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode)) {
    Info("decode()","Event %d -->  Unpacker TrbHeader: %08x  Ver:%d  Ext:%d Pair:%d \n",nEvt,*data_c1,trbDataVer,trbExtensionSize,(Int_t)trbDataPairFlag);
  }
  if( trbExtensionSize != 0 ) {
    for ( UInt_t ii = 0; ii < trbExtensionSize; ii++ ) {

      data_c1++;    
      if(trbExtensionSize < 128) {
	trbDataExtension[ii] = *data_c1; 
      } else{
	Error("decode()","Event %d -->  SubEventId = 0x%x (TRB ID 0x%x ) too many data words (%d), limit is 128",nEvt,subEvtId,uTrbNetAdress,trbExtensionSize);
	return 0;
      }
    }
  }
  */
  // ------------------------- 
  // 4. set resolution mode 
  // data_c1 =  AATTXXX
  // AA = 0   high resolution (100 ps)
  // AA = 1   very  high resolution (25 ps)
  // AA = 2   very  high resolution (25 ps) - calibration
  // ------------------------- 
  if( trbDataVer == 0 ){
    highResModeOn = kFALSE;
  }
  else if( trbDataVer == 1 ){
    highResModeOn = kTRUE;
  }    
  else{
    Warning("decode()","Event %d -->  Unknown trbDataVer ( 0x%x ) for TRB net adress 0x%x",nEvt,trbDataVer,uTrbNetAdress);
  }
  
  if( debugFlag > 0 ) cout<<" High Resolution mode: "<<highResModeOn<< hex << " for TRB net adress: "<<uTrbNetAdress<<dec<<endl;


  //nSizeCounter++;
  // ------------------------- 
  // 5. Loop over TDCs
  // ------------------------- 
  //data_c2+=(2+trbExtensionSize);// jump to last TrbHeader dataword
  data_c2+=(2);// jump to last TrbHeader dataword
  //for ( UInt_t ii = 0; ii < uSubBlockSize-trbExtensionSize-2; ii++ ) {
  //for ( UInt_t ii = 0; ii < uSubBlockSize-2; ii++ ) {
  while (data_c2 < end) {

      nSizeCounter++;

      UInt_t dataword = *data_c2;
      if(debugFlag > 0) cout<<"--EvNb "<<nEvt<<" data_"<<nSizeCounter<<": "<<hex<<dataword<<dec<<endl;

      if(dataword==0xDEADFACE) {
	  if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
	      Info("decode()","Event %d -->  Found DEADFACE -> break %08X %08X\n",nEvt,*data,*end);
	  }
	  break;
      }

      TdcId = (dataword>>24)&0xF;      // might be wrong for TRB board
      if(TdcId>3) TdcId = nCountTDC;
      if(TdcId>3)  {
	  cout<<"TDcId larger that 3 !! skipping TDC!!   "<<endl;
	  printf("%08x\t", dataword);
	  printDebugSpecial = kTRUE;
	  break;
      }
      if(TdcDataLen>0) TdcDataLen++;

      // -------------------------
      // 6. check data word TYPE
      // -------------------------
      switch(dataword>>28) {
      case 0: { // Group Header
	  if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
	      Info("decode()","Event %d --> Found GLOBAL Header $%08X\n",nEvt,dataword);
	  }
	  if(!quietMode) Error("decode()","Event %d --> Global Header not expected!",nEvt);
	  if(wasRecovered) return kDsSkip; // in case of a repaired event we skip if data are problematic
	  nCountTDC=0;
	  break;
      }

      case 1: { // Group Trailer
	  if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
	      Info("decode()","Event %d -->  Found GLOBAL Trailer $%08X\n",nEvt,dataword);
	  }
	  if(!quietMode) Error("decode()","Event %d --> Global Trailer not expected!",nEvt);
	  if(wasRecovered) return kDsSkip; // in case of a repaired event we skip if data are problematic
	  break;
      }

      case 2: { // TDC Header
	  if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
	      Info("decode()","Event %d -->  Found TDC %d Header $%04X $%04X\n",nEvt,TdcId,(dataword>>12)&0xFFF,dataword&0xFFF);
	  }
	  if( TdcId>0 && nTdcEvtId!=((dataword>>12)&0xFFF)) {
	      if(!quietMode) {
		  Error("TRB unpack","Event %d -->  TDCs have different EventIds ******* Event Mixing *******",nEvt);
		  printf("nTdcEvtId: %06X   dataword:  %06X  nEvt: %02X\n" , nTdcEvtId ,((dataword>>12)&0xFFF),nEvt);
	      }
	      if(wasRecovered) return kDsSkip; // in case of a repaired event we skip if data are problematic
	  }
	  if( nEvt!=((dataword>>12)&0xFF)) {
	      if(!quietMode) {
		  //      TDC Trigger TAG problem, To be checked !!!!
		  //      Error("TRB unpack","Event %d -->  TDC EventIds != Main EventId ******* Event Mixing *******");
		  //      printf("SubEvtId: %06X nTdcEvtId: %06X   dataword:  %06X  nEvt: %02X\n" ,subEvtId, nTdcEvtId ,((dataword>>12)&0xFFF),nEvt);
	      }
	      //          exit();
	      //          return(kFALSE);
	  }
	  if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
	      Info("decode()","nTdcEvtId: %06X   dataword:  %06X  nEvt: %02X\n" , nTdcEvtId ,((dataword>>12)&0xFFF),nEvt);
	  }
	  nTdcEvtId=(dataword>>12)&0xFFF;

	  TdcDataLen=1;
	  break;
      }

      case 3: { // TDC Trailer
	  if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
	      Info("decode()","Found TDC %d Trailer $%04X $%04X\n",TdcId,(dataword>>12)&0xFFF,dataword&0xFFF);
	  }
	  if(TdcDataLen!=(dataword&0xFFF)) {
	      if(!quietMode) Error("decode()","Event %d -->  TRB unpack : TdcDataLen!= length in Trailer!",nEvt);
	      if(!quietMode) printf("TRB unpack: TdcDataLen %d != %d ",TdcDataLen,dataword&0xFFF);
	      if(wasRecovered) return kDsSkip; // in case of a repaired event we skip if data are problematic
	      return kDsSkip;
	  }
	  TdcDataLen=0;
	  if(nTdcEvtId!=((dataword>>12)&0xFFF)) {
	      if(!quietMode) Error("TRB unpack","Event %d -->  TDC Header and Trailer have different EventIds",nEvt);
	      //          exit();
	      //          return(kFALSE);
	      return kDsSkip;
	  }
	  nCountTDC++;
	  break;
      }

      case 4: { // TDC DATA Leading
	  if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
	      Info("decode()","Event %d -->  Found TDC %d Lead Data $%08X\n",nEvt,TdcId,dataword);
	  }
	  Int_t nData, nChannel;
	  Int_t nDataLow, nDataHigh;

	  if(!highResModeOn) { // 100ps binning
	      nChannel=(dataword>>19)&0x1f; // decode channel
	      nChannel+=TdcId*32;
	      nData=dataword&0x7ffff;   // decode 19bit data
	      if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
		  Info("decode()","Event %d --> (Chan,Data,InlCorrData) %3d, %d, %f \n",nEvt,nChannel,nData, doINLCorrection(nChannel, nData));
	      }
	      if(nChannel > 127){
		  Error("decode","Event %d -->  channel %i number bigger than 128",nEvt,nChannel);
		  if(wasRecovered) return kDsSkip; // in case of a repaired event we skip if data are problematic
		  return kDsSkip;
	      }
	  }
	  else {  // 25ps binning mode
	      nChannel=(dataword>>21)&0x7; // decode channel
	      nChannel+=TdcId*8;
	      nDataHigh = (dataword&0x7ffff)<<2;
	      nDataLow  = (dataword>>19)&0x3;
	      nData     = nDataHigh + nDataLow;

	      if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
		  Info("decode()","Event %d -->  (Chan,Data) %3d, %d\n",nEvt,nChannel,nData);
	      }
	      if(nChannel > 31){
		  Error("decode","Event %d -->  channel number bigger than 31",nEvt);
		  if(wasRecovered) return kDsSkip; // in case of a repaired event we skip if data are problematic
		  return 0;
	      }
	  }

	  // this is for SINGLE LEADING/TRAILING EDGE measurements only!!!
	  // No check if the order is correct!!!
	  // i am depending on the TDC to deliver the right order

	  if(trbDataPairFlag) {
	      if(!fill_pair(nChannel,(Float_t)(nData&0xFFFF),(Float_t)((nData>>12)&0x7F))) {

	      }
	  }
	  else { // Leading as usual
	      // fill INL corrected time
	      if(!fill_lead(nChannel, doINLCorrection(nChannel, nData))) {

	      }
	  }

	  break;
      }
      case 5: { // TDC DATA Trailing
	  if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
	      Info("decode()","Event %d --> Found TDC %d Trail Data $%08X\n",nEvt,TdcId,dataword);
	  }
	  Int_t nData, nChannel;
	  Int_t nDataLow, nDataHigh;
	  if(!highResModeOn){ //100ps binning
	      nChannel=(dataword>>19)&0x1f; // decode channel
	      nChannel+=TdcId*32;
	      nData=dataword&0x7ffff;   // decode 19bit data
	      if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
		  Info("decode()","Event %d --> (Chan,Data,InlCorrData) %3d, %d, %f \n",nEvt,nChannel,nData, doINLCorrection(nChannel, nData));
	      }
	      if(nChannel > 127){
		  Error("decode","Event %d -->  channel %i number bigger than 128",nEvt,nChannel);
		  if(wasRecovered) return kDsSkip; // in case of a repaired event we skip if data are problematic
		  return kDsSkip;
	      }

	  } else {  // 25ps binning mode
	      nChannel=(dataword>>21)&0x7; // decode channel
	      nChannel+=TdcId*8;
	      nDataHigh = (dataword&0x7ffff)<<2;
	      nDataLow  = (dataword>>19)&0x3;
	      nData     = nDataHigh + nDataLow;
	      if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
		  Info("decode()","Event %d --> (Chan,Data) %3d, %d\n",nEvt,nChannel,nData);
	      }
	      if(nChannel > 127){
		  Error("decode","Event %d -->  channel %i number bigger than 128",nEvt,nChannel);
		  if(wasRecovered) return kDsSkip; // in case of a repaired event we skip if data are problematic
		  return kDsSkip;
	      }


	  }

	  // this is for SINGLE LEADING/TRAILING EDGE measurements only!!!
	  // No check if the order is correct!!!
	  // i am depending on the TDC to deliver the right order

	  // fill INL corrected time
	  if(!fill_trail(nChannel, doINLCorrection(nChannel, nData) )) {

	  }
	  break;
      }
      case 6:{// TDC ERROR
	  if((dataword&0x7FFF)==0x1000) {// special case for non fatal errors
	      if(!quietMode)Error("TRB unpack","Event %d -->  TDC Event Size Limit exceeded!\n",nEvt);
	      if(!quietMode)printf("(TDC %d Error Event Size Limit: $%08X)\n",TdcId,dataword);
	  } else {
	      if(!quietMode) {
		  Error("TRB unpack","Event %d -->  Found TDC Error(s)!\n",nEvt);
		  printf("TDC %d Error $%04X ($%08X)\n",TdcId,dataword&0x7FFF,dataword);
		  printTdcError(dataword&0x7FFF,subEvtId);

	      } /*else if(reportCritical && (dataword&0x6000)!=0) {
	      Error("TRB unpack","Event %d -->  Found CRITICAL error",nEvt);
	      if((dataword&0x2000)!=0) {
	      Error("TRB unpack","Event %d -->  Event lost (trigger FIFO overflow)",nEvt);
	      printf("TDC %d Error $%04X ($%08X)\n",TdcId,dataword&0x7FFF,dataword);
	      } else if((dataword&0x4000)!=0) {
	      Error("TRB unpack","Event %d -->  Internal fatal chip error has been detected",nEvt);
	      printf("TDC %d Error $%04X ($%08X)\n",TdcId,dataword&0x7FFF,dataword);
	      }
	      }   */
	  }
	  break;
      }
      case 7: { // Debug Info
	  if(!quietMode)Error("TRB unpack","Event %d -->  Found DEBUG Info",nEvt);
	  if(!quietMode)printf("TRB unpack: TDC %d: Found Debug Info $%08X",TdcId,dataword);
	  break;
      }
      default: {// not defined!
	  if(!quietMode) Error("TRB unpack","Event %d -->  Found UNDEFINED data",nEvt);
	  if(!quietMode) printf("TRB unpack: TDC %d: Found undefined $%08X",TdcId,dataword);
	  break;
      }
      }
      data_c2++;
  }// end of for-loop over TRB datawords

  /*
  if(uSubBlockSize-trbExtensionSize-2!=nSizeCounter) {
    cout<<" --Block size ="<<uSubBlockSize<<" , counter Size: "<<nSizeCounter <<endl; 
    if(!quietMode)Error("TRB unpack","Event %d -->  Blocksize!=Counted words!!!",nEvt);
    }
    */

  if(printDebugSpecial== kTRUE) {
	//if(true) {
	data = pSubEvt->getData();
	UInt_t *start = data;
	while (data < end)
	{
	    UInt_t dataword;
	    dataword = *data;

	    printf("%08x\t", dataword);

	    if (((data - start) % 4 == 3) && (data != start)) printf("\n");

	    data++;
	}

	printf("\n============================ END OF SUB EVENT ======================================\n\n");
    }







  if( ((debugFlag > 0 || DEBUG_LEVEL>4) && (!quietMode))) {
    Info("decode()","Event %d --> ===Unpacker end (SubEvtID 0x%x , TRB ID 0x%x ) ===================\n",nEvt,subEvtId,uTrbNetAdress);
  }
  return(kTRUE);
}

// correct INL of HPTDC for a single DATA word  
Float_t HldEvent::doINLCorrection(Int_t nTrbCh, Int_t nRawTime)
{  
  Int_t nTimeHigh, nTimeLow;
  Float_t fCorrectedTime;
  /*
  if(correctINL == kTRUE  && correctINLboard == kTRUE ){
    if( highResModeOn == kTRUE ){            //correct INL for 25ps/channel
      
      nTimeHigh = ( nRawTime & 0x1FFC00);  // 11 MSB
      nTimeLow  = ( nRawTime & 0x003ff );  // 10 LSB 
      
      fCorrectedTime = (Float_t)( (Float_t)nTimeLow  +
				  ( trbinlcorr->getCorrection( (nTrbCh*4)+(Int_t)(nTimeLow / 256.0), (nTimeLow % 256) ) ));
      fCorrectedTime = (Float_t)( (Float_t)nTimeHigh + fCorrectedTime);
      
      if(debugFlag > 0){
	cout<<"==TRB_Res_25ps, chan: "<<nTrbCh<<", measured time: "<<nRawTime<<", corrected:"<<fCorrectedTime<<
	  ", correction:"<<( trbinlcorr->getCorrection( (nTrbCh*4)+(Int_t)(nTimeLow / 256.0), (nTimeLow % 256) ) )
	    <<"  --BIN: "<<nTimeLow<<" newCh: "<<(nTrbCh*4)+(Int_t)(nTimeLow / 256.0)<<" binInParHist: "<<( (nTimeLow % 256)  )<< endl;
      }
    } 
    
    else{      //correct INL for 100ps/channel
      
      nTimeHigh = ( nRawTime & 0x7ff00 );
      nTimeLow  = ( nRawTime & 0x000ff );  
      fCorrectedTime = (Float_t)( (Float_t)nTimeLow  +
				  ( trbinlcorr->getCorrection( nTrbCh , nTimeLow ) ));
      fCorrectedTime = (Float_t)( (Float_t)nTimeHigh + fCorrectedTime);
      
      if(debugFlag > 0){
	cout<<"==chan: "<<nTrbCh<<", measured time: "<<nRawTime<<", corrected:"<<fCorrectedTime<<
	  ", correction:"<<( trbinlcorr->getCorrection( nTrbCh , nTimeLow ) )<<endl;
      }
    }
  }

  */
  //no INL correction will be done
  //else{
  fCorrectedTime = (Float_t)nRawTime;
  //}
  
  return fCorrectedTime;
}
Bool_t HldEvent::fill_pair(Int_t ch,Float_t time,Float_t length) {

  ///////////////////////////////////////////
  // Stores the given time in the next data element
  // and sets the multiplicity.
  // Return kFALSE if 10 hits are already stored.
  ///////////////////////////////////////////

  if( trbLeadingMult[ch]<10) {
    trbLeadingTime[ch][trbLeadingMult[ch]]=time;
    trbTrailingTime[ch][trbLeadingMult[ch]]=time+length;
  }
  trbLeadingMult[ch]++;
  trbTrailingMult[ch]=trbLeadingMult[ch];
  return (trbLeadingMult[ch]<=10);

}

Bool_t HldEvent::fill_lead(Int_t ch, Float_t time) {

  ///////////////////////////////////////////
  // Stores the given time in the next data element
  // and sets the multiplicity.
  // Return kFALSE if 10 hits are already stored.
  ///////////////////////////////////////////

  if( trbLeadingMult[ch]<10){
    trbLeadingTime[ch][trbLeadingMult[ch]]=time;
  }
  trbLeadingMult[ch]++;
  return(trbLeadingMult[ch]<=10);
}


Bool_t HldEvent::fill_trail(Int_t ch, Float_t time) {

  ///////////////////////////////////////////
  // Calculates the time between trailing and LAST(!) leading hit.
  // No other check if its really the right one,
  // i am depending on the TDC to deliver the right order
  // Return kFALSE if no leading yet or more than 4 Hits
  ///////////////////////////////////////////

  trbTrailingTotalMult[ch]++;
   
  Int_t m = -1;     // Leading Multiplicity
  m=trbLeadingMult[ch];
  if(m==0) return kFALSE;

  if( m<=10) {
    if( trbTrailingMult[ch]!=m) {
      trbTrailingTime[ch][m-1]=time;
      trbADC[ch][m-1] = time - trbLeadingTime[ch][m-1];
      if(trbADC[ch][m-1]<0) { // either a big error or an overflow because off too big time windows
	trbADC[ch][m-1]+=0x80000; // correct for overflow
      }
    } else {
      return kFALSE; // In this case we already have one trailing
    }
  }
  trbTrailingMult[ch]=m;

  return(trbTrailingMult[ch]<=10);
}
void HldEvent::printTdcError(UInt_t e, UInt_t subEvntId) {
  const Char_t *e_str[15]={
      "Hit lost in group 0 from read-out FIFO overflow",
      "Hit lost in group 0 from L1 buffer overflow",
      "Hit error have been detected in group 0",
      "Hit lost in group 1 from read-out FIFO overflow",
      "Hit lost in group 1 from L1 buffer overflow",
      "Hit error have been detected in group 1",
      "Hit lost in group 2 from read-out FIFO overflow",
      "Hit lost in group 2 from L1 buffer overflow",
      "Hit error have been detected in group 2",
      "Hit lost in group 3 from read-out FIFO overflow",
      "Hit lost in group 3 from L1 buffer overflow",
      "Hit error have been detected in group 3",
      "Hits rejected because of programmed event size limit",
      "Event lost (trigger FIFO overflow)",
      "Internal fatal chip error has been detected"
  };

  if(e==0) return;// No Error

  cout << "=== TRB/TDC Error analysis: TRB id = " <<subEvntId<< endl;
  for(Int_t i=0; i<15; i++){
    if( e&0x1){
      cout << e_str[i] << endl;
    }
    e>>=1;
  }
  cout << "===" << endl;
}

