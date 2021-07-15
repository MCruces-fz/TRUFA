#include "trpclookuptable.h"
#include "TString.h"
#include <fstream>
#include <cmath>
#include <stdlib.h>
//#include <sstream>

using namespace std;

ClassImp(RpcLookupTable);

RpcLookupTable::RpcLookupTable() {

    for(Int_t i=0;i<4;i++) {
        for(Int_t j=0;j<4;j++) {
            for(Int_t k=0;k<32;k++) {
                for(Int_t l=0;l<6;l++) {
                    pars[i][j][k][l] = -1;
                }
            }
        }
    }
    //Int_t ind[4] = {871,899,877,878};
    Int_t ind[4] = {837,871,888,899}; // T3-T4-T1 if you change this -> change also the luptab order
    //Int_t ind[4] = {888,837,871,899}; // T1-T3-T4 sorted top-bottom
    //Int_t ind[4] = {899,871,888,811};

    for(Int_t i=0;i<4;i++) {
        detInd[i] = ind[i];
    }

}
RpcLookupTable::RpcLookupTable(TString filename) {

    for(Int_t i=0;i<4;i++) {
        for(Int_t j=0;j<4;j++) {
            for(Int_t k=0;k<32;k++) {
                for(Int_t l=0;l<6;l++) {
                    pars[i][j][k][l] = -1;
                }
            }
        }
    }

    Int_t ind[4] = {837,871,888,899}; // T3-T4-T1 if you change this -> change also the luptab order
    //Int_t ind[4] = {888,837,871,899}; // T1-T3-T4 sorted top-bottom
    //Int_t ind[4] = {899,871,888,811};

    for(Int_t i=0;i<4;i++) {
        detInd[i] = ind[i];
    }

    readParams(filename);
}
void RpcLookupTable::getParams(Int_t num,Int_t mbch, Int_t ch, Int_t& cell, Int_t& col, Int_t& row, Float_t& x, Float_t& y, Float_t& z) {
    cell  = pars[num][mbch][ch][0];
    col   = pars[num][mbch][ch][1];
    row   = pars[num][mbch][ch][2];
    x     = pars[num][mbch][ch][3];
    y     = pars[num][mbch][ch][4];
    z     = pars[num][mbch][ch][5];
}
void RpcLookupTable::setParams(Int_t num,Int_t mbch, Int_t ch, Int_t cell, Int_t col, Int_t row, Float_t x, Float_t y, Float_t z) {

    pars[num][mbch][ch][0] = cell;
    pars[num][mbch][ch][1] = col;
    pars[num][mbch][ch][2] = row;
    pars[num][mbch][ch][3] = x;
    pars[num][mbch][ch][4] = y;
    pars[num][mbch][ch][5] = z;

}
Bool_t RpcLookupTable::readParams(TString filename) {
    std::ifstream file(filename.Data());
    TString name;
    Int_t detid;
    Char_t c;
    Int_t trbc, mbch, ch, cell, col, row;
    Float_t xc, yc, zc;
    Bool_t reading = kTRUE;
    Int_t num;
    while(file.get(c)) {
        if(c=='#') {
            file>>name;
            file>>name;
            //printf("name %s",name.Data());
            if(name.Contains("LUPTABLE")) {
                //
                //cout<<" LUPTABLE found "<<endl;
                file>>detid;
                file>>zc;
                //printf("id %i and z %f\n\n\n\n\n\n\n",detid,zc);
                //cout<<" detector id "<<detid<<endl;
                num = getDetNum(detid);

                if(num!=-1) {
                    //file>>zc;
                    // Now the loop over data until all
                    // the parameters for a ginven detectors are inside
                    //printf("allOk %i det id %i\n",detid,num);
                    while (reading && !file.eof()) {
                        //for(Int_t p=0;p<10;p++) {
                        //printf("allOk2 %i det id %i\n",detid,num);
                        file>>name;
                        //printf("name %s",name.Data());
                        //printf("c is %c\n",c);
                        //printf("allOk3 %i det id %i\n",detid,num);
                        if(name.Contains("#")) {
                            //printf("allOk4 %i det id %i\n",detid,num);
                            // we reached end of the loop
                            //reading = kFALSE;
                            break;
                        }
                        //printf("allOk5 %i det id %i\n",detid,num);
                        //file.unget();
                        //printf("allOk6 %i det id %i\n",detid,num);

                        trbc = atoi(name.Data());
                        //trbc = sprintf(atoi(name.Data());
                        //file>>trbc>>mbch>>ch>>cell>>col>>row>>xc>>yc;
                        file>>mbch>>ch>>cell>>col>>row>>xc>>yc;
                        //printf("num %i mbch %i ch %i cell %i col %i row %i xc %f yc %f zc %f\n",num,mbch,ch,cell,col,row,xc,yc,zc);

                        //if(num==1 &&mbch==3 && ch==31) break;
                        setParams(num,mbch-1,ch-1,cell,col,row,xc,yc,zc);
                    }
                }
            }
        }
    }
    return kTRUE;
    file.close();
}
Bool_t RpcLookupTable::writeParams(TString filename) {
    return kTRUE;
}
Int_t RpcLookupTable::getDetNum(Int_t trbId) {
    for(Int_t i=0;i<4;i++) {
        if(trbId == detInd[i]) return i;
    }
    //cout<<" TRB ID DOES NOT EXIST IN THE CONTAINER!!! ADD IT"<<endl;
    return -1;
}

