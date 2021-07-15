#include "tactivecells.h"
#include <fstream>
#include <iostream>

using namespace std;

ClassImp(TActiveCells);
TActiveCells::TActiveCells( ){

}
TActiveCells::TActiveCells(TString filename) {

    //We get from files the parameters
    fInPar = TFile::Open(filename);

    //cout<<"finpar "<<fInPar<<endl;

    for(Int_t j=0;j<4;j++) {
	h_active_cells[0][j] = (TH2D*)fInPar->Get(Form("h_h%i_2015",j));
	h_active_cells[1][j] = (TH2D*)fInPar->Get(Form("h_h%i_2016",j));
	h_active_cells[2][j] = (TH2D*)fInPar->Get(Form("h_h%i_2017",j));
        h_active_cells[3][j] = (TH2D*)fInPar->Get(Form("h_h%i_2018",j));
        h_active_cells[4][j] = (TH2D*)fInPar->Get(Form("h_h%i_2019",j));
        h_active_cells[5][j] = (TH2D*)fInPar->Get(Form("h_h%i_2020",j));
	h_active_cells[6][j] = (TH2D*)fInPar->Get(Form("h_h%i_2021",j));
        //cout<<"j "<<j<<" "<<h_active_cells[0][j]<<" "<<h_active_cells[1][j]<<endl;
    }
}

Bool_t TActiveCells::isCellActive(Int_t year, Float_t doy, Int_t detector, Int_t cell) {

    Bool_t isActive = kFALSE;
    Int_t y = 0;
    if(year==2015) y = 0;
    if(year==2016) y = 1;
    if(year==2017) y = 2;
    if(year==2018) y = 3;
    if(year==2019) y = 4;
    if(year==2020) y = 5;
    if(year==2021) y = 6;
    if(year>=2022) return kTRUE;

    if(!h_active_cells[y][detector])  {
        //cout<<"Pars not loaded "<<h_active_cells[y][detector]<<endl;
        return kTRUE;
    }

    Int_t bin = h_active_cells[y][detector]->GetXaxis()->FindBin(doy);
    isActive = h_active_cells[y][detector]->GetBinContent(bin,cell+1);
    //cout<<"bin "<<bin<<" cell "<<cell<<" "<<isActive<<endl;
    return isActive;

}
