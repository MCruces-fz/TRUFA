{

    TFile*  fIn  = TFile::Open("clean_cells_2015_2021.root");
    TFile*  fOut = new TFile("clean_cells_2015_2021_c.root", "recreate");

    TH2D*  h_active_cells[7][4];

    for(Int_t j=0;j<3;j++) {
	h_active_cells[0][j] = (TH2D*)fIn->Get(Form("h_h%i_2015",j));
	h_active_cells[1][j] = (TH2D*)fIn->Get(Form("h_h%i_2016",j));
	h_active_cells[2][j] = (TH2D*)fIn->Get(Form("h_h%i_2017",j));
        h_active_cells[3][j] = (TH2D*)fIn->Get(Form("h_h%i_2018",j));
	h_active_cells[4][j] = (TH2D*)fIn->Get(Form("h_h%i_2019",j));
	h_active_cells[5][j] = (TH2D*)fIn->Get(Form("h_h%i_2020",j));
	h_active_cells[6][j] = (TH2D*)fIn->Get(Form("h_h%i_2021",j));
	//h_active_cells[7][j] = (TH2D*)h_active_cells[6][j]->Clone(Form("h_h%i_2022", j));
        //cout<<"j "<<j<<" "<<h_active_cells[0][j]<<" "<<h_active_cells[1][j]<<endl;
    }
    /*
    for(Int_t i=2;i<4;i++){
        for(Int_t j=0;j<3;j++){
            for(Int_t bi=1;bi<=h_active_cells[i][j]->GetNbinsX();bi++){
                for(Int_t bj=0;bj<=h_active_cells[i][j]->GetNbinsY();bj++){
                   h_active_cells[i][j]->SetBinContent(bi, bj, 1);


                }
            }
        }
    }
    */
    //Int_t UnactiveCells_h0[10]={1,3,10,11,20,36,41,46,55,112};
    //Int_t UnactiveCells_h1[14]={1,3,6,8,9,10,20,30,60,100,112,113,119,120};
    //Int_t UnactiveCells_h2[14]={1,4,6,7,8,11,16,17,18,29,47,57,66,111};

    
    // Removing manually all bad cells per plane - h0=T3 / h1=T4 / h2=T1 -> Check lookuptable as this is the order of TRBs
    Int_t UnactiveCells_h0[32]={7,61,62,63,64,65,66,73,74,75,76,77,78,82,85,86,87,88,89,90,97,98,99,100,101,102,109,110,111,112,113,114};
    Int_t UnactiveCells_h1[3]={1,2,120};
    Int_t UnactiveCells_h2[1]={17};

     for(Int_t i=6;i<7;i++){
        for(Int_t j=0;j<3;j++){
            for(Int_t bi=1;bi<=h_active_cells[i][j]->GetNbinsX();bi++){
                if(j==0){for(Int_t bj=0;bj<32;bj++){h_active_cells[i][j]->SetBinContent(bi, UnactiveCells_h0[bj], 0);}}
                if(j==1){for(Int_t bj=0;bj<3;bj++){h_active_cells[i][j]->SetBinContent(bi, UnactiveCells_h1[bj], 0);}}
                if(j==2){for(Int_t bj=0;bj<1;bj++){h_active_cells[i][j]->SetBinContent(bi, UnactiveCells_h2[bj], 0);}}
            }
        }
    }
    
    /*
    // Setting all cells to 1 to start a new year table
    for(Int_t i=6;i<7;i++){  
        for(Int_t j=0;j<3;j++){
            for(Int_t bi=1;bi<=h_active_cells[i][j]->GetNbinsX();bi++){
                for(Int_t bj=1;bj<120;bj++){
                    h_active_cells[i][j]->SetBinContent(bi, bj, 1);
		}
	    }
        }
    }
    */
    fOut->cd();

    for(Int_t j=0;j<3;j++) {
        h_active_cells[0][j]->Write();
        h_active_cells[1][j]->Write();
        h_active_cells[2][j]->Write();
        h_active_cells[3][j]->Write();
	h_active_cells[4][j]->Write();
	h_active_cells[5][j]->Write();
	h_active_cells[6][j]->Write();
        //h_active_cells[7][j]->Write();
	//cout<<"j "<<j<<" "<<

    }

    fOut->Write();
    fOut->Close();

}






















