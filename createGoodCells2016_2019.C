#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TStyle.h"



#define nyears 7
#define ndetectors 3

TRandom3 randN;

void cleanPeriodThresholdY(TH2D* hIn, Float_t con) {
    TH1D* proj = hIn->ProjectionX(Form("%i",randN.Rndm()));
    for(int i=1;i<=proj->GetNbinsX();i++) {
        if(proj->GetBinContent(i)<120.*con) {
            for(int j=1;j<=120;j++) {
                hIn->SetBinContent(i,j,0);
                hIn->SetBinError(i,j,0);
            }
        }
    }
}

void cleanPeriodThresholdX(TH2D* hIn) {
    TH1D* proj;
    for(int i=0;i<=365;i++) {
        Int_t binfirst = hIn->GetXaxis()->FindBin(i+0.);
        Int_t binlast  = hIn->GetXaxis()->FindBin(i+1.);
        proj = hIn->ProjectionY("",binfirst,binlast-1);
        for(int j=1;j<=120;j++) {
            if(proj->GetBinContent(j)<24) {
                for(int bi=binfirst;bi<binlast;bi++) {
                    hIn->SetBinContent(bi,j,0);
                    hIn->SetBinError(bi,j,0);
                }
            }
        }
    }
}


void cleanCellsTimeRange(TH2D* hIn, Int_t cell, Float_t tmin, Float_t tmax ) {
    Int_t binfirst = hIn->GetXaxis()->FindBin(tmin);
    Int_t binlast  = hIn->GetXaxis()->FindBin(tmax);
    for(int i=binfirst;i<=binlast;i++) {
        hIn->SetBinContent(i,cell,0);
        hIn->SetBinError(i,cell,0);
    }
}

TH2D* createCellsBinary(TH2D* hIn0, TH2D* hIn1, TH2D* hIn2, TString name) {

    TH2D* h = (TH2D*)hIn0->Clone(name.Data());
    for(int i=1;i<=hIn0->GetNbinsX();i++) {
        for(int j=1;j<=hIn0->GetNbinsY();j++) {
            if(hIn0->GetBinContent(i,j)!=0 && hIn0->GetBinContent(i,j)!=0 && hIn0->GetBinContent(i,j)!=0 ) {
                h->SetBinContent(i,j,1);
                h->SetBinError(i,j,0.);
            }else{
                h->SetBinContent(i,j,0);
                h->SetBinError(i,j,0.);
            }
        }
    }
    return h;
}




TH1D* makeMostProbableProfile(TH2D* hIn,TString name, Float_t tmin, Float_t tmax) {

    Int_t binfirst = hIn->GetXaxis()->FindBin(tmin);
    Int_t binlast  = hIn->GetXaxis()->FindBin(tmax);

    TH2D* hIn_clone = (TH2D*)hIn->Clone("temp");

    TH2D h_pf2("hpf2d","hpf2d",120,0,120,100,0, 4.);

    TH1D* h_pf = new TH1D("hpf","hpf",120,0,120);


    TH1D* h_mean = new TH1D(name,Form("%s;cell;mean",name.Data()),120,0,120);
    TH1D* hAll_x = hIn->ProjectionX("_px",0,-1);
    TH1D* hAll_y = hIn->ProjectionY("_py",0,-1);


    TH1D* proj;
    // First loop to get the mean and std
    //for(int i=1;i<=hAll_x->GetNbinsX();i++) {
    for(int i=binfirst;i<=binlast;i++) {
        proj = hIn->ProjectionY("",i,i);
        if(proj->Integral()<100) continue;
        proj->Scale(100./proj->Integral());
        for(int nc = 0; nc<120; nc++) {
            if( proj->GetBinContent(nc+1)>0. ) {
                h_pf2.Fill(nc+0.5, proj->GetBinContent(nc+1) );
            }
        }
    }

    TF1 fg("fg","gaus",0.,5);
    TH1D* pg;
    for(int i=1;i<=120;i++) {

        //Fit bin by bin the projections!
        pg = h_pf2.ProjectionY("",i,i);
        if(pg->Integral()>50) {
            fg.SetParameters(pg->GetMaximum(),pg->GetMean(),pg->GetRMS());
            //pg->Fit(&fg,"NQWW");
            //pg->Fit(&fg,"NQWWR","",
            //        fg.GetParameter(1)-fg.GetParameter(2)*3.,
            //        fg.GetParameter(1)+fg.GetParameter(2)*3.);

            h_pf->SetBinContent(i,pg->GetMean());
            h_pf->SetBinError(i,pg->GetRMS());
        } else {
            h_pf->SetBinContent(i,0);
            h_pf->SetBinError(i,0);
        }
    }

    // First filtered loop
    for(int i=binfirst;i<=binlast;i++) {
        proj = hIn->ProjectionY("",i,i);
        if(proj->Integral()<100) {
            for(int nc = 0; nc<120; nc++) {
                hIn_clone -> SetBinContent(i,nc+1,0.);
                hIn -> SetBinContent(i,nc+1,0.);
                hIn -> SetBinError(i,nc+1,0.);
            }
            continue;
        }
        // Scale by the expected average: Integral * Sum/100.
        Double_t integral = proj->Integral();
        Double_t scaling = 0.;
        for(int nc = 1; nc<=120; nc++) {
            // If the cell is disconnected
            if( proj->GetBinContent(nc) == 0 ) {
                scaling += h_pf->GetBinContent(nc)/100.;
            }
        }
        //if(scaling>0.3) cout<<"sc problem "<<scaling<<endl;
        integral*=1./(1.-scaling);
        proj->Scale(100. / integral );
        for(int nc = 0; nc<120; nc++) {
            if( fabs(proj->GetBinContent(nc+1) - h_pf->GetBinContent(nc+1) ) < 4.5 * h_pf->GetBinError(nc+1) ) {
                // everything is OK
            }
            else {
                hIn_clone -> SetBinContent(i,nc+1,0);
                hIn -> SetBinContent(i,nc+1,0);
                hIn -> SetBinError(i,nc+1,0);
            }
        }
    }


    // Second loop
    // Reset the histos
    h_pf2.Reset();

    for(int i=binfirst;i<=binlast;i++) {
        proj = hIn_clone->ProjectionY("",i,i);
        if(proj->Integral()<100) continue;

        // Scale by the expected average: Integral * Sum/100.
        Double_t integral = proj->Integral();
        Double_t scaling = 0.;
        for(int nc = 1; nc<=120; nc++) {
            // If the cell is disconnected
            if( proj->GetBinContent(nc) == 0 ) {
                scaling += h_pf->GetBinContent(nc)/100.;
            }
        }
        //if(scaling>0.3) cout<<"sc problem2 "<<scaling<<endl;
        integral*=1./(1.-scaling);
        proj->Scale(100. / integral );
        for(int nc = 0; nc<120; nc++) {
            if( proj->GetBinContent(nc+1)>0. ) {
                h_pf2.Fill(nc+0.5, proj->GetBinContent(nc+1) );
            }
        }
    }
    for(int i=1;i<=120;i++) {
        //Fit bin by bin the projections!
        pg = h_pf2.ProjectionY("",i,i);
        if(pg->Integral()>50) {
            fg.SetParameters(pg->GetMaximum(),pg->GetMean(),pg->GetRMS());
//        pg->Fit(&fg,"NQWW");
//        pg->Fit(&fg,"NQWWR","",
//                fg.GetParameter(1)-fg.GetParameter(2)*3.,
//                fg.GetParameter(1)+fg.GetParameter(2)*3.);
//        h_pf->SetBinContent(i,fg.GetParameter(1));
//        h_pf->SetBinError(i,fg.GetParameter(2));
            h_pf->SetBinContent(i,pg->GetMean());
            h_pf->SetBinError(i,pg->GetRMS());
        } else {
            h_pf->SetBinContent(i,0);
            h_pf->SetBinError(i,0);
        }
    }

    for(int i=binfirst;i<=binlast;i++) {
        proj = hIn_clone->ProjectionY("",i,i);
        if(proj->Integral()<100) continue;
        // Scale by the expected average: Integral * Sum/100.

        Double_t integral = proj->Integral();
        Double_t scaling = 0.;
        for(int nc = 1; nc<=120; nc++) {
            // If the cell is disconnected
            if( proj->GetBinContent(nc) == 0 ) {
                scaling += h_pf->GetBinContent(nc)/100.;
            }
        }

        //if(scaling>0.3) cout<<"sc problem3 "<<scaling<<endl;
        integral*=1./(1.-scaling);
        proj->Scale(100. / integral );

        for(int nc = 0; nc<120; nc++) {
            if( fabs(proj->GetBinContent(nc+1) - h_pf->GetBinContent(nc+1) ) < 4.5 * h_pf->GetBinError(nc+1) ) {
                // everything is OK
            }
            else {
                hIn_clone -> SetBinContent(i,nc+1,0);
                hIn -> SetBinContent(i,nc+1,0);
                hIn -> SetBinError(i,nc+1,0);
            }
        }
    }
    // Third loop
    // Reset the histos
    h_pf2.Reset();

    for(int i=binfirst;i<=binlast;i++) {
        proj = hIn_clone->ProjectionY("",i,i);
        if(proj->Integral()<100) continue;
        // Scale by the expected average: Integral * Sum/100.

        Double_t integral = proj->Integral();
        Double_t scaling = 0.;
        for(int nc = 1; nc<=120; nc++) {
            // If the cell is disconnected
            if( proj->GetBinContent(nc) == 0 ) {
                scaling += h_pf->GetBinContent(nc)/100.;
            }
        }
        //if(scaling>0.5) cout<<"sc problem4 "<<scaling<<endl;

        integral*=1./(1.-scaling);
        proj->Scale(100. / integral );
        for(int nc = 0; nc<120; nc++) {
            if( proj->GetBinContent(nc+1)>0. ) {
                h_pf2.Fill(nc+0.5, proj->GetBinContent(nc+1) );
            }
        }
    }
    for(int i=1;i<=120;i++) {
        //Fit bin by bin the projections!
        pg = h_pf2.ProjectionY("",i,i);
        if(pg->Integral()>50) {
            fg.SetParameters(pg->GetMaximum(),pg->GetMean(),pg->GetRMS());
            pg->Fit(&fg,"NQ");
//            pg->Fit(&fg,"NQ","",
//                    fg.GetParameter(1)-fg.GetParameter(2)*3.,
//                    fg.GetParameter(1)+fg.GetParameter(2)*3.);
//            h_pf->SetBinContent(i,fg.GetParameter(1));
//            h_pf->SetBinError(i,fg.GetParameter(2));
            h_pf->SetBinContent(i,pg->GetMean());
            h_pf->SetBinError(i,pg->GetRMS());

        } else {
            h_pf->SetBinContent(i,0);
            h_pf->SetBinError(i,0);
        }

    }
    for(int i=binfirst;i<=binlast;i++) {
        proj = hIn->ProjectionY("",i,i);
        if(proj->Integral()<100) {
            for(int nc = 1; nc<=120; nc++) {
                hIn -> SetBinContent(i,nc+1,0);
                hIn -> SetBinError(i,nc+1,0);
            }
            continue;
        }
        // Scale by the expected average: Integral * Sum/100.
        Double_t integral = proj->Integral();
        Double_t scaling = 0.;
        for(int nc = 1; nc<=120; nc++) {
            // If the cell is disconnected
            if( proj->GetBinContent(nc) == 0 ) {
                scaling += h_pf->GetBinContent(nc)/100.;
            }
        }

        if(scaling>0.2) {
            for(int nc = 1; nc<=120; nc++) {
                hIn -> SetBinContent(i,nc+1,0);
                hIn -> SetBinError(i,nc+1,0);
            }
            //cout<<"sc problem too few in "<<scaling<<endl;
        }
        integral*=1./(1.-scaling);
        proj->Scale(100. / integral );


        for(int nc = 0; nc<120; nc++) {
            if( fabs(proj->GetBinContent(nc+1) - h_pf->GetBinContent(nc+1) ) < 4. * h_pf->GetBinError(nc+1) ) {
                // everything is OK
            }
            else {
                hIn -> SetBinContent(i,nc+1,0);
                hIn -> SetBinError(i,nc+1,0);
            }
        }
    }


    return h_pf;


}

void clean_with_daq(TH2D* h, TH1D* hdaq, Float_t thresholdMin,  Float_t thresholdMax) {
    for(int i=1;i<=h->GetNbinsX();i++) {
        if( hdaq -> GetBinContent(i) > thresholdMax ||
           hdaq -> GetBinContent(i) < thresholdMin) {
            for(int j=1;j<=120;j++) {
                h->SetBinContent(i,0);
                h->SetBinError(i,0);
            }
        }
    }
}

void createGoodCells2016_2019(){

    // Periods of data taking

    //std::vector<Int_t <Int_t>> periods;
    //std::vector<Int_t> temp;



    // Periods
    // 2015 :full
    //      cut values for cells in the period
    //      h0 cut  0.07 / 0.08
    //      h1 cut  0.07 / 0.05
    //      h2 cut  0.24 / 0.24
    // 2016 :full
    //      cut values for cells in the period
    //      h0 cut  0.05 / 0.05
    //      h1 cut  0.06 / 0.05
    //      h2 cut  0.08 / 0.07
    // 2017 :0 - 214
    //      cut values for cells in the period
    //      h0 cut
    //      h1 cut
    //      h2 cut
    // 2017 :268 - 366
    //      cut values for cells in the period
    //      h0 cut
    //      h1 cut
    //      h2 cut
    // 2018 : 0 - 204.4
    //      cut values for cells in the period
    //      h0 cut
    //      h1 cut
    //      h2 cut
    // 2018 : 207 - 220.4
    //      cut values for cells in the period
    //      h0 cut
    //      h1 cut
    //      h2 cut
    // 2018 : 220.4 - 366
    //      cut values for cells in the period
    //      h0 cut
    //      h1 cut
    //      h2 cut
    // 2019 : full
    //      cut values for cells in the period
    //      h0 cut
    //      h1 cut
    //      h2 cut
    // 2020 : full
    //      cut values for cells in the period
    //      h0 cut
    //      h1 cut
    //      h2 cut




    // open all the original files
    //int nyears = 6;
    TFile* fIn[nyears];
    for(int i=0;i<nyears;i++) {
        fIn[i] = TFile::Open(Form("cells_dst_%i.root",2015+i));
    }
    //for(int i=0;i<4;i++) {
    //    //fIn[0] = TFile::Open(Form("dst_20%i_full.root",16+i));
    //    cout<<fIn[i]<<Form(" dst_20%i_full.root",16+i)<<(TH2D*)fIn[i]->Get("h_cell_counts_h0")<<endl;
    //}

    TH2D* h_good_cells[nyears][ndetectors][3];
    TH2D* h_daq_cells[nyears];

    TH1D* h_daq[nyears][4];
    TString name_daq[4] = {"h_daq","h_daq_total","h_daq_seconds","h_daq_sync"};
    for(int i=0;i<nyears;i++) {
        for(int j=0;j<4;j++) {
            h_daq[i][j] = (TH1D*)fIn[i]->Get( Form("%s",name_daq[j].Data()) );
            // Rebin in order to acoomodate to the active cells histograms
            if(h_daq[i][j]) {
                h_daq[i][j]->Sumw2();
                h_daq[i][j]->Rebin(6);
                h_daq[i][j]->Scale(1./6);
            }
        }
        // remove all intervals with not sync
        for(int bi=1;bi<= h_daq[i][0]-> GetNbinsX();bi++) {
            if(h_daq[i][3] -> GetBinContent(bi)!=0) {
                for(int j=0;j<3;j++) {
                    h_daq[i][j] -> SetBinContent(bi,0);
                }

            }
        }
    }





    TCanvas* can_daq_total[nyears];// = new TCanvas("can_daq_total"),"can_daq_total"),)


    for(int i=0;i<nyears;i++) {
        can_daq_total[i] = new TCanvas(Form("can_daq_total_%i",2015+i),
                                       Form("can_daq_total_%i",2015+i));
        can_daq_total[i] ->Divide(2,2,0.001,0.001);

        h_daq[i][1] -> Divide(h_daq[i][2]);
        h_daq[i][1] -> SetMaximum(200);

        can_daq_total[i] -> cd(1);
        h_daq[i][0] -> DrawClone();
        can_daq_total[i] -> cd(2);
        h_daq[i][1] -> DrawClone();
        can_daq_total[i] -> cd(3);
        h_daq[i][2] -> DrawClone();
        can_daq_total[i] -> cd(4);
        h_daq[i][3] -> DrawClone();
        //for(int j=0;j<4;j++) {
            //
            // correct for the number of active seconds.
            //
            //
        //}
    }


    TCanvas* can_good_cells[nyears];
    TH1D* proj;
    for(int i=0;i<nyears;i++) {
        can_good_cells[i] = new TCanvas(Form("can_good_cells_%i",2015+i),
                                        Form("can_good_cells_%i",2015+i));

        can_good_cells[i] -> Divide(2,1);

        for(int j=0;j<3;j++) {

            h_good_cells[i][j][0] = (TH2D*)fIn[i]->Get(Form("h_cell_counts_h%i",j)   );
            h_good_cells[i][j][1] = (TH2D*)fIn[i]->Get(Form("h_cell_counts_n1_h%i",j));
            h_good_cells[i][j][2] = (TH2D*)fIn[i]->Get(Form("h_cell_counts_n2_h%i",j));

            h_good_cells[i][j][0] -> SetName( Form("h_cell_counts_y%i_h%i",i+2015,j)    );
            h_good_cells[i][j][1] -> SetName( Form("h_cell_counts_y%i_n1_h%i",i+2015,j) );
            h_good_cells[i][j][2] -> SetName( Form("h_cell_counts_y%i_n2_h%i",i+2015,j) );




            for(int k=0;k<3;k++) {
                if(i<3)
                    clean_with_daq(h_good_cells[i][j][k],h_daq[i][1], 50. , 100. );
                else
                    clean_with_daq(h_good_cells[i][j][k],h_daq[i][1], 60. , 120. );
            }

            proj = h_good_cells[i][j][0] -> ProjectionX("");

            can_good_cells[i] -> cd(1);
            h_good_cells[i][j][0] -> DrawClone("colz");

            can_good_cells[i] -> cd(2);
            proj->Divide(h_daq[i][2]);
            proj->DrawClone();

        }
        // divide by daq;
        h_daq_cells[i] = (TH2D*) h_good_cells[i][0][0] -> Clone(Form("h_daq_cells_y%i",i));
        h_daq_cells[i] -> SetNameTitle(Form("h_daq_cells_y%i",i),
                                       Form("h_daq_cells_y%i",i));
        h_daq_cells[i] -> Reset();
        for(int bi=1;bi<=h_daq_cells[i]->GetNbinsX();bi++) {
            for(int bj =1;bj<=120;bj++) {
                h_daq_cells[i] -> SetBinContent( bi,bj,h_daq[i][2]->GetBinContent(bi) );
                h_daq_cells[i] -> SetBinError( bi,bj,h_daq[i][2]->GetBinError(bi) );
            }
        }

        for(int j=0;j<3;j++) {
            for(int k=0;k<3;k++) {
                //h_good_cells[i][j][k] -> Divide( h_daq_cells[i] );
            }
        }
    }






    // h_good_cells[0][0] = (TH2D*)fIn[2]->Get("h_cell_counts_h2");


    TF1* fpol0 = new TF1("fpol0","pol0",0,400);
    /*
    for(int i=0;i<4;i++) {
        for(int j=0;j<3;j++) {
            if(h_daq[i][0])
            clean_with_daq(h_good_cells[i][j], h_daq[i][2], 500. , 620. );

            if(h_daq[i][2])
            clean_with_daq(h_good_cells[i][j], h_daq[i][2], 500. , 620. );

            if(h_daq[i][1]) {
                h_daq[i][1] -> Fit(fpol0,"NQRW");
                clean_with_daq(h_good_cells[i][j], h_daq[i][1],
                               fpol0->GetParameter(0)*0.65,
                               fpol0->GetParameter(0)*1.35);
            }

            if(h_daq[i][3])
            clean_with_daq(h_good_cells[i][j], h_daq[i][3], -1 , 50 );

        }
        }*/
    TH1D* h_profile[9][3][3];
    TH1D* h_rms_cell_prof[9][3][3];

    // Periods
    // 2015 :full
    //      cut values for cells in the period
    //      h0 cut  0.07 / 0.08
    //      h1 cut  0.07 / 0.05
    //      h2 cut  0.24 / 0.24
    Double_t thresCut0[3][2] = {{0.07,0.08},{0.07,0.05},{0.24,0.24}};
    int i = 0;
    for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++) {
            h_profile[0][j][k] = makeMostProbableProfile(h_good_cells[0][j][k],Form("profile_y%i_d%i_n%i",0,j,k),0.,366.0);
            h_rms_cell_prof[0][j][k] = new TH1D(Form("h_rms_cells_prof_y%i_d%i_n%i",0,j,k),
                                                Form("h_rms_cells_prof_y%i_d%i_n%i",0,j,k),
                                                250,0,1.);
            for(int bi=1;bi<=120;bi++) {
                h_rms_cell_prof[i][j][k] -> Fill(h_profile[i][j][k]->GetBinError(bi));
                if(k<2 && thresCut0[j][k]<h_profile[i][j][k]->GetBinError(bi) ) {
                    cleanCellsTimeRange(h_good_cells[i][j][k], bi , 0.,366. );
                }
            }
        }
    }
    // 2016 :full
    //      cut values for cells in the period
    //      h0 cut  0.05 / 0.05
    //      h1 cut  0.06 / 0.05
    //      h2 cut  0.08 / 0.07
    Double_t thresCut1[3][2] = {{0.05,0.05},{0.06,0.05},{0.08,0.07}};
    i=1;
    int yearblock = 1;
    for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++) {
            h_profile[i][j][k] = makeMostProbableProfile(h_good_cells[yearblock][j][k],Form("profile_y%i_d%i_n%i",i,j,k),0.,366.0);
            h_rms_cell_prof[i][j][k] = new TH1D(Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                250,0,1.);
            for(int bi=1;bi<=120;bi++) {
                h_rms_cell_prof[i][j][k] -> Fill(h_profile[i][j][k]->GetBinError(bi));
                if(k<2 && thresCut1[j][k]<h_profile[i][j][k]->GetBinError(bi) ) {
                    cleanCellsTimeRange(h_good_cells[yearblock][j][k], bi , 0.,366. );
                }
            }
        }
    }
    // 2017 :0 - 214
    //      cut values for cells in the period    std>0
    //      h0 cut 0.07 / 0.07
    //      h1 cut 0.07 / 0.07
    //      h2 cut 0.4 / 0.35
    Double_t thresCut2[3][2] = {{0.07,0.07},{0.07,0.07},{0.4,0.35}};
    i=2;
    yearblock = 2;
    for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++) {
            h_profile[i][j][k] = makeMostProbableProfile(h_good_cells[yearblock][j][k],Form("profile_y%i_d%i_n%i",i,j,k),0.,214.0);
            h_rms_cell_prof[i][j][k] = new TH1D(Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                250,0,1.);
            for(int bi=1;bi<=120;bi++) {
                h_rms_cell_prof[i][j][k] -> Fill(h_profile[i][j][k]->GetBinError(bi));
                if(k<2 && thresCut2[j][k]<h_profile[i][j][k]->GetBinError(bi) ) {
                    cleanCellsTimeRange(h_good_cells[yearblock][j][k], bi , 0.,214.0 );
                }
            }
        }
    }
    // 2017 :268 - 366
    //      cut values for cells in the period
    //      h0 cut 0.42 / 0.4
    //      h1 cut 0.1 /  0.1
    //      h2 cut 0.14 / 0.12
    Double_t thresCut3[3][2] = {{0.42,0.4},{0.1,0.1},{0.14,0.12}};
    i=3;
    yearblock = 2;
    for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++) {
            h_profile[i][j][k] = makeMostProbableProfile(h_good_cells[yearblock][j][k],Form("profile_y%i_d%i_n%i",i,j,k),268,366.0);
            h_rms_cell_prof[i][j][k] = new TH1D(Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                250,0,1.);
            for(int bi=1;bi<=120;bi++) {
                h_rms_cell_prof[i][j][k] -> Fill(h_profile[i][j][k]->GetBinError(bi));
                if(k<2 && thresCut3[j][k]<h_profile[i][j][k]->GetBinError(bi) ) {
                    cleanCellsTimeRange(h_good_cells[yearblock][j][k], bi , 268 , 366.0 );
                }
            }
        }
    }
    // 2018 : 0 - 204.4
    //      cut values for cells in the period
    //      h0 cut  0.2 / 0.2
    //      h1 cut  0.15 / 0.15
    //      h2 cut  0.3 / 0.3
    Double_t thresCut4[3][2] = {{0.2,0.2},{0.15,0.15},{0.3,0.3}};
    i=4;
    yearblock = 3;
    for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++) {
            h_profile[i][j][k] = makeMostProbableProfile(h_good_cells[yearblock][j][k],Form("profile_y%i_d%i_n%i",i,j,k),0.,204.4);
            h_rms_cell_prof[i][j][k] = new TH1D(Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                250,0,1.);
            for(int bi=1;bi<=120;bi++) {
                h_rms_cell_prof[i][j][k] -> Fill(h_profile[i][j][k]->GetBinError(bi));
                if(k<2 && thresCut4[j][k]<h_profile[i][j][k]->GetBinError(bi) ) {
                    cleanCellsTimeRange(h_good_cells[yearblock][j][k], bi , 0. , 204.4 );
                }
            }
        }
    }
    // 2018 : 207 - 220.4
    //      cut values for cells in the period
    //      h0 cut  0.06 / 0.06
    //      h1 cut  0.05 / 0.04
    //      h2 cut  0.04 / 0.05
    Double_t thresCut5[3][2] = {{0.06,0.06},{0.05,0.04},{0.04,0.05}};
    i=5;
    yearblock = 3;
    for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++) {
            h_profile[i][j][k] = makeMostProbableProfile(h_good_cells[yearblock][j][k],Form("profile_y%i_d%i_n%i",i,j,k),207.0,220.4);
            h_rms_cell_prof[i][j][k] = new TH1D(Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                250,0,1.);
            for(int bi=1;bi<=120;bi++) {
                h_rms_cell_prof[i][j][k] -> Fill(h_profile[i][j][k]->GetBinError(bi));
                if(k<2 && thresCut5[j][k]<h_profile[i][j][k]->GetBinError(bi) ) {
                    cleanCellsTimeRange(h_good_cells[yearblock][j][k], bi , 207.0 , 220.4 );
                }
            }
        }
    }
    // 2018 : 220.4 - 366
    //      cut values for cells in the period
    //      h0 cut  0.2 / 0.15
    //      h1 cut  0.12 / 0.12
    //      h2 cut  0.25 / 0.22
    Double_t thresCut6[3][2] = {{0.2,0.15},{0.12,0.12},{0.25,0.22}};
    i=6;
    yearblock = 3;
    for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++) {
            h_profile[i][j][k] = makeMostProbableProfile(h_good_cells[yearblock][j][k],Form("profile_y%i_d%i_n%i",i,j,k),220.4,366.);
            h_rms_cell_prof[i][j][k] = new TH1D(Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                250,0,1.);
            for(int bi=1;bi<=120;bi++) {
                h_rms_cell_prof[i][j][k] -> Fill(h_profile[i][j][k]->GetBinError(bi));
                if(k<2 && thresCut6[j][k]<h_profile[i][j][k]->GetBinError(bi) ) {
                    cleanCellsTimeRange(h_good_cells[yearblock][j][k], bi , 220.4, 366.0);
                }
            }
        }
    }
    // 2019 : full
    //      cut values for cells in the period
    //      h0 cut  0.1/ 0.1
    //      h1 cut  0.08 / 0.08
    //      h2 cut  0.1 / 0.1
    Double_t thresCut7[3][2] = {{0.1,0.1},{0.08,0.08},{0.1,0.1}};
    i=7;
    yearblock = 4;
    for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++) {
            h_profile[i][j][k] = makeMostProbableProfile(h_good_cells[yearblock][j][k],Form("profile_y%i_d%i_n%i",i,j,k),0.   ,366.);
            h_rms_cell_prof[i][j][k] = new TH1D(Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                250,0,1.);
            for(int bi=1;bi<=120;bi++) {
                h_rms_cell_prof[i][j][k] -> Fill(h_profile[i][j][k]->GetBinError(bi));
                if(k<2 && thresCut7[j][k]<h_profile[i][j][k]->GetBinError(bi) ) {
                    cleanCellsTimeRange(h_good_cells[yearblock][j][k], bi , 0, 366.0);
                }
            }
        }
    }
    // 2020 : full
    //      cut values for cells in the period
    //      h0 cut  0.1 / 0.1
    //      h1 cut  0.08 / 0.08
    //      h2 cut  0.1 / 0.09
    Double_t thresCut8[3][2] = {{0.1,0.1},{0.08,0.08},{0.1,0.1}};
    i=8;
    yearblock = 5;
    for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++) {
            h_profile[i][j][k] = makeMostProbableProfile(h_good_cells[yearblock][j][k],Form("profile_y%i_d%i_n%i",i,j,k),0.   ,366.);
            h_rms_cell_prof[i][j][k] = new TH1D(Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                Form("h_rms_cells_prof_y%i_d%i_n%i",i,j,k),
                                                250,0,1.);
            for(int bi=1;bi<=120;bi++) {
                h_rms_cell_prof[i][j][k] -> Fill(h_profile[i][j][k]->GetBinError(bi));
                if(k<2 && thresCut8[j][k]<h_profile[i][j][k]->GetBinError(bi) ) {
                    cleanCellsTimeRange(h_good_cells[yearblock][j][k], bi , 0, 366.0);
                }
            }
        }
    }

    /*
    TCanvas* can_rms_prof[9];
    for(int i=0;i<9;i++) {
        can_rms_prof[i] = new TCanvas(Form("can_rms_y_%i",i),Form("can_rms_y_%i",i));
        can_rms_prof[i] -> Divide(3,3,0.001,0.001);
        for(int j=0;j<3;j++) {
            for(int k=0;k<3;k++) {
                can_rms_prof[i] -> cd(j*3+k+1);
                h_rms_cell_prof[i][j][k] -> Draw();
            }
        }
    }
    */


    // Fully disconnect the cells with very large fluctuations!
    TH2D* h_clean_cells[nyears][3];
    for(int i=0;i<nyears;i++) {
        for(int j=0;j<3;j++) {
            h_clean_cells[i][j] = createCellsBinary(h_good_cells[i][j][0],h_good_cells[i][j][1],h_good_cells[i][j][2],Form("h_h%i_%i",j,i+2015));
            cleanPeriodThresholdX(h_clean_cells[i][j]);
            cleanPeriodThresholdY(h_clean_cells[i][j], 0.65);
        }
    }

    TCanvas* can_cells[nyears];
    for(int i=0;i<nyears;i++) {
        can_cells[i] = new TCanvas(Form("can_cells_total_%i",2015+i),
                                   Form("can_cells_total_%i",2015+i));

        can_cells[i] -> Divide(1,3);
        for(int j=0;j<3;j++) {
            can_cells[i]->cd(j+1);
            h_clean_cells[i][j] ->SetMaximum(1);
            h_clean_cells[i][j] ->SetMinimum(0);
            h_clean_cells[i][j] -> DrawClone("col");
        }
    }

    //
    // Make 2020 equal to 1 for the not measured interval!
    //
    for(int j=0;j<3;j++) {

        Int_t start_bin = h_clean_cells[5][j] -> GetXaxis() -> FindBin(97.0);
        cout<<"Bins: "<<start_bin<<endl;
        for(int bi=start_bin; bi<= h_clean_cells[5][j] -> GetNbinsX(); bi++) {
            for(int bj=1; bj<=120;bj++) {
                h_clean_cells[5][j] -> SetBinContent(bi,bj,1);
            }
        }
    }

    TFile* fOut = new TFile("clean_cells_2015_2021.root","recreate");
    for(int i=0;i<nyears;i++) {
        for(int j=0;j<3;j++) {
            h_clean_cells[i][j] -> Write();
        }
    }
    fOut->Close();


    //TCanvas* can1 = new TCanvas("can1","can1");
    //h_good_cells[0][0] ->Draw("colz");

    //TCanvas* can = new TCanvas();
    //h->DrawClone();

    /*
    TCanvas* can_clean[6][3];// = new TCanvas("can_daq_total"),"can_daq_total"),)
    for(int i=0;i<nyears;i++) {
        for(int j=0;j<3;j++) {
        can_clean[i][j] = new TCanvas(Form("can_clean_total_%i_%i",2015+i,j),
                                       Form("can_clean_total_%i_%i",2015+i,j));
        can_clean[i][j] ->Divide(3,2,0.001,0.001);
            for(int k=0;k<3;k++) {
                can_clean[i][j]->cd(k*2+1);
                h_good_cells[i][j][k] -> DrawClone("colz");
                can_clean[i][j]->cd(k*2+2);
                h_profile[i][j][k] -> DrawClone();
            }
        }
    }
    */




}
