void makeChargeReportFromCalFile(TString dir, TString filename){

    TString fullname = dir+filename;
    TFile* fIn = TFile::Open(fullname.Data());
    TCanvas* can = new TCanvas("canvas_charge","canvas",625,1000);
    can->Divide(2,4);

    can->SaveAs(Form("charge_report_%s.pdf[",filename.Data()) ,"pdf");

    TH1D* h;
    int counter=0;
    for(int i=0;i<3;i++) {
        for(int j=0;j<10;j++) {
            for(int k=0;k<12;k++) {
                counter++;
                can->cd(counter);
                h = (TH1D*)fIn->Get(Form("h_q_d%i_r%i_c%i",i,j,k));
                h->GetXaxis()->SetRangeUser(20,80);
                h->DrawClone();
                if(counter==8) {
                    can->SaveAs(Form("charge_report_%s.pdf",filename.Data()),"pdf");
                    counter=0;
                }
            }
        }
    }
    can->SaveAs(Form("charge_report_%s.pdf]",filename.Data()),"pdf");


}