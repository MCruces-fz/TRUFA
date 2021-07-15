{

gSystem->Load("libGraf");
gSystem->Load("libtunpacker.so");
printf("Unpacker for stand alone TRB loaded\n");
gStyle->SetPalette(55); 
gROOT->SetStyle("Plain"); 
//gROOT->ProcessLine(".L makeTdcHistos.C");
//	   gROOT->ProcessLine(".L analysis.C");
//	   gROOT->ProcessLine(".L DriftTime.C");
//	   gROOT->ProcessLine(".L Planes.C");
}
