#include "utils.h"

void Multiplicity(const char *inFile = "LHC180_var_0", const char* outFile = "Mult"){
  TFile *_file0 = TFile::Open(Form("%s/%s.root", kResDir, inFile));
  TFile *fout = TFile::Open(Form("%s/%s.root", kCalibDir, outFile), "recreate");

  TH1D *hNtrkl = (TH1D*)_file0->Get("hNtrkl_0");
  TH1D *hEstim = new TH1D("hEstim", ";#it{N}_{tracklets}^{0.7 < |#eta| < 1.2};Pecentile (%)", 100, 0, 100);
  double ntrkl_intg = hNtrkl->Integral(0, 100);
  for (int i{0}; i < 100; ++i){
    double intg_tmp = hNtrkl->Integral(i, 100);
    hEstim->SetBinContent(i + 1, intg_tmp / ntrkl_intg * 100.);
  }

  fout->cd();
  hEstim->Write();

  fout->Close();
}