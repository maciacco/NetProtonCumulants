void RatioHist(){
  TFile *_file0 = TFile::Open("out_sys_18_finalBinning_k6.root");
  TGraphErrors* rec = (TGraphErrors*)_file0->Get("g_364");
  TGraphErrors* gen = (TGraphErrors*)_file0->Get("g_gen_364");
  double multBins[] = {0, 10, 20, 30, 40, 50, 70, 100};
  TH1D *hgen = new TH1D("hgen", ";V0M Multiplicity (%);#kappa_{6}", 7, multBins);
  TH1D *hrec = new TH1D("hrec", ";V0M Multiplicity (%);#kappa_{6} Rec - Gen", 7, multBins);
  for (int i{1}; i < 8; ++i) {
    hgen->SetBinContent(i, gen->GetPointY(i - 1));
    hgen->SetBinError(i, 0., gen->GetErrorY(i - 1));
    hrec->SetBinContent(i, rec->GetPointY(i - 1));
    hrec->SetBinError(i, 0., rec->GetErrorY(i - 1));
    hrec->SetBinContent(i, hrec->GetBinContent(i) - hgen->GetBinContent(i));
    hrec->SetBinError(i, sqrt(pow(hrec->GetBinError(i), 2.)) - pow(hgen->GetBinError(i), 2.));
  }

  TFile *out = TFile::Open("ratio_k5.root", "recreate");
  out->cd();
  hrec->Write();
  out->Close();
}
