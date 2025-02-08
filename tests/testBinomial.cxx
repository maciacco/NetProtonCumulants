void testBinomial(){
  TFile *f = TFile::Open("NetProtonCumulants/results/MB/LHC22f30_var_364.root");
  TFile *fo = TFile::Open("CheckBinom.root", "recreate");
  TH3D* hGR3D = (TH3D*)f->Get("hMGenRecProton_364");
  for (int iG{2}; iG < 7; ++iG) {
    TH1D* hR = (TH1D*)hGR3D->ProjectionZ("hRec", 1, 7, iG, iG);
    TH1D* hR_shift = new TH1D("hR_shift", Form("N_gen = %d;#it{N}_{rec};#it{P}(#it{N}_{rec})", iG - 1), 8, -0.5, 7.5);
    hR_shift->SetLineWidth(2);
    hR_shift->SetLineColor(kBlack);
    for (int iB{1}; iB < hR_shift->GetNbinsX() + 1; ++iB) {
      hR_shift->SetBinContent(iB, hR->GetBinContent(iB));
      hR_shift->SetBinError(iB, hR->GetBinError(iB));
    }
    hR_shift->Scale(1./hR_shift->Integral());
    TF1 binom("binom", "ROOT::Math::binomial_pdf(x+0.5,[1],[0])", -0.5, 7.5);
    //binom.SetNpx(10000);
    //binom.SetParLimits(1, 0, 1000.);
    binom.SetParLimits(1, 0., 1.);
    binom.SetParameter(0, (iG - 1) + 0.01);
    if (iG == 6) binom.SetParLimits(0, iG - 1, iG - 1 + 0.99);
    hR_shift->Fit("binom", "M+");
    //TF1 binom1("binom1", "TMath::Binomial(TMath::Nint([1]), TMath::Nint(x)) * TMath::Power([0], TMath::Nint(x)) * TMath::Power(1-[0], TMath::Nint([1])-TMath::Nint(x))", -0.5, 6.5);
    //binom1.SetParameter(0, binom.GetParameter(0));
    //binom1.SetParameter(1, binom.GetParameter(1));
    fo->cd();
    //binom1.Write();
    hR_shift->Write();
  }
  fo->Close();
  f->Close();
}
