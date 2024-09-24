// TODO: study the purity of the selected sample

#include "utils.h"

void Purity(const char *inFile = "LHC180_var_0", const char* outFile = "outPID"){
  TFile *_file0 = TFile::Open(Form("%s/%s.root", kResDir, inFile));
  TFile *fout = TFile::Open(Form("%s.root", outFile), "recreate");
  TH3D *outerPID[2]{nullptr};
  TH1D *purity[2]{nullptr};
  float ptBins[kNBinsPt + 1];
  for (int iB = 0; iB < kNBinsPt + 1; ++iB){
    ptBins[iB] = kMinPt + kDeltaPt * iB;
  }
  for (int iC{0}; iC < 2; ++iC){
    outerPID[iC] = (TH3D*)_file0->Get(Form("h%sOuterPID_0", kAntiMatterLabel[iC]));
    purity[iC] = new TH1D(Form("Purity_%s", kAntiMatterLabel[iC]), ";#it{p}_{T} (GeV/#it{c});S / (S + B)", kNBinsPt, ptBins);
    for (int iP{1}; iP < kNBinsPt; ++iP){
      TH1D* outerPID_proj = (TH1D*)outerPID[iC]->ProjectionZ(Form("outerPID_%s_%.2f_%.2f", kAntiMatterLabel[iC], outerPID[iC]->GetYaxis()->GetBinLowEdge(iP), outerPID[iC]->GetYaxis()->GetBinUpEdge(iP)), 1, kNCentBins, iP, iP);
      TF1 fit("fit", "gaus(0)+expo(3)");
      fit.SetParLimits(0, 0., 1.e7);
      fit.SetParLimits(1, -2., 2.);
      fit.SetParLimits(2, 0.1, 2.);
      fit.SetParLimits(3, -10., 10.);
      fit.SetParLimits(4, -10., 10.);
      outerPID_proj->Fit("fit", "ML+");
      TF1 sig("sig", "gaus");
      TF1 bkg("bkg", "expo");
      sig.SetParameter(0, fit.GetParameter(0));
      sig.SetParameter(1, fit.GetParameter(1));
      sig.SetParameter(2, fit.GetParameter(2));
      bkg.SetParameter(0, fit.GetParameter(3));
      bkg.SetParameter(1, fit.GetParameter(4));
      double purity_ = (sig.Integral(-3., 3.)) / (sig.Integral(-3., 3.) + bkg.Integral(-3., 3.));
      purity[iC]->SetBinContent(iP, purity_);
      fout->cd();
      outerPID_proj->Write();
      std::cout << "Purity = " << purity_ << std::endl;
    }
    fout->cd();
    purity[iC]->Write();
  }
  fout->Close();
}