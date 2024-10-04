// TODO: study the purity of the selected sample

#include "utils.h"

Double_t GausDExpExp(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  Double_t u = (xx - par[1]) / par[2];
  Double_t f = 0;
  if (u < par[3])
    f = TMath::Exp(-par[3] * (u - 0.5 * par[3]));
  else if (u <= par[4])
    f = TMath::Exp(-u * u * 0.5);
  else
    f = TMath::Exp(-par[4] * (u - 0.5 * par[4]));
  f = par[0] * f + TMath::Exp(-par[5] * u + par[6]);
  return f;
}

Double_t GausDExp(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  Double_t u = (xx - par[1]) / par[2];
  Double_t f = 0;
  if (u < par[3])
    f = TMath::Exp(-par[3] * (u - 0.5 * par[3]));
  else if (u <= par[4])
    f = TMath::Exp(-u * u * 0.5);
  else
    f = TMath::Exp(-par[4] * (u - 0.5 * par[4]));
  f = par[0] * f;
  return f;
}

void Purity(const char *inFile = "LHC21d30_var_364", const char* outFile = "outPID"){
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

      TF1 fit("fit", GausDExpExp, -10., 10., 7);
      fit.SetParLimits(0, 0., 1.e7);
      fit.SetParLimits(1, -2., 2.);
      fit.SetParLimits(2, 0.1, 2.);
      fit.SetParLimits(3, -2., .8);
      fit.SetParLimits(4, 0.8, 2.);
      fit.SetParLimits(5, -10., 10.);
      fit.SetParLimits(6, -10., 10.);
      outerPID_proj->Fit("fit", "MRL+", "", -5., 6.);
      TF1 sig("sig", GausDExp, -10., 10., 5);
      TF1 bkg("bkg", "TMath::Exp(-[0]*x + [1])");
      sig.SetParameter(0, fit.GetParameter(0));
      sig.SetParameter(1, fit.GetParameter(1));
      sig.SetParameter(2, fit.GetParameter(2));
      sig.SetParameter(3, fit.GetParameter(3));
      sig.SetParameter(4, fit.GetParameter(4));
      bkg.SetParameter(0, fit.GetParameter(5));
      bkg.SetParameter(1, fit.GetParameter(6));
      double purity_ = (sig.Integral(-3., 3., 1.e-7)) / (sig.Integral(-3., 3., 1.e-7) + bkg.Integral(-3., 3., 1.e-7));
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