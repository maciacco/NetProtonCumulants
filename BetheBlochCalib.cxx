#include "utils.h"

// ported from o2
constexpr double MassProton = 0.9832720894329;
Double_t alephBetheBloch(Double_t *x, Double_t *par)
{
   Float_t bg =x[0];
   Double_t beta = bg / std::sqrt(1. + bg * bg);
   Double_t aa = std::pow(beta, par[3]);
   Double_t bb = std::pow(1. / bg, par[4]);
   bb = std::log(par[2] + bb);
   Double_t f = (par[1] - aa - bb) * par[0] / aa;
   return f;
}

void BetheBlochCalib(const char *period = "LHC18_dat", const double low = 0.25, const double up = 2.){
  auto inFile = TFile::Open(Form("%s/AnalysisResults_%s.root", kDataDir, period));
  auto outFile = TFile::Open(Form("%s/out_file_%s.root", kCalibDir, period), "recreate");
  outFile->mkdir("fits");
  auto hTPC = dynamic_cast<TH2F*>(inFile->Get("ebye-maker/QA/tpcSignalPr"));
  // hTPC->RebinX(5);
  auto xmin = hTPC->GetXaxis()->FindBin(low);
  auto xmax = hTPC->GetXaxis()->FindBin(up);
  // std::cout << xmax - xmin << std::endl;
  TH1D* proj = new TH1D[xmax - xmin];
  TH1D hTPCPr("hTPCPr", ";#beta#gamma;d#it{E}/d#it{x}_{TPC} (a.u.)", xmax - 1, 0., up / MassProton);
  for (int iP{xmin}; iP < xmax; ++iP) {
    // if (hTPC->GetXaxis()->GetBinCenter(iP) > 3. && hTPC->GetXaxis()->GetBinCenter(iP) < 4.5) continue;
    TH1D* proj_tmp = dynamic_cast<TH1D*>(hTPC->ProjectionY(Form("h_%d", iP), iP, iP));
    proj_tmp->Copy(proj[iP - xmin]);
    double mean = proj[iP - xmin].GetBinCenter(proj[iP - xmin].GetMaximumBin());
    proj[iP - xmin].GetXaxis()->SetRangeUser(mean - 0.12 * mean, mean + 0.12 * mean);
    double rms = proj[iP - xmin].GetRMS();
    proj[iP - xmin].GetXaxis()->SetRangeUser(0, 1.e4);
    proj[iP - xmin].Fit("gaus", "MRQL+", "", mean - 2.5 * rms, mean + 2.5 * rms);
    auto fitFunc = proj[iP - xmin].GetFunction("gaus");
    outFile->cd("fits");
    proj[iP - xmin].Write();
    hTPCPr.SetBinContent(iP, fitFunc->GetParameter(1));
    hTPCPr.SetBinError(iP, fitFunc->GetParameter(2));
  }
  TF1 alephBB("alephBB", alephBetheBloch, low, up, 5);
  hTPCPr.Fit("alephBB");
  outFile->cd();
  hTPCPr.Write();
  alephBB.Write();
  outFile->Close();
  inFile->Close();
}