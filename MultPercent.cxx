#include <TFile.h>
#include <TH2D.h>
#include "utils.h"

constexpr double kOffset = 0.005;

void MultPercent(const char* fname = "mulPlotsOut.root", const char* ofname = "multPercentOut.root"){
  TFile *_file0 = TFile::Open(fname);
  TFile *ofile = TFile::Open(ofname, "recreate");
  auto hTracksVsMult = static_cast<TH2D*>(_file0->Get("hTracksVsMult"));
  for (int iB{1}; iB < kNCentBins + 1; ++iB){
    int iB_low = hTracksVsMult->GetXaxis()->FindBin(kCentBins[iB - 1] + kOffset);
    int iB_up = hTracksVsMult->GetXaxis()->FindBin(kCentBins[iB] - kOffset);
    TH1D* proj = static_cast<TH1D*>(hTracksVsMult->ProjectionY(Form("proj_%d", iB - 1), iB_low, iB_up));
    TH1D* perc = new TH1D(Form("perc_%d", iB - 1), ";#it{N_{Tracks}};", proj->GetNbinsX(), 0, proj->GetNbinsX());
    double thr = 0.10;
    for (int iBB{1}; iBB < proj->GetNbinsX() + 1; ++iBB) {
      double percent = proj->Integral(1, iBB) / proj->Integral(1, proj->GetNbinsX());
      perc->SetBinContent(iBB, percent);
      if (percent > thr && thr < 1) {
        std::cout << iBB /* << "\t" << static_cast<int>(percent / 0.1) * 0.1 */ << ", ";
        thr = 0.10 + (static_cast<int>(percent / 0.10) * 0.10);
      }
    }
    std::cout << "\n";
    ofile->cd();
    perc->Write();
  }
  ofile->Close();
}
