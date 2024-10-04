#include "Func.h"

void Eff(const char* inFile = "AnalysisResults_Mult_mc", const char* outFile = "eff"){
  TFile *_file0 = TFile::Open(Form("../data/%s.root", inFile));
  TFile *fout = TFile::Open(Form("%s.root", outFile), "recreate");
  TH3D* hRec = (TH3D*)_file0->Get("ebye-mult/RecPart");
  TH3D* hGen = (TH3D*)_file0->Get("ebye-mult/GenPart");

  int multBins[][2] = {{1, 4}, {5, 7}, {8, 10}, {11, 13}, {14, 16}, {17, 20}, {20, 27}, {28, 36}, {37, 100}};
  for (int iM{0}; iM < 9; ++iM) {
    auto [hA, hM] = projectAM(hRec, multBins[iM][0], multBins[iM][1], 1, 5);
    auto [hAG, hMG] = projectAM(hGen, multBins[iM][0], multBins[iM][1], 1, 5);
    fout->cd();
    // hA->Write();
    // hM->Write();
    // hAG->Write();
    // hMG->Write();
    hA->Divide(hA, hAG, 1., 1., "B");
    hM->Divide(hM, hMG, 1., 1., "B");
    hA->SetName(Form("hAEff_%d", iM));
    hM->SetName(Form("hMEff_%d", iM));
    hA->Write();
    hM->Write();
  }
  fout->Close();
  _file0->Close();
}