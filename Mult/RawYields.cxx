#include "Func.h"

double LevyTsallis_Func(const double *x, const double *p) {
  /* dN/dpt */

  double pt = x[0];
  double mass = p[0];
  double mt = TMath::Sqrt(pt * pt + mass * mass);
  double n = p[1];
  double C = p[2];
  double norm = p[3];

  double part1 = (n - 1.) * (n - 2.);
  double part2 = n * C * (n * C + mass * (n - 2.));
  double part3 = part1 / part2;
  double part4 = 1. + (mt - mass) / n / C;
  double part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 * LevyTsallis(const Char_t *name, double mass, double n = 5., double C = 0.1, double norm = 1.) {

  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}

void RawYields(const char* inFile = "AnalysisResults_Mult_data", const char* outFile = "rawyields"){
  TFile *_file0 = TFile::Open(Form("../data/%s.root", inFile));
  TFile *_fileEff = TFile::Open("eff.root");
  TFile *_filePrim = TFile::Open("primaryfrac.root");
  TFile *fout = TFile::Open(Form("%s.root", outFile), "recreate");
  TH3D* hRec = (TH3D*)_file0->Get("ebye-mult/RecTracks");
  TH2D* hEv = (TH2D*)_file0->Get("ebye-mult/QA/nTrklCorrelation");

  int multBins[][2] = {{1, 4}, {5, 7}, {8, 10}, {11, 13}, {14, 16}, {17, 20}, {20, 27}, {28, 36}, {37, 100}};
  double multBinsAxis[] = {0, 4, 7, 10, 13, 16, 20, 27, 36, 100};
  TH1D hMultTrkl("hMultTrkl", ";SPD tracklets 0.7 < |#eta| < 1.2;#LT d#it{N}/d#eta #GT_{INEL>0}", 9, multBinsAxis);
  for (int iM{0}; iM < 9; ++iM) {
    TH1D* hEvProj = (TH1D*)hEv->ProjectionX("projEv", multBins[iM][0], multBins[iM][1]);
    auto [hA, hM] = projectAM(hRec, multBins[iM][0], multBins[iM][1], hRec->GetZaxis()->FindBin(-0.1), hRec->GetZaxis()->FindBin(0.1));
    TH1D* hEffA = (TH1D*)_fileEff->Get(Form("hAEff_%d", iM));
    TH1D* hEffM = (TH1D*)_fileEff->Get(Form("hMEff_%d", iM));
    TH1D* hFprimA = (TH1D*)_filePrim->Get(Form("hAPrimFrac_%d", iM));
    TH1D* hFprimM = (TH1D*)_filePrim->Get(Form("hMPrimFrac_%d", iM));
    fout->cd();
    hA->Divide(hEffA);
    hM->Divide(hEffM);
    hA->Multiply(hFprimA);
    hM->Multiply(hFprimM);
    std::cout << hEvProj->Integral(1, hEvProj->GetNbinsX()) << std::endl;
    hA->Scale(1. / 1.2 / hEvProj->Integral(1, hEvProj->GetNbinsX()));
    hM->Scale(1. / 1.2 / hEvProj->Integral(1, hEvProj->GetNbinsX()));
    double errA = 0;
    double errM = 0;
    double integA = hA->IntegralAndError(1, hA->GetNbinsX(), errA);
    double integM = hM->IntegralAndError(1, hM->GetNbinsX(), errM);
    TF1* fitFunA = LevyTsallis("fitFuncA", 0.139);
    hA->Fit("fitFuncA", "RM+", "", 0.05, .8);

    TF1* fitFunM = LevyTsallis("fitFuncM", 0.139);
    hM->Fit("fitFuncM", "RM+", "", 0.05, .8);

    double extrapA = fitFunA->Integral(0., 0.05);
    double extrapM = fitFunM->Integral(0., 0.05);

    integA += extrapA; // TODO: propagate uncertainty
    integM += extrapM; // TODO: propagate uncertainty

    hA->Write();
    hM->Write();


    std::cout << integA + integM << " +/- " << std::hypot(errA, errM) << std::endl;
    hMultTrkl.SetBinContent(iM + 1, integA + integM);
    hMultTrkl.SetBinError(iM + 1, std::hypot(errA, errM));
  }

  hMultTrkl.Write();
  fout->Close();
  _file0->Close();
}