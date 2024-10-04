#include "Func.h"

constexpr double DCAxyCut = 0.1;

void PrimaryFraction(const char* inFile = "AnalysisResults_Mult_data", const char* inFileMC = "AnalysisResults_Mult_mc", const char* outFile = "primaryfrac"){
  TFile *_file0 = TFile::Open(Form("../data/%s.root", inFile));
  TFile *_file0_mc = TFile::Open(Form("../data/%s.root", inFileMC));
  TFile *_fileEff = TFile::Open("eff.root");
  TFile *fout = TFile::Open(Form("%s.root", outFile), "recreate");
  TH3D* hRec = (TH3D*)_file0->Get("ebye-mult/RecTracks");
  TH3D* hPrim = (TH3D*)_file0_mc->Get("ebye-mult/PrimTracks");
  TH3D* hSecWD = (TH3D*)_file0_mc->Get("ebye-mult/SecWDTracks");
  TH3D* hSec = (TH3D*)_file0_mc->Get("ebye-mult/SecTracks");

  int multBins[][2] = {{1, 4}, {5, 7}, {8, 10}, {11, 13}, {14, 16}, {17, 20}, {20, 27}, {28, 36}, {37, 100}};
  const char* AM[] = {"A", "M"};
  for (int iMult{0}; iMult < 9; ++iMult) {
    TH1D* hPrimFrac[2];
    int iM = 0;
    auto [hA, hM] = projectAM2D(hRec, multBins[iMult][0], multBins[iMult][1]);
    auto [hAPrim, hMPrim] = projectAM2D(hPrim, multBins[iMult][0], multBins[iMult][1]);
    auto [hASecWD, hMSecWD] = projectAM2D(hSecWD, multBins[iMult][0], multBins[iMult][1]);
    auto [hASec, hMSec] = projectAM2D(hSec, multBins[iMult][0], multBins[iMult][1]);
    for (auto hM : {hA, hM}) {
      hPrimFrac[iM] = new TH1D(Form("h%sPrimFrac_%d", AM[iM], iMult), ";#it{p}_{T} (GeV/#it{c});#it{f}_{prim}", hM->GetNbinsX(), hM->GetXaxis()->GetBinLowEdge(1), hM->GetXaxis()->GetBinUpEdge(hM->GetNbinsX()));
      for (int iB{1}; iB < hM->GetNbinsX(); ++iB) {
        if (hM->GetXaxis()->GetBinLowEdge(iB) > 2.425 || hM->GetXaxis()->GetBinLowEdge(iB) < .04) continue;

        TH1D* hMProj = (TH1D*)hM->ProjectionY("hMProj", iB, iB);
        TH1D* hMPrimProj = (TH1D*)hMPrim->ProjectionY("hMPrimProj", iB, iB);
        TH1D* hMSecWDProj = (TH1D*)hMSecWD->ProjectionY("hMSecWDProj", iB, iB);
        TH1D* hMSecProj = (TH1D*)hMSec->ProjectionY("hMSecProj", iB, iB);
        // hMProj->Rebin(2);
        // hMPrimProj->Rebin(2);
        // hMSecWDProj->Rebin(2);
        // hMSecProj->Rebin(2);
        // const int kNDCABins = 34;
        // const double kDCABins[kNDCABins + 1] = {-1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.07, -0.05, -0.04, -0.02, 0.00, 0.02, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00};

        // hMProj = (TH1D*)hMProj->Rebin(kNDCABins, "hMProj", kDCABins);
        // hMPrimProj = (TH1D*)hMPrimProj->Rebin(kNDCABins, "hMPrimProj", kDCABins);
        // hMSecWDProj = (TH1D*)hMSecWDProj->Rebin(kNDCABins, "hMSecWDProj", kDCABins);
        // hMSecProj = (TH1D*)hMSecProj->Rebin(kNDCABins, "hMSecProj", kDCABins);

        TObjArray *mc = new TObjArray(2); // MC histograms are put in this array
        mc->Add(hMPrimProj);
        mc->Add(hMSecWDProj);
        // mc->Add(hMSecProj);
        TFractionFitter *fit = new TFractionFitter(hMProj, mc, "Q"); // initialise
        ROOT::Fit::Fitter *fitter = fit->GetFitter();
        int status = -100;

        fit->Constrain(0, 0., 1.);
        fit->Constrain(1, 0., 1.);
        double fitRange = 1.;
        fit->SetRangeX(hMProj->FindBin(-fitRange), hMProj->FindBin(fitRange - 0.001));

        for(int i{0}; i < 2; ++i) status = fit->Fit();
        TH1F *result = (TH1F *)fit->GetPlot();
        TH1F *mc1 = (TH1F *)fit->GetMCPrediction(0);
        TH1F *mc2 = (TH1F *)fit->GetMCPrediction(1);
        mc1->SetName(Form("%s_Prediction", hMPrimProj->GetName()));
        mc2->SetName(Form("%s_Prediction", hMSecWDProj->GetName()));
        double integralData = hMProj->Integral();
        double fracMc1, fracMc2, fracMc3, errFracMc1, errFracMc2, errFracMc3;
        fit->GetResult(0, fracMc1, errFracMc1);
        fit->GetResult(1, fracMc2, errFracMc2);
        mc1->Scale(fracMc1 * integralData / mc1->Integral(), "width");
        mc2->Scale(fracMc2 * integralData / mc2->Integral(), "width");
        hMProj->Scale(1, "width");
        result->Scale(1, "width");

        double intPrimDCAcutError = 0.;
        double intPrimDCAcut = mc1->IntegralAndError(result->FindBin(-DCAxyCut), result->FindBin(DCAxyCut - 0.001), intPrimDCAcutError);
        double intResDCAcutError = 0.;
        double intResDCAcut = result->IntegralAndError(result->FindBin(-DCAxyCut), result->FindBin(DCAxyCut - 0.001), intResDCAcutError);
        double primaryRatio = intPrimDCAcut / intResDCAcut;
        double primaryRatioError = errFracMc1 / fracMc1 * primaryRatio;

        // std::cout << primaryRatio << std::endl;
        hPrimFrac[iM]->SetBinContent(iB, primaryRatio);
        hPrimFrac[iM]->SetBinError(iB, primaryRatioError);

        fout->cd();
        // hA->Write();
        // hM->Write();
        result->Write();
        // hAPrim->Write();
        // hMPrim->Write();
        // hASecWD->Write();
        // hMSecWD->Write();
        // hASec->Write();
        // hMSec->Write();
        hMProj->Write();
        mc1->Write();
        mc2->Write();
      }
      hPrimFrac[iM]->Write();
      ++iM;
    }
  }
  // (TH1F *)fit->GetMCPrediction(2)->Write();
  fout->Close();
  _file0->Close();
}