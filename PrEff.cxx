#include "utils.h"
#include "Utils.h"

using namespace utils;

void PrEff(const char* inFileName = "LHC18_var", const char* outFileName = "prEff"){
  gStyle->SetOptStat(0);
  TFile *fOut = TFile::Open(Form("%s/%s.root", kResDir, outFileName), "recreate");
  for (int iVar{0}; iVar < 3; ++iVar){
    TFile *fMC = TFile::Open(Form("%s/%s_%d.root", kResDir, inFileName, iVar));
    for (int iS = 0; iS < 1; ++iS){
      fOut->mkdir(Form("subsample_%d_var_%d", iS + 1, iVar));
      fOut->cd(Form("subsample_%d_var_%d", iS + 1, iVar));
      for (int iP = 0; iP < 1; ++iP){
        for (int iM = 0; iM < 2; ++iM){
          // TTList *l = (TTList*)fMC->Get("nuclei_kaon_mcTrue_");
          // auto hGenKaon = (TH2D*)l->Get(Form("f%sTotal", kAntiMatterLabel[iM]));
          // auto hRecKaon = (TH2D*)l->Get(Form("f%sITS_TPC", kAntiMatterLabel[iM]));
          auto hGen = (TH2D*)fMC->Get(Form(/* subsample__%d/ */"h%sGen%s_%d",/*  iS + 1,  */kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          auto hRec = (TH2D*)fMC->Get(Form(/* subsample__%d/ */"h%sRec%s_%d",/*  iS + 1,  */kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          TCanvas ce(Form("c%sEff_%s_%d", kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar), Form("c%sEff_%s", kAntiMatterLabel[iM], kPartLabelExtend[iP]), 500, 500);
          TLegend leg(0.2, 0.4, 0.3, 0.8);
          leg.SetTextFont(44);
          leg.SetTextSize(20);
          TH1D *hEff[kNCentBins][kNEtaBins];
          for (int iC = 0; iC < kNCentBins; ++iC){
            for (int iE = 0; iE < kNEtaBins; ++iE){
              auto hGenProj = (TH1D*)hGen->ProjectionY(Form("h%sGen%sProj_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), iC + 1, iC + 1/* , iE + 1, iE + 1 */);
              auto hRecProj = (TH1D*)hRec->ProjectionY(Form("h%sRec%sProj_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), iC + 1, iC + 1/* , iE + 1, iE + 1 */);
              hEff[iC][iE] = new TH1D(*hGenProj);
              hEff[iC][iE]->SetName(Form("h%sEff%s_%d_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE, iVar));
              hEff[iC][iE]->SetTitle(";#it{p}_{T} (GeV/#it{c});Efficiency");
              hEff[iC][iE]->Divide(hRecProj, hGenProj, 1, 1, "B");
              hEff[iC][iE]->SetMarkerColor(colors[iC]);
              hEff[iC][iE]->SetLineColor(colors[iC]);
              hEff[iC][iE]->SetLineWidth(2);
              fOut->cd();
              fOut->cd(Form("subsample_%d_var_%d", iS + 1, iVar));
              hEff[iC][iE]->Write();
              ce.cd();
              hEff[iC][iE]->GetYaxis()->SetRangeUser(0., iP == 1 ? 0.15 : 1.);
              hEff[iC][iE]->GetXaxis()->SetRangeUser(iP == 1 ? 1. : .2, iP == 1 ? 3. : 1.);
              leg.AddEntry(hEff[iC][iE], Form("%.0f-%.0f%%", kCentBins[iC], kCentBins[iC + 1]));
              hEff[iC][iE]->Draw(iC == 0 && iE == 0 ? "pe" : "pesame");
              //hGenProj->Write();
              //hRecProj->Write();
            }
          }
          ce.cd();
          leg.Draw("same");
          ce.Write();
          //ce.Print(Form("c%sEff_%s_18qr.pdf", kAntiMatterLabel[iM], kPartLabel[iP]));
        }
      }
    }
    fMC->Close();
  }
  fOut->Close();
}
