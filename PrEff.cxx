#include "utils.h"
#include "Utils.h"

//using namespace utils;

void PrEff(const char* inFileName = "LHC21pp0_var", const char* outFileName = "prEff", const int iVarMin = 364, const int iVarMax = 365){
  gStyle->SetOptStat(0);
  TFile *fOut = TFile::Open(Form("%s/%s.root", kResDir, outFileName), "recreate");
  for (int iVar{iVarMin}; iVar < iVarMax; ++iVar){
    TFile *fMC = TFile::Open(Form("%s/%s_%d.root", kResDir, inFileName, iVar));
    for (int iS = 0; iS < 1; ++iS){
      fOut->mkdir(Form("subsample_%d_var_%d", iS + 1, iVar));
      fOut->cd(Form("subsample_%d_var_%d", iS + 1, iVar));
      for (int iP = 0; iP < 1; ++iP){
        for (int iM = 0; iM < 2; ++iM){
          // TTList *l = (TTList*)fMC->Get("nuclei_kaon_mcTrue_");
          // auto hGenKaon = (TH2D*)l->Get(Form("f%sTotal", kAntiMatterLabel[iM]));
          // auto hRecKaon = (TH2D*)l->Get(Form("f%sITS_TPC", kAntiMatterLabel[iM]));
          auto hGen = (TH3D*)fMC->Get(Form(/* subsample__%d/ */"h%sGen%s_%d",/*  iS + 1,  */kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          auto hRec = (TH3D*)fMC->Get(Form(/* subsample__%d/ */"h%sRec%s_%d",/*  iS + 1,  */kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          auto hGenTrkl = (TH3D*)fMC->Get(Form(/* subsample__%d/ */"h%sGen%sTrkl_%d",/*  iS + 1,  */kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          auto hRecTrkl = (TH3D*)fMC->Get(Form(/* subsample__%d/ */"h%sRec%sTrkl_%d",/*  iS + 1,  */kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          TCanvas ce(Form("c%sEff_%s_%d", kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar), Form("c%sEff_%s", kAntiMatterLabel[iM], kPartLabelExtend[iP]), 500, 500);
          TLegend leg(0.2, 0.4, 0.3, 0.8);
          leg.SetTextFont(44);
          leg.SetTextSize(20);
          TH1D *hEff[kNCentBins][kNEtaBins];
          TH1D *hEffMult[kNCentBins][kNEtaBins];
          TH1D *hEffTrkl[kNTrklBins][kNEtaBins];
          TH1D *hEffMultTrkl[kNTrklBins][kNEtaBins];
          for (int iC = 0; iC < kNCentBins; ++iC){
            for (int iE = 0; iE < kNEtaBins; ++iE){
              auto hGenProjMult = (TH1D*)hGen->ProjectionY(Form("h%sGen%sProjMult_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), iC + 1, iC + 1, iE + 1, iE + 1);
              auto hRecProjMult = (TH1D*)hRec->ProjectionY(Form("h%sRec%sProjMult_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), iC + 1, iC + 1, iE + 1, iE + 1);
              auto hGenProj = (TH1D*)hGen->ProjectionY(Form("h%sGen%sProj_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), 1, kNCentBins, iE + 1, iE + 1);
              auto hRecProj = (TH1D*)hRec->ProjectionY(Form("h%sRec%sProj_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), 1, kNCentBins, iE + 1, iE + 1);
              hEffMult[iC][iE] = new TH1D(*hGenProjMult);
              hEff[iC][iE] = new TH1D(*hGenProj);
              hEff[iC][iE]->SetName(Form("h%sEff%s_%d_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE, iVar));
              hEff[iC][iE]->SetTitle(";#it{p}_{T} (GeV/#it{c});Efficiency");
              hEff[iC][iE]->Divide(hRecProj, hGenProj, 1, 1, "B");
              hEffMult[iC][iE]->Divide(hRecProjMult, hGenProjMult, 1, 1, "B");
              hEff[iC][iE]->SetMarkerColor(colors[iC]);
              hEff[iC][iE]->SetLineColor(colors[iC]);
              hEff[iC][iE]->SetLineWidth(2);
              fOut->cd();
              fOut->cd(Form("subsample_%d_var_%d", iS + 1, iVar));
              hEff[iC][iE]->Write();
              ce.cd();
              //hEff[iC][iE]->GetYaxis()->SetRangeUser(0., iP == 1 ? 0.15 : 1.);
              //hEff[iC][iE]->GetXaxis()->SetRangeUser(iP == 1 ? 1. : .2, iP == 1 ? 3. : 1.2);
              leg.AddEntry(hEff[iC][iE], Form("%.0f-%.0f%%", kCentBins[iC], kCentBins[iC + 1]));
              hEffMult[iC][iE]->SetMarkerColor(colors[iC]);
              hEffMult[iC][iE]->SetLineColor(colors[iC]);
              hEffMult[iC][iE]->SetLineWidth(2);
              // hEff[iC][iE]->Draw(iC == 0 && iE == 0 ? "pe" : "pesame");
              hEffMult[iC][iE]->Divide(hEff[iC][iE]);
              hEffMult[iC][iE]->Write();
              hEffMult[iC][iE]->Draw(iC == 0 && iE == 0 ? "pe" : "pesame");
              //hGenProj->Write();
              //hRecProj->Write();
            }
          }
          for (int iTrkl = 0; iTrkl < kNTrklBins; ++iTrkl){
            for (int iE = 0; iE < kNEtaBins; ++iE){
              auto hGenProjMult = (TH1D*)hGenTrkl->ProjectionY(Form("h%sGen%sProjMultTrkl_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iTrkl, iE), iTrkl + 1, iTrkl + 1, iE + 1, iE + 1);
              auto hRecProjMult = (TH1D*)hRecTrkl->ProjectionY(Form("h%sRec%sProjMultTrkl_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iTrkl, iE), iTrkl + 1, iTrkl + 1, iE + 1, iE + 1);
              auto hGenProj = (TH1D*)hGenTrkl->ProjectionY(Form("h%sGen%sProjTrkl_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iTrkl, iE), 1, kNTrklBins, iE + 1, iE + 1);
              auto hRecProj = (TH1D*)hRecTrkl->ProjectionY(Form("h%sRec%sProjTrkl_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iTrkl, iE), 1, kNTrklBins, iE + 1, iE + 1);
              hEffMultTrkl[iTrkl][iE] = new TH1D(*hGenProjMult);
              hEffTrkl[iTrkl][iE] = new TH1D(*hGenProj);
              hEffTrkl[iTrkl][iE]->SetName(Form("h%sEffTrkl%s_%d_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iTrkl, iE, iVar));
              hEffTrkl[iTrkl][iE]->SetTitle(";#it{p}_{T} (GeV/#it{c});Efficiency");
              hEffTrkl[iTrkl][iE]->Divide(hRecProj, hGenProj, 1, 1, "B");
              hEffMultTrkl[iTrkl][iE]->Divide(hRecProjMult, hGenProjMult, 1, 1, "B");
              // hEffTrkl[iTrkl][iE]->SetMarkerColor(colors[iTrkl]);
              // hEffTrkl[iTrkl][iE]->SetLineColor(colors[iTrkl]);
              hEffTrkl[iTrkl][iE]->SetLineWidth(2);
              fOut->cd();
              fOut->cd(Form("subsample_%d_var_%d", iS + 1, iVar));
              hEffTrkl[iTrkl][iE]->Write();
              ce.cd();
              // leg.AddEntry(hEffTrkl[iTrkl][iE], Form("%.0f-%.0f%%", kCentBins[iTrkl], kCentBins[iTrkl + 1]));
              //hEffMultTrkl[iTrkl][iE]->SetMarkerColor(colors[iTrkl]);
              //hEffMultTrkl[iTrkl][iE]->SetLineColor(colors[iTrkl]);
              hEffMultTrkl[iTrkl][iE]->SetLineWidth(2);
              // hEffTrkl[iTrkl][iE]->Draw(iTrkl == 0 && iE == 0 ? "pe" : "pesame");
              hEffMultTrkl[iTrkl][iE]->Divide(hEffTrkl[iTrkl][iE]);
              hEffMultTrkl[iTrkl][iE]->Write();
              // hEffMultTrkl[iTrkl][iE]->Draw(iTrkl == 0 && iE == 0 ? "pe" : "pesame");
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
