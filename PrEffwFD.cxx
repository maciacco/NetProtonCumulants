#include "utils.h"
#include "Utils.h"

//using namespace utils;

void PrEffwFD(const char* inFileName = "LHC22f30_var", const char* inFileNameWD = "LHC22f3_WD0_var", const char* outFileName = "prEff", const int iVarMin = 364, const int iVarMax = 365, const bool MB = true){
  gStyle->SetOptStat(0);
  const char* dir = kTriggerSel == 0x1 ? "MB" : "HM";

  TFile *fOut = TFile::Open(Form("%s/%s/%s_%d_%d.root", kResDir, dir, outFileName, iVarMin, iVarMax), "recreate");
  for (int iVar{iVarMin}; iVar < iVarMax; ++iVar){
    bool inVars = false;
    for (int iV{0}; iV < kNVar; ++iV) {
      if (iVar == kVar[iV]) {
        inVars = true;
        break;
      }
    }
    if (!inVars) continue;
    const char* dir = kTriggerSel ? "MB" : "HM";
    TFile *fPt = TFile::Open(MB ? "data/la2prPt.root" : "data/la2prPtHM.root");
    TFile *fMC = TFile::Open(Form("%s/%s/%s_%d.root", kResDir, dir, inFileName, iVar));
    TFile *fMCWD = TFile::Open(Form("%s/%s/%s_%d.root", kResDir, dir, inFileNameWD, iVar));
    for (int iS = 0; iS < 1; ++iS){
      fOut->mkdir(Form("subsample_%d_var_%d", iS + 1, iVar));
      fOut->cd(Form("subsample_%d_var_%d", iS + 1, iVar));
      for (int iP = 0; iP < 1; ++iP){
        for (int iM = 0; iM < 2; ++iM){
          auto hGenP = (TH3D*)fMC->Get(Form("h%sGen%s_%d",kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          auto hRecP = (TH3D*)fMC->Get(Form("h%sRec%s_%d",kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          auto hGenWD = (TH3D*)fMCWD->Get(Form("h%sGen%s_%d",kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          auto hRecWD = (TH3D*)fMCWD->Get(Form("h%sRec%s_%d",kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          auto hGenTrkl = (TH3D*)fMC->Get(Form("h%sGen%sTrkl_%d",kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
          auto hRecTrkl = (TH3D*)fMC->Get(Form("h%sRec%sTrkl_%d",kAntiMatterLabel[iM], kPartLabelExtend[iP], iVar));
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
              auto ptPrWD = (TH1D*)fPt->Get(Form("hPrPt_%d", iC));
              auto ptPrWDSig = (TH1D*)fPt->Get(Form("hPrPtSig_%d", iC));
              // ptPrWD->Add(ptPrWDSig);
              auto ptPrPrim = (TH1D*)fPt->Get(Form("hPrPtPrim_%d", iC));
              auto hGenProjMult = (TH1D*)hGenP->ProjectionY(Form("h%sGen%sProjMult_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), iC + 1, iC + 1, iE + 1, iE + 1);
              auto hRecProjMult = (TH1D*)hRecP->ProjectionY(Form("h%sRec%sProjMult_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), iC + 1, iC + 1, iE + 1, iE + 1);
              auto hGenProjMultWD = (TH1D*)hGenWD->ProjectionY(Form("h%sGenWD%sProjMult_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), iC + 1, iC + 1, iE + 1, iE + 1);
              auto hRecProjMultWD = (TH1D*)hRecWD->ProjectionY(Form("h%sRecWD%sProjMult_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), iC + 1, iC + 1, iE + 1, iE + 1);
              auto hGenProj = (TH1D*)hGenP->ProjectionY(Form("h%sGen%sProj_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), 1, kNCentBins, iE + 1, iE + 1);
              auto hRecProj = (TH1D*)hRecP->ProjectionY(Form("h%sRec%sProj_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), 1, kNCentBins, iE + 1, iE + 1);
              auto hGenProjWD = (TH1D*)hGenWD->ProjectionY(Form("h%sGenWD%sProj_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), 1, kNCentBins, iE + 1, iE + 1);
              auto hRecProjWD = (TH1D*)hRecWD->ProjectionY(Form("h%sRecWD%sProj_%d_%d", kAntiMatterLabel[iM], kPartLabel[iP], iC, iE), 1, kNCentBins, iE + 1, iE + 1);

              hRecProj->Divide(hRecProj, hGenProj, 1., 1., "B");
              hRecProjWD->Divide(hRecProjWD, hGenProjWD, 1., 1., "B");
              hRecProj->Multiply(ptPrPrim);
              hRecProjWD->Multiply(ptPrWD);
              hRecProj->Add(hRecProjWD);

              hGenProj->Add(ptPrPrim, ptPrWD);

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
