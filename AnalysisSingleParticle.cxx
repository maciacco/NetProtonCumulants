#include "utils.h"

void AnalysisSingleParticle()
{
  TFile f("out_sys_finalBinning_singleParticle.root", "recreate");
  TH1D *hSys[kNCentBins];
  TGraphErrors gObs;
  gObs.SetName("gObs");
  gObs.SetTitle(";Centrality (%);#kappa_{2}/#kappa_{1}");
  TGraphErrors gObsSys;
  gObsSys.SetName("gObsSys");
  TCanvas cSys("cSys", "cSys");
  cSys.Divide(3, 3);
  for (int i{0}; i < kNCentBins; ++i){
    hSys[i] = new TH1D(Form("hSys_%d", i), ";#kappa_{2}/#kappa_{1};Entries", 10, -1, -2);
  }

  for(int iVar = 364; iVar < 365; ++iVar)
  {
    // std::cout << "var = " << iVar << "..." << std::endl;
    double nSkip = 0;
    double k11_pn_prpr[kNSample][10], k1pr_p[kNSample][10], k1pr_n[kNSample][10], k2pr_p[kNSample][10], k2pr_n[kNSample][10], k2pr_pn[kNSample][10];
    for(int sample = 0; sample < kNSample; sample++)
    {
      TFile *fin = new TFile(Form("%s/output_sys_singleParticle_%d_%d.root", kResDir, sample, iVar));

      if (!fin) {nSkip++; fin->Close(); delete fin; continue;}

      TProfile *q1prpr_pn = (TProfile*)fin->Get(Form("var_%d/q1_pr_pr_pn", iVar));
      TProfile *q1kk_pn = (TProfile*)fin->Get(Form("var_%d/q1_kaon_kaon_pn", iVar));
      TProfile *q1pr_p = (TProfile*)fin->Get(Form("var_%d/q1_pr_p", iVar));
      TProfile *q2pr_p = (TProfile*)fin->Get(Form("var_%d/q2_pr_p", iVar));
      TProfile *q1pr_n = (TProfile*)fin->Get(Form("var_%d/q1_pr_n", iVar));
      TProfile *q2pr_n = (TProfile*)fin->Get(Form("var_%d/q2_pr_n", iVar));
      TProfile *q1squarepr_p = (TProfile*)fin->Get(Form("var_%d/q1square_pr_p", iVar));
      TProfile *q1squarepr_n = (TProfile*)fin->Get(Form("var_%d/q1square_pr_n", iVar));

      if (!q1prpr_pn || !q1pr_n || !q1pr_p || !q2pr_n || !q2pr_p || !q1squarepr_n || !q1squarepr_p){
        nSkip++; fin->Close(); delete fin; continue;
      }

      for(int i = 1; i <= kNCentBins; i++)
      {
        k11_pn_prpr[sample][i - 1] = q1prpr_pn->GetBinContent(i) - q1pr_p->GetBinContent(i) * q1pr_n->GetBinContent(i);
        k1pr_p[sample][i - 1] = q1pr_p->GetBinContent(i);
        k1pr_n[sample][i - 1] = q1pr_n->GetBinContent(i);
        k2pr_p[sample][i - 1] = q1squarepr_p->GetBinContent(i) - TMath::Power(q1pr_p->GetBinContent(i), 2.0) + q1pr_p->GetBinContent(i) - q2pr_p->GetBinContent(i);
        k2pr_n[sample][i - 1] = q1squarepr_n->GetBinContent(i) - TMath::Power(q1pr_n->GetBinContent(i), 2.0) + q1pr_n->GetBinContent(i) - q2pr_n->GetBinContent(i);
        k2pr_pn[sample][i - 1] = k2pr_p[sample][i - 1] + k2pr_n[sample][i - 1] - 2 * k11_pn_prpr[sample][i - 1];
      }

      fin->Close();
      delete fin;
    }

    TGraphErrors g;
    g.SetName(Form("g_%d", iVar));
    for(int i = 1; i < kNCentBins; i++)
    {
      double rhomean = 0.0;
      double rhorms = 0.0;
      for(int sample = 0; sample < kNSample; sample++)
      {
        rhomean = rhomean + ( ( k2pr_n[sample][i - 1] / k1pr_n[sample][i - 1] ));
      }
      rhomean = rhomean / ( kNSample - nSkip);

      for(int sample = 0; sample < kNSample; sample++)
      {
        rhorms = rhorms + TMath::Power(rhomean - ( k2pr_n[sample][i - 1] / k1pr_n[sample][i - 1] ), 2.0);
      }

      g.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), rhomean);
      hSys[i - 1]->Fill(rhomean);

      g.SetPointError(i - 1, 0, TMath::Sqrt(rhorms / (( kNSample - nSkip) * (( kNSample - nSkip) - 1))));
      if (iVar == 364){
        gObs.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), rhomean);
        gObsSys.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), rhomean);
        gObs.SetPointError(i - 1, 0, TMath::Sqrt(rhorms / (( kNSample - nSkip) * (( kNSample - nSkip) - 1))));
      }
    }

    TCanvas c(Form("c_%d", iVar), Form("c_%d", iVar));
    c.cd();
    f.cd();
    g.Write();
  }

  for (int i{0}; i < kNCentBins - 1; ++i){
    hSys[i]->Write();
    hSys[i]->SetFillStyle(3004);
    hSys[i]->SetLineColor(kBlue);
    hSys[i]->SetFillColor(kBlue);
    hSys[i]->GetXaxis()->SetRangeUser(hSys[i]->GetMean() - 5.*hSys[i]->GetStdDev(), hSys[i]->GetMean() + 5.*hSys[i]->GetStdDev());
    gObsSys.SetPointError(i, 2, hSys[i]->GetStdDev());
    hSys[i]->SetTitle(Form("mult. class %d", i));
    cSys.cd(i + 1);
    hSys[i]->Draw("histo");
  }
  gObs.Write();
  gObsSys.Write();
  TCanvas c("c", "c");
  c.cd();
  gObsSys.SetLineWidth(2);
  gObs.SetLineWidth(2);
  gObs.SetMarkerStyle(20);
  gObs.SetMarkerSize(1.1);
  gObsSys.SetLineColor(kRed);
  gObs.SetMarkerColor(kRed);
  gObs.SetLineColor(kRed);
  gObsSys.Draw("ape5");
  gObsSys.GetXaxis()->SetTitle("Centrality (%)");
  gObsSys.GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}");
  gObs.Draw("pesame");
  c.Write();
  cSys.Write();
  f.Close();
}