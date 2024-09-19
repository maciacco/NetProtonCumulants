#include "utils.h"

void remove_outlier(TH1D* h, double reject_level = 3.){
  double mean = h->GetMean();
  double std_dev = h->GetStdDev();
  for (int iB{1}; iB <= h->GetNbinsX(); ++iB){
    if (std::abs(h->GetBinCenter(iB) - mean) > reject_level * std_dev){
      h->SetBinContent(iB, 0.);
    }
  }
}

double cbwc(const double *array, const int centbin, const TH1D* hCent){
  double lowEdge = kCentBins[centbin];
  double upEdge = kCentBins[centbin + 1];
  int iL = 0; int iU = kNCentBinsSmall - 1;
  while (kCentBinsSmall[iL] + 0.05 < lowEdge) {
    iL++;
  }
  while (kCentBinsSmall[iU] + 0.05 > upEdge) {
    iU--;
  }
  // std::cout << iL << "\t" << iU << "\t" << lowEdge << "\t" << upEdge << std::endl;
  double avg = 0.;
  double cc = 0.;
  for (int iB{iL}; iB <= iU; ++iB){
    double nc = hCent->GetBinContent(iB + 1);
    avg += (nc * array[iB]);
    cc += nc;
  }

  return cc > 0 ? avg / cc : 0.;
}

void Analysis(const char* period = "18")
{
  TFile f(Form("out_sys_%s_finalBinning.root", period), "recreate");
  TH1D *hSys[kNCentBins];
  TGraphErrors gk2k1;
  gk2k1.SetName("gk2k1");
  gk2k1.SetTitle(";Centrality (%);#k2k1_{#Delta#XiDeltaK}");
  TGraphErrors gk2k1Sys;
  gk2k1Sys.SetName("gk2k1Sys");
  TCanvas cSys("cSys", "cSys");
  cSys.Divide(3, 3);
  for (int i{0}; i < kNCentBins; ++i){
    hSys[i] = new TH1D(Form("hSys_%d", i), ";#kappa_{2}/#kappa_{1};Entries", 500, .0, 2.);
  }

  for(int iVar = 0; iVar < 3; ++iVar)
  {
    std::cout << "var = " << iVar << "..." << std::endl;

    double iB0 = 0.;
    double iB1 = 0.;
    double nSkip = 0;
    double c11_pn_pp_small[N_SAMPLE][100], c1pr_p_small[N_SAMPLE][100], c1pr_n_small[N_SAMPLE][100], c2pr_p_small[N_SAMPLE][100], c2pr_n_small[N_SAMPLE][100], c2pr_pn_small[N_SAMPLE][100];
    double c11_pn_pp[N_SAMPLE][10], c1pr_p[N_SAMPLE][10], c1pr_n[N_SAMPLE][10], c2pr_p[N_SAMPLE][10], c2pr_n[N_SAMPLE][10], c2pr_pn[N_SAMPLE][10];
    for(int sample = 0; sample < N_SAMPLE; sample++)
    {
      TFile *fin = new TFile(Form("output_sys_%d_%d.root", sample, iVar));
      TFile *fCent = TFile::Open(Form("LHC18%d_var_%d.root", sample, iVar));

      TH1D *hCent = (TH1D*)fCent->Get("hCent");

      if (!fin) {nSkip++; fin->Close(); delete fin; continue;}

      TProfile *q1pp_pn = (TProfile*)fin->Get(Form("var_%d/q1_pr_pr_pn", iVar));
      TProfile *q1pr_p = (TProfile*)fin->Get(Form("var_%d/q1_pr_p", iVar));
      TProfile *q2pr_p = (TProfile*)fin->Get(Form("var_%d/q2_pr_p", iVar));
      TProfile *q1pr_n = (TProfile*)fin->Get(Form("var_%d/q1_pr_n", iVar));
      TProfile *q2pr_n = (TProfile*)fin->Get(Form("var_%d/q2_pr_n", iVar));
      TProfile *q1squarepr_p = (TProfile*)fin->Get(Form("var_%d/q1square_pr_p", iVar));
      TProfile *q1squarepr_n = (TProfile*)fin->Get(Form("var_%d/q1square_pr_n", iVar));

      // TH1D *q1pp_pn = (TH1D*)fCent->Get(Form("subsample__%d/hAProtonQ11_Gen_%d", sample + 1, iVar));
      // TH1D *q1pr_p = (TH1D*)fCent->Get(Form("subsample__%d/hMProtonQ1_Gen_%d", sample + 1, iVar));
      // TH1D *q2pr_p = (TH1D*)fCent->Get(Form("subsample__%d/hMProtonQ1_Gen_%d", sample + 1, iVar));
      // TH1D *q1pr_n = (TH1D*)fCent->Get(Form("subsample__%d/hAProtonQ1_Gen_%d", sample + 1, iVar));
      // TH1D *q2pr_n = (TH1D*)fCent->Get(Form("subsample__%d/hAProtonQ1_Gen_%d", sample + 1, iVar));
      // TH1D *q1squarepr_p = (TH1D*)fCent->Get(Form("subsample__%d/hMProtonQ1Sq_Gen_%d", sample + 1, iVar));
      // TH1D *q1squarepr_n = (TH1D*)fCent->Get(Form("subsample__%d/hAProtonQ1Sq_Gen_%d", sample + 1, iVar));

      if (!q1pp_pn || !q1pr_n || !q1pr_p || !q2pr_n || !q2pr_p || !q1squarepr_n || !q1squarepr_p){
        std::cout << "skip..." << std::endl;
        nSkip++; fin->Close(); delete fin; continue;
      }

      for(int i = 1; i <= kNCentBinsSmall; i++)
      {
        c11_pn_pp_small[sample][i - 1] = q1pp_pn->GetBinContent(i) - q1pr_p->GetBinContent(i) * q1pr_n->GetBinContent(i);
        c1pr_p_small[sample][i - 1] = q1pr_p->GetBinContent(i);
        c1pr_n_small[sample][i - 1] = q1pr_n->GetBinContent(i);
        c2pr_p_small[sample][i - 1] = q1squarepr_p->GetBinContent(i) - TMath::Power(q1pr_p->GetBinContent(i), 2.0) + q1pr_p->GetBinContent(i) - q2pr_p->GetBinContent(i);
        c2pr_n_small[sample][i - 1] = q1squarepr_n->GetBinContent(i) - TMath::Power(q1pr_n->GetBinContent(i), 2.0) + q1pr_n->GetBinContent(i) - q2pr_n->GetBinContent(i);
        c2pr_pn_small[sample][i - 1] = c2pr_p_small[sample][i - 1] + c2pr_n_small[sample][i - 1] - 2 * c11_pn_pp_small[sample][i - 1];
        // std::cout << c11_pn_pp_small[sample][i - 1] << std::endl;
      }

      for(int i = 1; i <= kNCentBins; i++)
      {
        c11_pn_pp[sample][i - 1] = cbwc(c11_pn_pp_small[sample], i - 1, hCent);
        c1pr_p[sample][i - 1] = cbwc(c1pr_p_small[sample], i - 1, hCent);
        c1pr_n[sample][i - 1] = cbwc(c1pr_n_small[sample], i - 1, hCent);
        c2pr_p[sample][i - 1] = cbwc(c2pr_p_small[sample], i - 1, hCent);
        c2pr_n[sample][i - 1] = cbwc(c2pr_n_small[sample], i - 1, hCent);
        c2pr_pn[sample][i - 1] = cbwc(c2pr_pn_small[sample], i - 1, hCent);
      }

      fin->Close();
      fCent->Close();
      delete fin;
      delete fCent;
    }

    TGraphErrors g;
    g.SetName(Form("g_%d", iVar));
    for(int i = 1; i < kNCentBins; i++)
    {
      double mean = 0.0;
      double rms = 0.0;
      for(int sample = 0; sample < N_SAMPLE; sample++)
      {
        if (c1pr_n[sample][i - 1] + c1pr_p[sample][i - 1] > 0) {
          mean = mean + ( ( c2pr_pn[sample][i - 1] ) / ( c1pr_n[sample][i - 1] + c1pr_p[sample][i - 1] ));
        }
      }
      mean = mean / ( N_SAMPLE - nSkip);
      for(int sample = 0; sample < N_SAMPLE; sample++)
      {
        if (c1pr_n[sample][i - 1] + c1pr_p[sample][i - 1] > 0) {
          rms = rms + TMath::Power(mean - ( c2pr_pn[sample][i - 1] ) / ( c1pr_n[sample][i - 1] + c1pr_p[sample][i - 1] ), 2.0);
        }
      }

      g.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), mean);
      hSys[i - 1]->Fill(mean);
      g.SetPointError(i - 1, 0, TMath::Sqrt(rms / (( N_SAMPLE - nSkip) * (( N_SAMPLE - nSkip) - 1))));

      if (iVar == 0){
        gk2k1.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), mean);
        gk2k1Sys.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), mean);
        gk2k1.SetPointError(i - 1, 0, TMath::Sqrt(rms / (( N_SAMPLE - nSkip) * (( N_SAMPLE - nSkip) - 1))));
      }
    }

    TCanvas c(Form("c_%d", iVar), Form("c_%d", iVar));
    c.cd();
    f.cd();
    g.Write();
  }

  for (int i{0}; i < kNCentBins - 1; ++i){
    // remove_outlier(hSys[i]);
    hSys[i]->Write();
    //hSys[i]->SetStats(0);
    hSys[i]->SetFillStyle(3004);
    hSys[i]->SetLineColor(kBlue);
    hSys[i]->SetFillColor(kBlue);
    hSys[i]->GetXaxis()->SetRangeUser(hSys[i]->GetMean() - 5.* hSys[i]->GetStdDev(), hSys[i]->GetMean() + 5.* hSys[i]->GetStdDev());
    // gk2k1Sys.SetPointError(i, 2, hSys[i]->GetStdDev());
    hSys[i]->SetTitle(Form("mult. class %d", i));
    cSys.cd(i + 1);
    //hSys[i]->Fit("gaus", "LM+");
    hSys[i]->Draw("histo");
  }
  gk2k1.Write();
  gk2k1Sys.Write();
  TCanvas c("c", "c");
  c.cd();
  gk2k1Sys.SetLineWidth(2);
  gk2k1.SetLineWidth(2);
  gk2k1.SetMarkerStyle(20);
  gk2k1.SetMarkerSize(1.1);
  gk2k1Sys.SetLineColor(kRed);
  gk2k1.SetMarkerColor(kRed);
  gk2k1.SetLineColor(kRed);
  gk2k1Sys.Draw("ape5");
  gk2k1Sys.GetXaxis()->SetTitle("Centrality (%)");
  gk2k1Sys.GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}");
  gk2k1.Draw("pesame");
  c.Print("k2k1_18.pdf");
  c.Write();
  cSys.Write();
  f.Close();
}