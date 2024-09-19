#include "utils.h"

#define FILL_MC

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

void k2_to_k2sk(double &mean, double &rms, const double *k2sk, const double *k2, const double nSkip = 0.){
  mean = 0.0;
  rms = 0.0;
  for(int sample = 0; sample < N_SAMPLE; sample++)
  {
    if (k2sk[sample] > 0) {
      mean = mean + (k2[sample] / k2sk[sample]);
    }
  }
  mean = mean / ( N_SAMPLE - nSkip);
  for(int sample = 0; sample < N_SAMPLE; sample++)
  {
    if (k2sk[sample] > 0) {
      rms = rms + powI(mean - (k2[sample] / k2sk[sample]), 2);
    }
  }
}

void Analysis(const char* period = "18")
{
  TFile f(Form("out_sys_%s_finalBinning.root", period), "recreate");
  TH1D *hSys[kNCentBins];
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

    #ifdef FILL_MC
      double k2sk_small_gen[N_SAMPLE][100];
      double k2_small_gen[N_SAMPLE][100];
      double k2sk_gen[10][N_SAMPLE];
      double k2_gen[10][N_SAMPLE];
    #endif

    double Q1_small;
    double Q2_small;
    double k2sk_small[N_SAMPLE][100];
    double k2_small[N_SAMPLE][100];
    double k2sk[10][N_SAMPLE];
    double k2[10][N_SAMPLE];

    for(int sample = 0; sample < N_SAMPLE; sample++)
    {
      TFile *fin = new TFile(Form("%s/output_sys_%d_%d.root", kResDir, sample, iVar));
      TFile *fCent = TFile::Open(Form("%s/LHC18%d_var_%d.root", kResDir, sample, iVar));

      TH1D *hCent = (TH1D*)fCent->Get("hCent");

      if (!fin) {nSkip++; fin->Close(); delete fin; continue;}

      TProfile *q1_1_1 = (TProfile*)fin->Get(Form("var_%d/q1_1_1", iVar));
      TProfile *q2_1_1 = (TProfile*)fin->Get(Form("var_%d/q2_1_1", iVar));
      TProfile *q1_2_1 = (TProfile*)fin->Get(Form("var_%d/q1_2_1", iVar));
      TProfile *q1_2_2 = (TProfile*)fin->Get(Form("var_%d/q1_2_2", iVar));

      #ifdef FILL_MC
        TProfile *N1p = (TProfile*)fin->Get(Form("var_%d/N1p", iVar));
        TProfile *N1 = (TProfile*)fin->Get(Form("var_%d/N1", iVar));
        TProfile *N2 = (TProfile*)fin->Get(Form("var_%d/N2", iVar));
        TProfile *N3 = (TProfile*)fin->Get(Form("var_%d/N3", iVar));
        TProfile *N4 = (TProfile*)fin->Get(Form("var_%d/N4", iVar));
        TProfile *N5 = (TProfile*)fin->Get(Form("var_%d/N5", iVar));
        TProfile *N6 = (TProfile*)fin->Get(Form("var_%d/N6", iVar));
      #endif // FILL_MC

      if (!q1_1_1 || !q2_1_1 || !q1_2_1 || !q1_2_2){
        std::cout << "skip..." << std::endl;
        nSkip++; fin->Close(); delete fin; continue;
      }

      for(int i = 1; i <= kNCentBinsSmall; i++)
      {
        Q1_small = q1_1_1->GetBinContent(i);
        Q2_small = q2_1_1->GetBinContent(i) + q1_2_1->GetBinContent(i) - q1_2_2->GetBinContent(i);

        k2sk_small[sample][i - 1] = q1_2_1->GetBinContent(i);
        k2_small[sample][i - 1] = Q2_small - powI(Q1_small, 2);

        #ifdef FILL_MC
          k2sk_small_gen[sample][i - 1] = N1p->GetBinContent(i);
          k2_small_gen[sample][i - 1] = N2->GetBinContent(i) - powI(N1->GetBinContent(i), 2);
        #endif // FILL_MC
      }

      for(int i = 1; i <= kNCentBins; i++)
      {
        k2sk[i - 1][sample] = cbwc(k2sk_small[sample], i - 1, hCent);
        k2[i - 1][sample] = cbwc(k2_small[sample], i - 1, hCent);

        #ifdef FILL_MC
          k2sk_gen[i - 1][sample] = cbwc(k2sk_small_gen[sample], i - 1, hCent);
          k2_gen[i - 1][sample] = cbwc(k2_small_gen[sample], i - 1, hCent);
        #endif // FILL_MC
      }

      fin->Close();
      fCent->Close();
      delete fin;
      delete fCent;
    }

    TGraphErrors g;
    g.SetName(Form("g_%d", iVar));
    #ifdef FILL_MC
      TGraphErrors g_gen;
      g_gen.SetName(Form("g_gen_%d", iVar));
    #endif // FILL_MC

    for(int i = 1; i < kNCentBins; i++)
    {
      double mean = 0.0;
      double rms = 0.0;
      k2_to_k2sk(mean, rms, k2sk[i - 1], k2[i - 1], nSkip);

      g.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), mean);
      hSys[i - 1]->Fill(mean);
      g.SetPointError(i - 1, 0, TMath::Sqrt(rms / (( N_SAMPLE - nSkip) * (( N_SAMPLE - nSkip) - 1))));

      #ifdef FILL_MC
        double mean_gen = 0.0;
        double rms_gen = 0.0;
        k2_to_k2sk(mean_gen, rms_gen, k2sk_gen[i - 1], k2_gen[i - 1], nSkip);

        g_gen.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), mean_gen);
        g_gen.SetPointError(i - 1, 0, TMath::Sqrt(rms_gen / (( N_SAMPLE - nSkip) * (( N_SAMPLE - nSkip) - 1))));
      #endif // FILL_MC
    }

    TCanvas c(Form("c_%d", iVar), Form("c_%d", iVar));
    c.cd();
    f.cd();
    g.Write();

    #ifdef FILL_MC
      g_gen.Write();
    #endif // FILL_MC
  }

  for (int i{0}; i < kNCentBins - 1; ++i){
    // remove_outlier(hSys[i]);
    hSys[i]->Write();
    hSys[i]->SetFillStyle(3004);
    hSys[i]->SetLineColor(kBlue);
    hSys[i]->SetFillColor(kBlue);
    hSys[i]->GetXaxis()->SetRangeUser(hSys[i]->GetMean() - 5.* hSys[i]->GetStdDev(), hSys[i]->GetMean() + 5.* hSys[i]->GetStdDev());
    hSys[i]->SetTitle(Form("mult. class %d", i));
    cSys.cd(i + 1);
    hSys[i]->Draw("histo");
  }

  cSys.Write();
  f.Close();
}