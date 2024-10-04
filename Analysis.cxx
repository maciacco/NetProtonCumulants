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

const int kNSample = 10;

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

void cumulant_ratio(double &mean, double &rms, const double *denom, const double *num, const double nSkip = 0.){
  mean = 0.0;
  rms = 0.0;
  for(int sample = 0; sample < kNSample; sample++)
  {
    if (std::abs(denom[sample]) > 1.e-9) {
      mean = mean + (num[sample] / denom[sample]);
    }
  }
  mean = mean / ( kNSample - nSkip);
  for(int sample = 0; sample < kNSample; sample++)
  {
    if (std::abs(denom[sample]) > 1.e-9) {
      rms = rms + powI(mean - (num[sample] / denom[sample]), 2);
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
    hSys[i] = new TH1D(Form("hSys_%d", i), ";#kappa_{2}/#kappa_{1};Entries", 500, -2., 2.);
  }

  for(int iVar = 364; iVar < 365; ++iVar)
  {
    std::cout << "var = " << iVar << "..." << std::endl;

    double iB0 = 0.;
    double iB1 = 0.;
    double nSkip = 0;

    #ifdef FILL_MC
      double k2sk_small_gen[kNSample][100];
      double k2_small_gen[kNSample][100];
      double k3_small_gen[kNSample][100];
      double k4_small_gen[kNSample][100];
      double k5_small_gen[kNSample][100];
      double k6_small_gen[kNSample][100];
      double k2sk_gen[10][kNSample];
      double k2_gen[10][kNSample];
      double k3_gen[10][kNSample];
      double k4_gen[10][kNSample];
      double k5_gen[10][kNSample];
      double k6_gen[10][kNSample];
    #endif

    double Q1_small;
    double Q2_small;
    double Q3_small;
    double Q4_small;
    double Q5_small;
    double Q6_small;
    double k2sk_small[kNSample][100];
    double k2_small[kNSample][100];
    double k3_small[kNSample][100];
    double k4_small[kNSample][100];
    double k5_small[kNSample][100];
    double k6_small[kNSample][100];
    double k2sk[10][kNSample];
    double k2[10][kNSample];
    double k3[10][kNSample];
    double k4[10][kNSample];
    double k5[10][kNSample];
    double k6[10][kNSample];

    for(int sample = 0; sample < kNSample; sample++)
    {
      // if (sample == 0) {nSkip++; continue;}
      TFile *fin = new TFile(Form("%s/output_sys_%d_%d.root", kResDir, sample, iVar));
      TFile *fCent = TFile::Open(Form("%s/LHC21d3%d_var_%d.root", kResDir, sample, iVar));

      TH1D *hCent = (TH1D*)fCent->Get(Form("hCent_%d", sample));

      if (!fin || sample == 4) {nSkip++; fin->Close(); delete fin; continue;}

      // 1st order
      TProfile *q1_1_1 = (TProfile*)fin->Get(Form("var_%d/q1_1_1", iVar));

      // 2nd order
      TProfile *q2_1_1 = (TProfile*)fin->Get(Form("var_%d/q2_1_1", iVar));
      TProfile *q1_2_1 = (TProfile*)fin->Get(Form("var_%d/q1_2_1", iVar));
      TProfile *q1_2_2 = (TProfile*)fin->Get(Form("var_%d/q1_2_2", iVar));

      // 3rd order
      TProfile *q3_1_1 = (TProfile*)fin->Get(Form("var_%d/q3_1_1", iVar));
      TProfile *q1_1_1_x_q1_2_1 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_2_1", iVar));
      TProfile *q1_1_1_x_q1_2_2 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_2_2", iVar));
      TProfile *q1_3_1 = (TProfile*)fin->Get(Form("var_%d/q1_3_1", iVar));
      TProfile *q1_3_2 = (TProfile*)fin->Get(Form("var_%d/q1_3_2", iVar));
      TProfile *q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q1_3_3", iVar));

      // 4th order
      TProfile *q4_1_1 = (TProfile*)fin->Get(Form("var_%d/q4_1_1", iVar));
      TProfile *q2_1_1_x_q1_2_1 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_2_1", iVar));
      TProfile *q2_1_1_x_q1_2_2 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_2_2", iVar));
      TProfile *q1_1_1_x_q1_3_1 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_3_1", iVar));
      TProfile *q2_2_1 = (TProfile*)fin->Get(Form("var_%d/q2_2_1", iVar));
      TProfile *q2_2_2 = (TProfile*)fin->Get(Form("var_%d/q2_2_2", iVar));
      TProfile *q1_1_1_x_q1_3_2 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_3_2", iVar));
      TProfile *q1_1_1_x_q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_3_3", iVar));
      TProfile *q1_2_1_x_q1_2_2 = (TProfile*)fin->Get(Form("var_%d/q1_2_1_x_q1_2_2", iVar));
      TProfile *q1_4_1 = (TProfile*)fin->Get(Form("var_%d/q1_4_1", iVar));
      TProfile *q1_4_2 = (TProfile*)fin->Get(Form("var_%d/q1_4_2", iVar));
      TProfile *q1_4_3 = (TProfile*)fin->Get(Form("var_%d/q1_4_3", iVar));
      TProfile *q1_4_4 = (TProfile*)fin->Get(Form("var_%d/q1_4_4", iVar));

      // 5th order
      TProfile *q5_1_1 = (TProfile*)fin->Get(Form("var_%d/q5_1_1", iVar));
      TProfile *q3_1_1_x_q1_2_1 = (TProfile*)fin->Get(Form("var_%d/q3_1_1_x_q1_2_1", iVar));
      TProfile *q3_1_1_x_q1_2_2 = (TProfile*)fin->Get(Form("var_%d/q3_1_1_x_q1_2_2", iVar));
      TProfile *q2_1_1_x_q1_3_1 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_3_1", iVar));
      TProfile *q2_1_1_x_q1_3_2 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_3_2", iVar));
      TProfile *q2_1_1_x_q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_3_3", iVar));
      TProfile *q2_2_2_x_q1_1_1 = (TProfile*)fin->Get(Form("var_%d/q2_2_2_x_q1_1_1", iVar));
      TProfile *q2_2_1_x_q1_1_1 = (TProfile*)fin->Get(Form("var_%d/q2_2_1_x_q1_1_1", iVar));
      TProfile *q1_1_1_x_q1_2_1_x_q1_2_2 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_2_1_x_q1_2_2", iVar));
      TProfile *q1_1_1_x_q1_4_1 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_4_1", iVar));
      TProfile *q1_1_1_x_q1_4_2 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_4_2", iVar));
      TProfile *q1_1_1_x_q1_4_3 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_4_3", iVar));
      TProfile *q1_1_1_x_q1_4_4 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_4_4", iVar));
      TProfile *q1_2_1_x_q1_3_1 = (TProfile*)fin->Get(Form("var_%d/q1_2_1_x_q1_3_1", iVar));
      TProfile *q1_2_1_x_q1_3_2 = (TProfile*)fin->Get(Form("var_%d/q1_2_1_x_q1_3_2", iVar));
      TProfile *q1_2_1_x_q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q1_2_1_x_q1_3_3", iVar));
      TProfile *q1_2_2_x_q1_3_1 = (TProfile*)fin->Get(Form("var_%d/q1_2_2_x_q1_3_1", iVar));
      TProfile *q1_2_2_x_q1_3_2 = (TProfile*)fin->Get(Form("var_%d/q1_2_2_x_q1_3_2", iVar));
      TProfile *q1_2_2_x_q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q1_2_2_x_q1_3_3", iVar));
      TProfile *q1_5_1 = (TProfile*)fin->Get(Form("var_%d/q1_5_1", iVar));
      TProfile *q1_5_2 = (TProfile*)fin->Get(Form("var_%d/q1_5_2", iVar));
      TProfile *q1_5_3 = (TProfile*)fin->Get(Form("var_%d/q1_5_3", iVar));
      TProfile *q1_5_4 = (TProfile*)fin->Get(Form("var_%d/q1_5_4", iVar));
      TProfile *q1_5_5 = (TProfile*)fin->Get(Form("var_%d/q1_5_5", iVar));

      // 6th order
      TProfile *q6_1_1 = (TProfile*)fin->Get(Form("var_%d/q6_1_1", iVar));
      TProfile *q4_1_1_x_q1_2_1 = (TProfile*)fin->Get(Form("var_%d/q4_1_1_x_q1_2_1", iVar));
      TProfile *q4_1_1_x_q1_2_2 = (TProfile*)fin->Get(Form("var_%d/q4_1_1_x_q1_2_2", iVar));
      TProfile *q3_1_1_x_q1_3_1 = (TProfile*)fin->Get(Form("var_%d/q3_1_1_x_q1_3_1", iVar));
      TProfile *q3_1_1_x_q1_3_2 = (TProfile*)fin->Get(Form("var_%d/q3_1_1_x_q1_3_2", iVar));
      TProfile *q3_1_1_x_q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q3_1_1_x_q1_3_3", iVar));
      TProfile *q2_1_1_x_q1_2_2_x_q1_2_1 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_2_2_x_q1_2_1", iVar));
      TProfile *q2_1_1_x_q2_2_1 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q2_2_1", iVar));
      TProfile *q2_1_1_x_q2_2_2 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q2_2_2", iVar));
      TProfile *q3_2_1 = (TProfile*)fin->Get(Form("var_%d/q3_2_1", iVar));
      TProfile *q3_2_2 = (TProfile*)fin->Get(Form("var_%d/q3_2_2", iVar));
      TProfile *q2_1_1_x_q1_4_1 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_4_1", iVar));
      TProfile *q2_1_1_x_q1_4_2 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_4_2", iVar));
      TProfile *q2_1_1_x_q1_4_3 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_4_3", iVar));
      TProfile *q2_1_1_x_q1_4_4 = (TProfile*)fin->Get(Form("var_%d/q2_1_1_x_q1_4_4", iVar));
      TProfile *q2_2_1_x_q1_2_2 = (TProfile*)fin->Get(Form("var_%d/q2_2_1_x_q1_2_2", iVar));
      TProfile *q2_2_2_x_q1_2_1 = (TProfile*)fin->Get(Form("var_%d/q2_2_2_x_q1_2_1", iVar));
      TProfile *q1_1_1_x_q1_2_1_x_q1_3_1 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_2_1_x_q1_3_1", iVar));
      TProfile *q1_1_1_x_q1_2_1_x_q1_3_2 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_2_1_x_q1_3_2", iVar));
      TProfile *q1_1_1_x_q1_2_1_x_q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_2_1_x_q1_3_3", iVar));
      TProfile *q1_1_1_x_q1_2_2_x_q1_3_1 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_2_2_x_q1_3_1", iVar));
      TProfile *q1_1_1_x_q1_2_2_x_q1_3_2 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_2_2_x_q1_3_2", iVar));
      TProfile *q1_1_1_x_q1_2_2_x_q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_2_2_x_q1_3_3", iVar));
      TProfile *q1_1_1_x_q1_5_1 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_5_1", iVar));
      TProfile *q1_1_1_x_q1_5_2 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_5_2", iVar));
      TProfile *q1_1_1_x_q1_5_3 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_5_3", iVar));
      TProfile *q1_1_1_x_q1_5_4 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_5_4", iVar));
      TProfile *q1_1_1_x_q1_5_5 = (TProfile*)fin->Get(Form("var_%d/q1_1_1_x_q1_5_5", iVar));
      TProfile *q1_2_1_x_q1_4_1 = (TProfile*)fin->Get(Form("var_%d/q1_2_1_x_q1_4_1", iVar));
      TProfile *q1_2_1_x_q1_4_2 = (TProfile*)fin->Get(Form("var_%d/q1_2_1_x_q1_4_2", iVar));
      TProfile *q1_2_1_x_q1_4_3 = (TProfile*)fin->Get(Form("var_%d/q1_2_1_x_q1_4_3", iVar));
      TProfile *q1_2_1_x_q1_4_4 = (TProfile*)fin->Get(Form("var_%d/q1_2_1_x_q1_4_4", iVar));
      TProfile *q1_2_2_x_q1_4_1 = (TProfile*)fin->Get(Form("var_%d/q1_2_2_x_q1_4_1", iVar));
      TProfile *q1_2_2_x_q1_4_2 = (TProfile*)fin->Get(Form("var_%d/q1_2_2_x_q1_4_2", iVar));
      TProfile *q1_2_2_x_q1_4_3 = (TProfile*)fin->Get(Form("var_%d/q1_2_2_x_q1_4_3", iVar));
      TProfile *q1_2_2_x_q1_4_4 = (TProfile*)fin->Get(Form("var_%d/q1_2_2_x_q1_4_4", iVar));
      TProfile *q2_3_1 = (TProfile*)fin->Get(Form("var_%d/q2_3_1", iVar));
      TProfile *q1_3_1_x_q1_3_2 = (TProfile*)fin->Get(Form("var_%d/q1_3_1_x_q1_3_2", iVar));
      TProfile *q1_3_1_x_q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q1_3_1_x_q1_3_3", iVar));
      TProfile *q2_3_2 = (TProfile*)fin->Get(Form("var_%d/q2_3_2", iVar));
      TProfile *q1_3_2_x_q1_3_3 = (TProfile*)fin->Get(Form("var_%d/q1_3_2_x_q1_3_3", iVar));
      TProfile *q2_3_3 = (TProfile*)fin->Get(Form("var_%d/q2_3_3", iVar));
      TProfile *q1_6_1 = (TProfile*)fin->Get(Form("var_%d/q1_6_1", iVar));
      TProfile *q1_6_2 = (TProfile*)fin->Get(Form("var_%d/q1_6_2", iVar));
      TProfile *q1_6_3 = (TProfile*)fin->Get(Form("var_%d/q1_6_3", iVar));
      TProfile *q1_6_4 = (TProfile*)fin->Get(Form("var_%d/q1_6_4", iVar));
      TProfile *q1_6_5 = (TProfile*)fin->Get(Form("var_%d/q1_6_5", iVar));
      TProfile *q1_6_6 = (TProfile*)fin->Get(Form("var_%d/q1_6_6", iVar));

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
        Q3_small = q3_1_1->GetBinContent(i) + (3. * q1_1_1_x_q1_2_1->GetBinContent(i)) - (3. * q1_1_1_x_q1_2_2->GetBinContent(i)) + q1_3_1->GetBinContent(i) - (3. * q1_3_2->GetBinContent(i)) + (2. * q1_3_3->GetBinContent(i));
        Q4_small = q4_1_1->GetBinContent(i) + (6. * q2_1_1_x_q1_2_1->GetBinContent(i)) - (6. * q2_1_1_x_q1_2_2->GetBinContent(i)) + (4. * q1_1_1_x_q1_3_1->GetBinContent(i)) + (3. * q2_2_1->GetBinContent(i))
                   + (3. * q2_2_2->GetBinContent(i)) - (12. * (q1_1_1_x_q1_3_2->GetBinContent(i))) + (8. * q1_1_1_x_q1_3_3->GetBinContent(i)) - (6. * q1_2_1_x_q1_2_2->GetBinContent(i))
                   + (q1_4_1->GetBinContent(i)) - (7. * q1_4_2->GetBinContent(i)) + (12. * q1_4_3->GetBinContent(i)) - (6. * q1_4_4->GetBinContent(i));
        Q5_small = q5_1_1->GetBinContent(i) + (10. * q3_1_1_x_q1_2_1->GetBinContent(i)) - (10. * q3_1_1_x_q1_2_2->GetBinContent(i)) + (10. * q2_1_1_x_q1_3_1->GetBinContent(i)) - (30. * q2_1_1_x_q1_3_2->GetBinContent(i))
                   + (20. * q2_1_1_x_q1_3_3->GetBinContent(i)) + (15. * q2_2_2_x_q1_1_1->GetBinContent(i)) + (15. * q2_2_1_x_q1_1_1->GetBinContent(i)) - (30. * q1_1_1_x_q1_2_1_x_q1_2_2->GetBinContent(i))
                   + (5. * q1_1_1_x_q1_4_1->GetBinContent(i)) - (35. * q1_1_1_x_q1_4_2->GetBinContent(i)) + (60. * q1_1_1_x_q1_4_3->GetBinContent(i)) - (30. * q1_1_1_x_q1_4_4->GetBinContent(i))
                   + (10. * q1_2_1_x_q1_3_1->GetBinContent(i)) - (30. * q1_2_1_x_q1_3_2->GetBinContent(i)) + (20. * q1_2_1_x_q1_3_3->GetBinContent(i))
                   - (10. * q1_2_2_x_q1_3_1->GetBinContent(i)) + (30. * q1_2_2_x_q1_3_2->GetBinContent(i)) - (20. * q1_2_2_x_q1_3_3->GetBinContent(i))
                   + q1_5_1->GetBinContent(i) - (15. * q1_5_2->GetBinContent(i)) + (50. * q1_5_3->GetBinContent(i)) - (60. * q1_5_4->GetBinContent(i)) + (24. * q1_5_5->GetBinContent(i));
        Q6_small = q6_1_1->GetBinContent(i) + (15. * q4_1_1_x_q1_2_1->GetBinContent(i)) - (15. * q4_1_1_x_q1_2_2->GetBinContent(i)) + (20. * q3_1_1_x_q1_3_1->GetBinContent(i)) - (60. * q3_1_1_x_q1_3_2->GetBinContent(i))
                   + (40. * q3_1_1_x_q1_3_3->GetBinContent(i)) - (90. * q2_1_1_x_q1_2_2_x_q1_2_1->GetBinContent(i)) + (45. * q2_1_1_x_q2_2_1->GetBinContent(i)) + (45. * q2_1_1_x_q2_2_2->GetBinContent(i))
                   + (15. * q3_2_1->GetBinContent(i)) - (15. * q3_2_2->GetBinContent(i)) + (15. * q2_1_1_x_q1_4_1->GetBinContent(i)) - (105. * q2_1_1_x_q1_4_2->GetBinContent(i)) + (180. * q2_1_1_x_q1_4_3->GetBinContent(i)) - (90. * q2_1_1_x_q1_4_4->GetBinContent(i))
                   - (45. * q2_2_1_x_q1_2_2->GetBinContent(i)) + (45. * q2_2_2_x_q1_2_1->GetBinContent(i)) + (60. * q1_1_1_x_q1_2_1_x_q1_3_1->GetBinContent(i)) - (180. * q1_1_1_x_q1_2_1_x_q1_3_2->GetBinContent(i))
                   + (120. * q1_1_1_x_q1_2_1_x_q1_3_3->GetBinContent(i)) - (60. * q1_1_1_x_q1_2_2_x_q1_3_1->GetBinContent(i)) + (180. * q1_1_1_x_q1_2_2_x_q1_3_2->GetBinContent(i)) - (120. * q1_1_1_x_q1_2_2_x_q1_3_3->GetBinContent(i))
                   + (6. * q1_1_1_x_q1_5_1->GetBinContent(i)) - (90. * q1_1_1_x_q1_5_2->GetBinContent(i)) + (300. * q1_1_1_x_q1_5_3->GetBinContent(i)) - (360. * q1_1_1_x_q1_5_4->GetBinContent(i)) + (144. * q1_1_1_x_q1_5_5->GetBinContent(i))
                   + (15. * q1_2_1_x_q1_4_1->GetBinContent(i)) - (105. * q1_2_1_x_q1_4_2->GetBinContent(i)) + (180. * q1_2_1_x_q1_4_3->GetBinContent(i)) - (90. * q1_2_1_x_q1_4_4->GetBinContent(i))
                   - (15. * q1_2_2_x_q1_4_1->GetBinContent(i)) + (105. * q1_2_2_x_q1_4_2->GetBinContent(i)) - (180. * q1_2_2_x_q1_4_3->GetBinContent(i)) + (90. * q1_2_2_x_q1_4_4->GetBinContent(i))
                   + (10. * q2_3_1->GetBinContent(i)) - (60. * q1_3_1_x_q1_3_2->GetBinContent(i)) + (40. * q1_3_1_x_q1_3_3->GetBinContent(i)) + (90. * q2_3_2->GetBinContent(i)) - (120. * q1_3_2_x_q1_3_3->GetBinContent(i)) + (40. * q2_3_3->GetBinContent(i))
                   + q1_6_1->GetBinContent(i) - (31. * q1_6_2->GetBinContent(i)) + (180. * q1_6_3->GetBinContent(i)) - (390. * q1_6_4->GetBinContent(i)) + (360. * q1_6_5->GetBinContent(i)) - (120. * q1_6_6->GetBinContent(i));
        // std::cout << Q3_small << std::endl;

        k2sk_small[sample][i - 1] = q1_2_1->GetBinContent(i);
        k2_small[sample][i - 1] = Q2_small - powI(Q1_small, 2);
        k3_small[sample][i - 1] = Q3_small - (3. * Q2_small * Q1_small) + (2. * powI(Q1_small, 3));
        k4_small[sample][i - 1] = Q4_small - (4. * Q3_small * Q1_small) - (3. * powI(Q2_small, 2)) + (12. * Q2_small * powI(Q1_small, 2)) - (6. * powI(Q1_small, 4));
        k5_small[sample][i - 1] = Q5_small - (5. * Q4_small * Q1_small) - (10. * Q3_small * Q2_small) + (20. * Q3_small * powI(Q1_small, 2)) + (30. * powI(Q2_small, 2) * Q1_small) - (60. * Q2_small * powI(Q1_small, 3)) + (24. * powI(Q1_small, 5));
        k6_small[sample][i - 1] = Q6_small - (6. * Q5_small * Q1_small) - (15. * Q4_small * Q2_small) + (30. * Q4_small * powI(Q1_small, 2)) - (10. * powI(Q3_small, 2)) + (120. * Q3_small * Q2_small * Q1_small)
                                  - (120. * Q3_small * powI(Q1_small, 3)) + (30. * powI(Q2_small, 3)) - (270. * powI(Q2_small, 2) * powI(Q1_small, 2)) + (360. * Q2_small * powI(Q1_small, 4)) - (120. * powI(Q1_small, 6));

        #ifdef FILL_MC
          k2sk_small_gen[sample][i - 1] = N1p->GetBinContent(i);
          k2_small_gen[sample][i - 1] = N2->GetBinContent(i) - powI(N1->GetBinContent(i), 2);
          k3_small_gen[sample][i - 1] = N3->GetBinContent(i) - (3. * N2->GetBinContent(i) * N1->GetBinContent(i)) + (2. * powI(N1->GetBinContent(i), 3));
          k4_small_gen[sample][i - 1] = N4->GetBinContent(i) - (4. * N3->GetBinContent(i) * N1->GetBinContent(i)) - (3. * powI(N2->GetBinContent(i), 2)) + (12. * N2->GetBinContent(i) * powI(N1->GetBinContent(i), 2)) - (6. * powI(N1->GetBinContent(i), 4));
          k5_small_gen[sample][i - 1] = N5->GetBinContent(i) - (5. * N4->GetBinContent(i) * N1->GetBinContent(i)) - (10. * N3->GetBinContent(i) * N2->GetBinContent(i)) + (20. * N3->GetBinContent(i) * powI(N1->GetBinContent(i), 2)) + (30. * powI(N2->GetBinContent(i), 2) * N1->GetBinContent(i)) - (60. * N2->GetBinContent(i) * powI(N1->GetBinContent(i), 3)) + (24. * powI(N1->GetBinContent(i), 5));
          k6_small_gen[sample][i - 1] = N6->GetBinContent(i) - (6. * N5->GetBinContent(i) * N1->GetBinContent(i)) - (15. * N4->GetBinContent(i) * N2->GetBinContent(i)) + (30. * N4->GetBinContent(i) * powI(N1->GetBinContent(i), 2)) - (10. * powI(N3->GetBinContent(i), 2)) + (120. * N3->GetBinContent(i) * N2->GetBinContent(i) * N1->GetBinContent(i))
                                        - (120. * N3->GetBinContent(i) * powI(N1->GetBinContent(i), 3)) + (30. * powI(N2->GetBinContent(i), 3)) - (270. * powI(N2->GetBinContent(i), 2) * powI(N1->GetBinContent(i), 2)) + (360. * N2->GetBinContent(i) * powI(N1->GetBinContent(i), 4)) - (120. * powI(N1->GetBinContent(i), 6));

        #endif // FILL_MC
      }

      for(int i = 1; i <= kNCentBins; i++)
      {
        k2sk[i - 1][sample] = cbwc(k2sk_small[sample], i - 1, hCent);
        k2[i - 1][sample] = cbwc(k2_small[sample], i - 1, hCent);
        k3[i - 1][sample] = cbwc(k3_small[sample], i - 1, hCent);
        k4[i - 1][sample] = cbwc(k4_small[sample], i - 1, hCent);
        k5[i - 1][sample] = cbwc(k5_small[sample], i - 1, hCent);
        k6[i - 1][sample] = cbwc(k6_small[sample], i - 1, hCent);

        #ifdef FILL_MC
          k2sk_gen[i - 1][sample] = cbwc(k2sk_small_gen[sample], i - 1, hCent);
          k2_gen[i - 1][sample] = cbwc(k2_small_gen[sample], i - 1, hCent);
          k3_gen[i - 1][sample] = cbwc(k3_small_gen[sample], i - 1, hCent);
          k4_gen[i - 1][sample] = cbwc(k4_small_gen[sample], i - 1, hCent);
          k5_gen[i - 1][sample] = cbwc(k5_small_gen[sample], i - 1, hCent);
          k6_gen[i - 1][sample] = cbwc(k6_small_gen[sample], i - 1, hCent);
        #endif // FILL_MC
      }

      fin->Close();
      fCent->Close();
      delete fin;
      delete fCent;
    }

    TGraphErrors g;
    g.SetName(Form("g_%d", iVar));
    g.SetLineWidth(2);
    g.SetLineColor(kRed);
    g.SetMarkerColor(kRed);
    #ifdef FILL_MC
      TGraphErrors g_gen;
      g_gen.SetName(Form("g_gen_%d", iVar));
      g_gen.SetLineWidth(2);
      g_gen.SetLineColor(kBlue);
      g_gen.SetMarkerColor(kBlue);
    #endif // FILL_MC

    for(int i = 1; i <= kNCentBins; i++)
    {
      double mean = 0.0;
      double rms = 0.0;
      //cumulant_ratio(mean, rms, k2sk[i - 1], k2[i - 1], nSkip);
      cumulant_ratio(mean, rms, k2sk[i - 1], k2[i - 1], nSkip);

      g.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), mean);
      hSys[i - 1]->Fill(mean);
      g.SetPointError(i - 1, 0, TMath::Sqrt(rms / (( kNSample - nSkip) * (( kNSample - nSkip) - 1))));

      #ifdef FILL_MC
        double mean_gen = 0.0;
        double rms_gen = 0.0;
        // cumulant_ratio(mean_gen, rms_gen, k2sk_gen[i - 1], k2_gen[i - 1], nSkip);
        cumulant_ratio(mean_gen, rms_gen, k2sk_gen[i - 1], k2_gen[i - 1], nSkip);

        g_gen.AddPoint(0.5 * (kCentBins[i - 1] + kCentBins[i]), mean_gen);
        g_gen.SetPointError(i - 1, 0, TMath::Sqrt(rms_gen / (( kNSample - nSkip) * (( kNSample - nSkip) - 1))));
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