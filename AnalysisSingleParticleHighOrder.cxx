#include "utils.h"

//#define FILL_MC

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
  double lowEdge = kMultV0M ? kCentBins[centbin] : kTrklBins[centbin];
  double upEdge = kMultV0M ? kCentBins[centbin + 1] : kTrklBins[centbin + 1];
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
//      std::cout << num[sample] << std::endl;
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

void cumulant(double &mean, double &rms, const double *denom, const double *num, const double nSkip = 0.){
  mean = 0.0;
  rms = 0.0;
  for(int sample = 0; sample < kNSample; sample++)
  {
    if (std::abs(denom[sample]) > 1.e-9) {
//      std::cout << num[sample] << std::endl;
      mean = mean + (num[sample]);
    }
  }
  mean = mean / ( kNSample - nSkip);
  for(int sample = 0; sample < kNSample; sample++)
  {
    if (std::abs(denom[sample]) > 1.e-9) {
      rms = rms + powI(mean - (num[sample]), 2);
    }
  }
}

void AnalysisSingleParticleHighOrder(const char* period = "18", const char* obs = "k2k1")
{
  TFile f(Form("out_sys_%s_finalBinning_singleParticleHighOrder_%s.root", period, obs), "recreate");
  std::string obs_str{obs};
  TH1D *hSys[(kMultV0M ? kNCentBins : kNTrklBins)];
  TCanvas cSys("cSys", "cSys");
  cSys.Divide(3, 3);
  for (int i{0}; i < (kMultV0M ? kNCentBins : kNTrklBins); ++i){
    hSys[i] = new TH1D(Form("hSys_%d", i), ";#kappa_{2}/#kappa_{1};Entries", 10000, -.5, 1.5);
  }

  for(int iVar = 0; iVar < 780; ++iVar)
  {
    bool inVars = false;
    for (int iV{0}; iV < kNVar; ++iV) {
      if (iVar == kVar[iV]) {
        inVars = true;
        break;
      }
    }
    if (!inVars) continue;

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
    double k1_small[kNSample][100];
    double k2sk_small[kNSample][100];
    double k2_small[kNSample][100];
    double k3_small[kNSample][100];
    double k4_small[kNSample][100];
    double k5_small[kNSample][100];
    double k6_small[kNSample][100];
    double k1[100][kNSample];
    double k2sk[100][kNSample];
    double k2[100][kNSample];
    double k3[100][kNSample];
    double k4[100][kNSample];
    double k5[100][kNSample];
    double k6[100][kNSample];

    for(int sample = 0; sample < kNSample; sample++)
    {
      // if (sample == 0) {nSkip++; continue;}
      const char* triggerDir = kTriggerSel == 0x1 ? "MB_08" : "HM_08";
      TFile *fin = new TFile(Form("%s/%s/output_sys_singleParticle_%d_%d.root", kResDir, triggerDir, sample, iVar));
      TFile *fCent = TFile::Open(Form("%s/%s/LHC18ppTrig_HM%d_var_%d.root", kResDir, triggerDir, sample, iVar));

      TH1D *hCent = nullptr, *hNtrkl = nullptr;
      if (fCent){
        hCent = (TH1D*)fCent->Get(Form("hCent_%d", sample));
        hNtrkl = (TH1D*)fCent->Get(Form("hNtrkl_%d", sample));
      }
      else {
        hCent = new TH1D("hCent", "hCent", kNCentBinsSmall, kCentBinsSmall);
        hNtrkl = new TH1D("hNtrkl", "hNtrkl", kNTrklBinsSmall, kTrklBinsSmall);
        for (int iB{1}; iB <= hCent->GetNbinsX(); ++iB) hCent->SetBinContent(iB, 1);
        for (int iB{1}; iB <= hNtrkl->GetNbinsX(); ++iB) hNtrkl->SetBinContent(iB, 1);
      }
      if (!fin /* || sample == 4 */) {nSkip++; fin->Close(); delete fin; continue;}

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

        k1_small[sample][i - 1] = q1_1_1->GetBinContent(i);
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

      for(int i = 1; i <= (kMultV0M ? kNCentBins : kNTrklBins); i++)
      {
        k1[i - 1][sample] = cbwc(k1_small[sample], i - 1, kMultV0M ? hCent : hNtrkl);
        k2sk[i - 1][sample] = cbwc(k2sk_small[sample], i - 1, kMultV0M ? hCent : hNtrkl);
        k2[i - 1][sample] = cbwc(k2_small[sample], i - 1, kMultV0M ? hCent : hNtrkl);
        k3[i - 1][sample] = cbwc(k3_small[sample], i - 1, kMultV0M ? hCent : hNtrkl);
        k4[i - 1][sample] = cbwc(k4_small[sample], i - 1, kMultV0M ? hCent : hNtrkl);
        k5[i - 1][sample] = cbwc(k5_small[sample], i - 1, kMultV0M ? hCent : hNtrkl);
        k6[i - 1][sample] = cbwc(k6_small[sample], i - 1, kMultV0M ? hCent : hNtrkl);

        #ifdef FILL_MC
          k2sk_gen[i - 1][sample] = cbwc(k2sk_small_gen[sample], i - 1, kMultV0M ? hCent : hNtrkl);
          k2_gen[i - 1][sample] = cbwc(k2_small_gen[sample], i - 1, kMultV0M ? hCent : hNtrkl);
          k3_gen[i - 1][sample] = cbwc(k3_small_gen[sample], i - 1, kMultV0M ? hCent : hNtrkl);
          k4_gen[i - 1][sample] = cbwc(k4_small_gen[sample], i - 1, kMultV0M ? hCent : hNtrkl);
          k5_gen[i - 1][sample] = cbwc(k5_small_gen[sample], i - 1, kMultV0M ? hCent : hNtrkl);
          k6_gen[i - 1][sample] = cbwc(k6_small_gen[sample], i - 1, kMultV0M ? hCent : hNtrkl);
        #endif // FILL_MC
      }

      fin->Close();
      if (fCent)
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
    #endif // FILL_MC;

    double multHM[]{31.25, -1};
    double mult[]{18.63, 12.90, 10.03, 7.95, 6.32, 4.49, 2.54};
    //double mult[]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32.5,37.5};
    int nPoints = kMultV0M ? kNCentBins : kNTrklBins;
    if (kMultV0M && kTriggerSel == 0x2) nPoints -= 1;
    for(int i = 1; i <= nPoints; i++)
    {
      double mean = 0.0;
      double rms = 0.0;
      //cumulant_ratio(mean, rms, k2sk[i - 1], k2[i - 1], nSkip);
      if (obs_str == "k2k1") cumulant_ratio(mean, rms, k1[i - 1], k2[i - 1], nSkip);
      else if (obs_str == "k3k1") cumulant_ratio(mean, rms, k1[i - 1], k3[i - 1], nSkip);
      else if (obs_str == "k1") cumulant(mean, rms, k1[i - 1], k1[i - 1], nSkip);
      else if (obs_str == "k2") cumulant(mean, rms, k1[i - 1], k2[i - 1], nSkip);
      else if (obs_str == "k3") cumulant(mean, rms, k1[i - 1], k3[i - 1], nSkip);
      else if (obs_str == "k4") cumulant(mean, rms, k1[i - 1], k4[i - 1], nSkip);
      else if (obs_str == "k5") cumulant(mean, rms, k1[i - 1], k5[i - 1], nSkip);
      else if (obs_str == "k6") cumulant(mean, rms, k1[i - 1], k6[i - 1], nSkip);

      g.AddPoint(kTriggerSel == 0x1 ? mult[i - 1] : multHM[i - 1]/* 0.5 * (kCentBins[i - 1] + kCentBins[i]) */, mean);
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

  int nPoints = kMultV0M ? kNCentBins : kNTrklBins;
  if (kMultV0M && kTriggerSel == 0x2) nPoints -= 1;
  for (int i{0}; i < nPoints; ++i){
    remove_outlier(hSys[i]);
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
