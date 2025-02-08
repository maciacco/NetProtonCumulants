#include "SpectraUtils.h"

double limits_pt[][2]{{0.5, 1}, {0.2, 1.5}, {0.3, 3.0}};

void CombinedBWFit(){
  TFile *f_pi = TFile::Open("data/HEPData-ins1784041-v1-Table_1.root");
  TFile *f_ka = TFile::Open("data/HEPData-ins1784041-v1-Table_3.root");
  TFile *f_pr = TFile::Open("data/HEPData-ins1784041-v1-Table_5.root");
  TObjArray *gr_dat = new TObjArray(3);
  TObjArray *gr_dat_2 = new TObjArray(3);
  TObjArray *gr_dat_3 = new TObjArray(3);
  TObjArray* gr_dat_cpy = new TObjArray(3);
  (*gr_dat)[0] = (TGraphErrors*)f_pi->Get("Table 1/Graph1D_y4");
  (*gr_dat)[1] = (TGraphErrors*)f_ka->Get("Table 3/Graph1D_y4");
  (*gr_dat)[2] = (TGraphErrors*)f_pr->Get("Table 5/Graph1D_y4");
  (*gr_dat_2)[0] = (TGraphErrors*)f_pi->Get("Table 1/Graph1D_y5");
  (*gr_dat_2)[1] = (TGraphErrors*)f_ka->Get("Table 3/Graph1D_y5");
  (*gr_dat_2)[2] = (TGraphErrors*)f_pr->Get("Table 5/Graph1D_y5");
//   (*gr_dat_3)[0] = (TGraphErrors*)f_pi->Get("Table 1/Graph1D_y3");
//   (*gr_dat_3)[1] = (TGraphErrors*)f_ka->Get("Table 3/Graph1D_y3");
//   (*gr_dat_3)[2] = (TGraphErrors*)f_pr->Get("Table 5/Graph1D_y3");

  // define ranges as in spectra paper
  for (int iPart{0}; iPart < 3; ++iPart) {
    TGraphAsymmErrors* g = dynamic_cast<TGraphAsymmErrors*>((*gr_dat)[iPart]);
    TGraphAsymmErrors* g_2 = dynamic_cast<TGraphAsymmErrors*>((*gr_dat_2)[iPart]);
    // TGraphAsymmErrors* g_3 = dynamic_cast<TGraphAsymmErrors*>((*gr_dat_3)[iPart]);
    (*gr_dat_cpy)[iPart] = new TGraphErrors;
    TGraphErrors* g_cpy = dynamic_cast<TGraphErrors*>((*gr_dat_cpy)[iPart]);
    for (int iP{0}; iP < g->GetN(); ++iP) { // pi
      if (g->GetPointX(iP) > limits_pt[iPart][0] && g->GetPointX(iP) < limits_pt[iPart][1]) {
        g_cpy->AddPoint(g->GetPointX(iP), 0.5 * g->GetPointY(iP) + 0.5 * g_2->GetPointY(iP) /* + 0.5 * g_3->GetPointY(iP) */);
        g_cpy->SetPointError(g_cpy->GetN() - 1, g->GetErrorXhigh(iP), g->GetErrorYhigh(iP));
      }
    }
  }
  double masses[]{0.139, 0.495, 0.938}; // GeV/c^2
  BGBlastWave_GlobalFit(gr_dat_cpy, masses);
}

void CombinedBWFit(bool isHM){
  if (!isHM) return;

  TObjArray* gr_dat[3];

  for (int iF{0}; iF < 3; ++iF) {
    gr_dat[iF] = new TObjArray(3);
    (*gr_dat[iF])[0] = new TGraphErrors;
    (*gr_dat[iF])[1] = new TGraphErrors;
    (*gr_dat[iF])[2] = new TGraphErrors;
    std::ifstream fin1(Form("data/HM%d.txt", iF + 1));

    string s;
    for (int iP{0}; iP < 3; ++iP) {
      TGraphErrors* g_ = dynamic_cast<TGraphErrors*>((*gr_dat[iF])[iP]);
      g_->Clear();
      fin1 >> s;
      while (s != "SystError" && !fin1.eof()) {
        fin1 >> s;
        std::cout << s <<std::endl;
      }
      fin1 >> s;
      while (s.at(1) != '*' && !fin1.eof()) {
        auto it = s.find("-");
        double subs1 = std::stod(s.substr(0, it));
        double subs2 = std::stod(s.substr(it + 1, s.size()));
        fin1 >> s;
        double p = std::stod(s);
        fin1 >> s;
        double est = std::stod(s);
        fin1 >> s;
        double esy = std::stod(s);

        double x = 0.5 * (subs1 + subs2);
        if (x > limits_pt[iP][0] && x < limits_pt[iP][1]){
          g_->AddPoint(x, p);
          g_->SetPointError(g_->GetN() - 1, 0.5 * std::abs(subs1 - subs2), std::hypot(est, esy));
          std::cout << subs1 << "\t" << subs2 << "\t" << p << "\t" << est << "\t" << esy << std::endl;
        }
        fin1 >> s;
      }
    }
  }

  TObjArray* gr_dat_cpy = new TObjArray(3);
  (*gr_dat_cpy)[0] = new TGraphErrors;
  (*gr_dat_cpy)[1] = new TGraphErrors;
  (*gr_dat_cpy)[2] = new TGraphErrors;
  for (int iP{0}; iP < 3; ++iP) {
    TGraphErrors* g_cpy = dynamic_cast<TGraphErrors*>((*gr_dat_cpy)[iP]);
    TGraphErrors* g_1 = dynamic_cast<TGraphErrors*>((*gr_dat[0])[iP]);
    TGraphErrors* g_2 = dynamic_cast<TGraphErrors*>((*gr_dat[1])[iP]);
    TGraphErrors* g_3 = dynamic_cast<TGraphErrors*>((*gr_dat[2])[iP]);
    for (int iP{0}; iP < g_1->GetN(); ++iP) {
      g_cpy->AddPoint(g_1->GetPointX(iP), 0.1 * g_1->GetPointY(iP) + 0.4 * g_2->GetPointY(iP) + 0.5 * g_3->GetPointY(iP));
      g_cpy->SetPointError(g_cpy->GetN() - 1, g_1->GetErrorX(iP), g_1->GetErrorY(iP));
    }
  }

  double masses[]{0.139, 0.495, 0.938}; // GeV/c^2
  BGBlastWave_GlobalFit(gr_dat_cpy, masses);
}
