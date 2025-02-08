#include "SpectraUtils.h"

double mass = 0.93827208816;

void LTFit(){
  TFile *f_pr = TFile::Open("data/HEPData-ins1784041-v1-Table_5.root");
  TGraphErrors* g_cpy = new TGraphErrors;
  TGraphAsymmErrors* g = (TGraphAsymmErrors*)f_pr->Get("Table 5/Graph1D_y4");
  TGraphAsymmErrors* g_2 = (TGraphAsymmErrors*)f_pr->Get("Table 5/Graph1D_y5");
  //TGraphErrors* g_3 = (TGraphErrors*)f_pr->Get("Table 5/Graph1D_y3");

  // define ranges as in spectra paper
  for (int iP{0}; iP < g->GetN(); ++iP) { // pi
    g_cpy->AddPoint(g->GetPointX(iP),  0.5 * g->GetPointY(iP) + 0.5 * g_2->GetPointY(iP)); // + 0.5 * g_3->GetPointY(iP));
    g_cpy->SetPointError(g_cpy->GetN() - 1, g->GetErrorXhigh(iP), g->GetErrorYhigh(iP));
  }
  TF1 *lt = LevyTsallis("lt", mass, 5, 0.1, 0.2);
  g_cpy->Fit("lt", "RM+");
  std::cout << lt->Integral(0.5, 1.5, 1.e-5) / lt->Integral(0., 10., 1.e-5) << std::endl;
}

void LTFit(bool isHM){
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
        g_->AddPoint(x, p);
        g_->SetPointError(g_->GetN() - 1, 0.5 * std::abs(subs1 - subs2), std::hypot(est, esy));
        std::cout << subs1 << "\t" << subs2 << "\t" << p << "\t" << est << "\t" << esy << std::endl;
        fin1 >> s;
      }
    }
  }

  TGraphErrors* g_cpy = new TGraphErrors;
  for (int iP{2}; iP < 3; ++iP) {
    TGraphErrors* g_1 = dynamic_cast<TGraphErrors*>((*gr_dat[0])[iP]);
    TGraphErrors* g_2 = dynamic_cast<TGraphErrors*>((*gr_dat[1])[iP]);
    TGraphErrors* g_3 = dynamic_cast<TGraphErrors*>((*gr_dat[2])[iP]);
    for (int iP{0}; iP < g_1->GetN(); ++iP) {
      g_cpy->AddPoint(g_1->GetPointX(iP), 0.1 * g_1->GetPointY(iP) + 0.4 * g_2->GetPointY(iP) + 0.5 * g_3->GetPointY(iP));
      g_cpy->SetPointError(g_cpy->GetN() - 1, g_1->GetErrorX(iP), g_1->GetErrorY(iP));
    }
  }

  TF1 *lt = LevyTsallis("lt", mass);
  g_cpy->Fit("lt", "M+");
  std::cout << lt->Integral(0.5, 1.5, 1.e-5) / lt->Integral(0., 10., 1.e-5) << std::endl;
}
