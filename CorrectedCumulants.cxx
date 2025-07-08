#include <TFile.h>
#include <TGraphErrors.h>

const char* dir = "../ResultsNetP/final_plots_08_hadPID";

const double v2TheFIST[]{0.0528284, 0.0637202, 0.0875567, 0.120773, 0.147558, 0.207523, 0.224011};
const double nch[]{18.68, 12.90, 10.03, 7.95, 6.32, 4.49, 2.54};

void CorrectedCumulants(){
  TFile *f_k4k2 = TFile::Open(Form("%s/out_sys_MB_18_08_finalBinning_singleParticleHighOrder_k2k1.root", dir));
  TFile *f_k2 = TFile::Open(Form("%s/out_sys_MB_18_08_finalBinning_singleParticleHighOrder_k1.root", dir));
  TFile *f_k2k1_ap = TFile::Open("/home/mciacco/out_sys_18_binning_fine_3_mix_finalBinning_singleParticleHighOrder_k2k1.root");
  TFile *f_k1_ap = TFile::Open("out_sys_18_binning_mix_finalBinning_singleParticleHighOrder_k1.root");

  TFile *f_k4k2_m = TFile::Open("out_sys_18_HM_mix_finalBinning_k4k2.root");
  TFile *f_k2sk_m = TFile::Open("out_sys_18_HM_mix_finalBinning_k2sk.root");

  TGraphErrors* k4k2 = static_cast<TGraphErrors*>(f_k4k2->Get("g_364"));
  TGraphErrors* k2 = static_cast<TGraphErrors*>(f_k2->Get("g_81"));
  TGraphErrors* k2k1_ap = static_cast<TGraphErrors*>(f_k2k1_ap->Get("g_81"));
  TGraphErrors* k1_ap = static_cast<TGraphErrors*>(f_k1_ap->Get("g_81"));

  TGraphErrors* k4k2_m = static_cast<TGraphErrors*>(f_k4k2_m->Get("g_81"));
  TGraphErrors* k2sk_m = static_cast<TGraphErrors*>(f_k2sk_m->Get("g_81"));

  TGraph v2fist(7, nch, v2TheFIST);

  for (int iP{0}; iP < k4k2->GetN(); ++iP){
    double k4k2_tmp = k4k2->GetPointY(iP);
    double k2_tmp = k2->GetPointY(iP);
    double k2k1_ap_tmp = k2k1_ap->GetPointY(iP);
    double k1_ap_tmp = k1_ap->GetPointY(iP);

    double k4k2_tmp_err = k4k2->GetErrorY(iP);
    double k2_tmp_err = k2->GetErrorY(iP);
    double k2k1_ap_tmp_err = k2k1_ap->GetErrorY(iP);
    double k1_ap_tmp_err = k1_ap->GetErrorY(iP);

    double k4k2_m_tmp = k4k2_m->GetPointY(iP);
    double k2sk_m_tmp = k2sk_m->GetPointY(iP);

    k4k2->SetPointY(iP,  k4k2_tmp -  k2_tmp *  (2. * k2k1_ap_tmp));
//    k4k2->SetPointError(iP, 0., k2k1_ap_tmp_err / k1_ap_tmp);
//    k4k2->SetPointY(iP, (k4k2_m_tmp - 1.) / k2sk_m_tmp / 3.);
  }
  TFile *fo = TFile::Open("fo_corr.root", "recreate");
  fo->cd();
  k4k2->Write();
  v2fist.Write();
  fo->Close();
}



/*
#include <TFile.h>
#include <TGraphErrors.h>

const char* dir = "../ResultsNetP/final_plots_08_hadPID";

void CorrectedCumulants(){
  TFile *f_k2k1 = TFile::Open(Form("%s/out_sys_MB_18_08_finalBinning_singleParticleHighOrder_k2k1.root", dir));
  TFile *f_k1 = TFile::Open(Form("%s/out_sys_MB_18_08_finalBinning_singleParticleHighOrder_k1.root", dir));
  TFile *f_k2k1_m = TFile::Open("out_sys_18_HM_mix_finalBinning_singleParticleHighOrder_k2k1.root");
  TFile *f_k1_m = TFile::Open("out_sys_18_HM_mix_finalBinning_singleParticleHighOrder_k1.root");

  TGraphErrors* k2k1 = static_cast<TGraphErrors*>(f_k2k1->Get("g_81"));
  TGraphErrors* k1 = static_cast<TGraphErrors*>(f_k1->Get("g_81"));
  TGraphErrors* k2k1_m = static_cast<TGraphErrors*>(f_k2k1_m->Get("g_81"));
  TGraphErrors* k1_m = static_cast<TGraphErrors*>(f_k1_m->Get("g_81"));

  for (int iP{0}; iP < k2k1->GetN(); ++iP){
    double k2k1_tmp = k2k1->GetPointY(iP);
    double k1_tmp = k1->GetPointY(iP);
    double k2k1_m_tmp = k2k1_m->GetPointY(iP);
    double k1_m_tmp = k1_m->GetPointY(iP);

    double k2k1_tmp_err = k2k1->GetErrorY(iP);
    double k1_tmp_err = k1->GetErrorY(iP);
    double k2k1_m_tmp_err = k2k1_m->GetErrorY(iP);
    double k1_m_tmp_err = k1_m->GetErrorY(iP);

    k2k1->SetPointY(iP, k2k1_tmp - (k2k1_m_tmp - 1.));
  }
  TFile *fo = TFile::Open("fo_corr.root", "recreate");
  fo->cd();
  k2k1->Write();
  fo->Close();
}

*/
