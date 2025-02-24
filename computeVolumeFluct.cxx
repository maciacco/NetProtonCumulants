void computeVolumeFluct(){
  TFile f("fo_compute_volF.root", "recreate");

  TFile fin_k2k1_n("../ResultsNetP/final_plots_08_hadPID/out_sys_MB_18_08_finalBinning_singleParticleHighOrder_k2k1.root");
  TFile fin_k3k1_n("../ResultsNetP/final_plots_08_hadPID/out_sys_MB_18_08_finalBinning_singleParticleHighOrder_k3k1.root");

  TFile fin_k2k1_n_hm("../ResultsNetP/final_plots_08_hadPID/out_sys_HM_18_08_finalBinning_singleParticleHighOrder_k2k1.root");
  TFile fin_k3k1_n_hm("../ResultsNetP/final_plots_08_hadPID/out_sys_HM_18_08_finalBinning_singleParticleHighOrder_k3k1.root");

  TFile fin_k2k1_n_model("/home/mciacco/Code/NetProtonCumulants/out_sys_18_278_finalBinning_k2k1_n.root");
  TFile fin_k3k1_n_model("/home/mciacco/Code/NetProtonCumulants/out_sys_18_278_finalBinning_k3k1_n.root");
  TFile fin_k1_n_model("/home/mciacco/Code/NetProtonCumulants/out_sys_18_278_finalBinning_k1_n.root");

  TGraphErrors *g_k2k1_n_ = static_cast<TGraphErrors*>(fin_k2k1_n.Get("g_364"));
  TGraphErrors *g_k3k1_n_ = static_cast<TGraphErrors*>(fin_k3k1_n.Get("g_364"));

  TGraphErrors *g_k2k1_n_hm = static_cast<TGraphErrors*>(fin_k2k1_n_hm.Get("g_364"));
  TGraphErrors *g_k3k1_n_hm = static_cast<TGraphErrors*>(fin_k3k1_n_hm.Get("g_364"));

  TGraphErrors *g_k2k1_n_model = static_cast<TGraphErrors*>(fin_k2k1_n_model.Get("g_364"));
  TGraphErrors *g_k3k1_n_model = static_cast<TGraphErrors*>(fin_k3k1_n_model.Get("g_364"));
  TGraphErrors *g_k1_n_model = static_cast<TGraphErrors*>(fin_k1_n_model.Get("g_364"));

  TGraphErrors* g_k2k1_n = new TGraphErrors();
  TGraphErrors* g_k3k1_n = new TGraphErrors();
  TH1D* sys_k2k1_hm = (TH1D*)fin_k2k1_n_hm.Get("hSys_0");
  TH1D* sys_k3k1_hm = (TH1D*)fin_k3k1_n_hm.Get("hSys_0");
  g_k2k1_n->AddPoint(g_k2k1_n_hm->GetPointX(0), g_k2k1_n_hm->GetPointY(0));
  g_k2k1_n->SetPointError(g_k2k1_n->GetN() - 1, g_k2k1_n_hm->GetErrorX(0), std::hypot(g_k2k1_n_hm->GetErrorY(0), sys_k2k1_hm->GetStdDev()));
  g_k3k1_n->AddPoint(g_k3k1_n_hm->GetPointX(0), g_k3k1_n_hm->GetPointY(0));
  g_k3k1_n->SetPointError(g_k3k1_n->GetN() - 1, g_k3k1_n_hm->GetErrorX(0), std::hypot(g_k3k1_n_hm->GetErrorY(0), sys_k3k1_hm->GetStdDev()));
  for (int i{0}; i < 7; ++i) {
    TH1D* sys_k2k1 = (TH1D*)fin_k2k1_n.Get(Form("hSys_%i", i));
    TH1D* sys_k3k1 = (TH1D*)fin_k3k1_n.Get(Form("hSys_%i", i));
    g_k2k1_n->AddPoint(g_k2k1_n_->GetPointX(i), g_k2k1_n_->GetPointY(i));
    g_k2k1_n->SetPointError(g_k2k1_n->GetN() - 1, g_k2k1_n_->GetErrorX(i), std::hypot(g_k2k1_n_->GetErrorY(i), sys_k2k1->GetStdDev()));
    g_k3k1_n->AddPoint(g_k3k1_n_->GetPointX(i), g_k3k1_n_->GetPointY(i));
    g_k3k1_n->SetPointError(g_k3k1_n->GetN() - 1, g_k3k1_n_->GetErrorX(i), std::hypot(g_k3k1_n_->GetErrorY(i), sys_k3k1->GetStdDev()));
  }

  for (int i{0}; i < 8; ++i) {
    double k2k1_n = g_k2k1_n->GetPointY(i);
    double k1_n_model = g_k1_n_model->GetPointY(i);
    double k2k1_n_model = g_k2k1_n_model->GetPointY(i);
    double k2k1_n_err = g_k2k1_n->GetErrorY(i);
    double k1_n_model_err = g_k1_n_model->GetErrorY(i);
    double k2k1_n_model_err = g_k2k1_n_model->GetErrorY(i);
    std::cout << "v_2 = " << ( k2k1_n - k2k1_n_model ) / k1_n_model << std::endl;
  }
  for (int i{0}; i < 8; ++i) {
    double k2k1_n = g_k2k1_n->GetPointY(i);
    double k1_n_model = g_k1_n_model->GetPointY(i);
    double k2k1_n_model = g_k2k1_n_model->GetPointY(i);
    double k2k1_n_err = g_k2k1_n->GetErrorY(i);
    double k1_n_model_err = g_k1_n_model->GetErrorY(i);
    double k2k1_n_model_err = g_k2k1_n_model->GetErrorY(i);
    std::cout << "v2_err = " << std::sqrt( std::pow(k2k1_n_err / k1_n_model, 2.)/* + std::pow(k2k1_n_model_err / k1_n_model, 2.) + std::pow(k1_n_model_err * (k2k1_n - k2k1_n_model) / std::pow(k1_n_model, 2.), 2.)*/) << std::endl;
  }
  for (int i{0}; i < 8; ++i) {
    double k2k1_n = g_k2k1_n->GetPointY(i);
    double k3k1_n = g_k3k1_n->GetPointY(i);
    double k2k1_n_model = g_k2k1_n_model->GetPointY(i);
    double k1_n_model = g_k1_n_model->GetPointY(i);
    double k3k1_n_model = g_k3k1_n_model->GetPointY(i);
    std::cout << "v_3 = " << ( k3k1_n  - k3k1_n_model - 3 * k2k1_n_model * ( k2k1_n - k2k1_n_model ) ) / std::pow(k1_n_model, 2.) << std::endl;
  }
  for (int i{0}; i < 8; ++i) {
    double k2k1_n = g_k2k1_n->GetPointY(i);
    double k3k1_n = g_k3k1_n->GetPointY(i);
    double k2k1_n_model = g_k2k1_n_model->GetPointY(i);
    double k3k1_n_model = g_k3k1_n_model->GetPointY(i);
    double k1_n_model = g_k1_n_model->GetPointY(i);

    double k2k1_n_err = g_k2k1_n->GetErrorY(i);
    double k3k1_n_err = g_k3k1_n->GetErrorY(i);
    double k2k1_n_model_err = g_k2k1_n_model->GetErrorY(i);
    double k3k1_n_model_err = g_k3k1_n_model->GetErrorY(i);
    double k1_n_model_err = g_k1_n_model->GetErrorY(i);

    double d_k3k1 = k3k1_n_err / std::pow(k1_n_model, 2.);
    double d_k3k1_mod = k3k1_n_model_err / std::pow(k1_n_model, 2.);
    double d_k2k1_mod = 3 * k2k1_n_model * (k2k1_n - 2) * k2k1_n_model_err;
    double d_k2k1 = 3 * k2k1_n_model * k2k1_n_err;
    double d_k1_mod = 2 * (k3k1_n - k3k1_n_model - 3. * k2k1_n_model * (k2k1_n - k2k1_n_model)) / std::pow(k1_n_model, 3.) * k1_n_model_err;

    std::cout << "v3_err = " << std::sqrt( std::pow(d_k3k1, 2.) /*+ std::pow(d_k3k1_mod, 2.) + std::pow(d_k2k1_mod, 2.)*/ + std::pow(d_k2k1, 2.)/* + std::pow(d_k1_mod, 2.)*/) << std::endl;
  }

  f.cd();
  g_k2k1_n->Write();
  g_k2k1_n_model->Write();
  g_k1_n_model->Write();
  f.Close();
}

