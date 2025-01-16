void computeVolumeFluct(){
  TFile fin_k2k1_n("out_sys_HM_18_finalBinning_singleParticleHighOrder_k2k1.root");
  TFile fin_k1_n("out_sys_HM_18_finalBinning_singleParticleHighOrder_k1.root");
  TFile fin_k2_n("out_sys_HM_18_finalBinning_singleParticleHighOrder_k2.root");
  TFile fin_k3_n("out_sys_HM_18_finalBinning_singleParticleHighOrder_k3.root");
  TFile fin_k2k1_n_model("out_sys_HM_18_finalBinning_n_k2k1.root");
  TFile fin_k1_n_model("out_sys_HM_18_finalBinning_n_k1.root");
  TFile fin_k2_n_model("out_sys_HM_18_finalBinning_n_k2.root");
  TFile fin_k3_n_model("out_sys_HM_18_finalBinning_n_k3.root");

  TGraphErrors *g_k2k1_n = static_cast<TGraphErrors*>(fin_k2k1_n.Get("g_364"));
  TGraphErrors *g_k1_n = static_cast<TGraphErrors*>(fin_k1_n.Get("g_364"));
  TGraphErrors *g_k2_n = static_cast<TGraphErrors*>(fin_k2_n.Get("g_364"));
  TGraphErrors *g_k3_n = static_cast<TGraphErrors*>(fin_k3_n.Get("g_364"));
  TGraphErrors *g_k2k1_n_model = static_cast<TGraphErrors*>(fin_k2k1_n_model.Get("g_gen_0"));
  TGraphErrors *g_k1_n_model = static_cast<TGraphErrors*>(fin_k1_n_model.Get("g_gen_0"));
  TGraphErrors *g_k2_n_model = static_cast<TGraphErrors*>(fin_k2_n_model.Get("g_gen_0"));
  TGraphErrors *g_k3_n_model = static_cast<TGraphErrors*>(fin_k3_n_model.Get("g_gen_0"));

  for (int i{0}; i < 1; ++i) {
    double k2k1_n = g_k2k1_n->GetPointY(i);
    double k1_n = g_k1_n->GetPointY(i);
    double k2k1_n_model = g_k2k1_n_model->GetPointY(i);
    double k2k1_n_err = g_k2k1_n->GetErrorY(i);
    double k1_n_err = g_k1_n->GetErrorY(i);
    double k2k1_n_model_err = g_k2k1_n_model->GetErrorY(i);
    //std::cout << k2k1_n << "\t" << k2k1_n_model << "\t" << k1_n << std::endl;
    std::cout << "v_2 = " << ( k2k1_n - k2k1_n_model ) / k1_n << std::endl;
    //std::cout << "v2_err = " << std::sqrt( std::pow(k2k1_n_err / k2k1_n, 2.) + std::pow(k2k1_n_model_err / k2k1_n_model, 2.) + std::pow((k2k1_n - k2k1_n_model) / k1_n / k1_n, 2.)) << std::endl;
  }
  for (int i{0}; i < 1; ++i) {
    double k2k1_n = g_k2k1_n->GetPointY(i);
    double k1_n = g_k1_n->GetPointY(i);
    double k2k1_n_model = g_k2k1_n_model->GetPointY(i);
    double k2k1_n_err = g_k2k1_n->GetErrorY(i);
    double k1_n_err = g_k1_n->GetErrorY(i);
    double k2k1_n_model_err = g_k2k1_n_model->GetErrorY(i);
    std::cout << k2k1_n_model_err << "\t" << k1_n_err << "\t" << k2k1_n_err << std::endl;
    //std::cout << "v_2 = " << ( k2k1_n - k2k1_n_model ) / k1_n << std::endl;
    std::cout << "v2_err = " << std::sqrt( std::pow(k2k1_n_err / k2k1_n, 2.) + std::pow(k2k1_n_model_err / k2k1_n_model, 2.) + std::pow(k1_n_err * (k2k1_n - k2k1_n_model) / k1_n / k1_n, 2.)) << std::endl;
  }
  for (int i{0}; i < 7; ++i) {
    double k2k1_n = g_k2k1_n->GetPointY(i);
    double k1_n = g_k1_n->GetPointY(i);
    double k2_n = g_k2_n->GetPointY(i);
    double k3_n = g_k3_n->GetPointY(i);
    double k2k1_n_model = g_k2k1_n_model->GetPointY(i);
    double k1_n_model = g_k1_n_model->GetPointY(i);
    double k2_n_model = g_k2_n_model->GetPointY(i);
    double k3_n_model = g_k3_n_model->GetPointY(i);
    std::cout << "v_3 = " << ( k3_n / k1_n  - k3_n_model / k1_n_model - 3 * k2_n * ( k2k1_n - k2k1_n_model ) / k1_n ) / std::pow(k1_n, 2.) << std::endl;
    //std::cout << "v2_err = " << std::sqrt( std::pow(k2k1_n_err / k2k1_n, 2.) + std::pow(k2k1_n_model_err / k2k1_n_model, 2.) + std::pow((k2k1_n - k2k1_n_model) / k1_n / k1_n, 2.)) << std::endl;
  }
}
