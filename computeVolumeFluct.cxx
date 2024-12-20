void computeVolumeFluct(){
  TFile fin_k2k1_n("out_sys_finalBinning_singleParticle.root");
  TFile fin_k1_n("out_sys_finalBinning_singleParticle_avg.root");
  TFile fin_k2k1_n_model("out_sys_18_finalBinning_n.root");

  TGraphErrors *g_k2k1_n = static_cast<TGraphErrors*>(fin_k2k1_n.Get("g_364"));
  TGraphErrors *g_k1_n = static_cast<TGraphErrors*>(fin_k1_n.Get("g_364"));
  TGraphErrors *g_k2k1_n_model = static_cast<TGraphErrors*>(fin_k2k1_n_model.Get("g_gen_0"));

  for (int i{0}; i < 8; ++i) {
    double k2k1_n = g_k2k1_n->GetPointY(i);
    double k1_n = g_k1_n->GetPointY(i);
    double k2k1_n_model = g_k2k1_n_model->GetPointY(i);
    std::cout << "v_2 = " << ( k2k1_n - k2k1_n_model ) / k1_n << std::endl;
  }
}