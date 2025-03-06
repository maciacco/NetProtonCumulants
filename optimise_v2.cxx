void optimise_v2() {
  TFile *fdat_mb = TFile::Open("../ResultsNetP/final_plots_08_hadPID/out_sys_MB_18_08_finalBinning_singleParticleHighOrder_k2k1.root");
  TFile *fdat_hm = TFile::Open("../ResultsNetP/final_plots_08_hadPID/out_sys_HM_18_08_finalBinning_singleParticleHighOrder_k2k1.root");

  TGraphErrors* g_mb = static_cast<TGraphErrors*>(fdat_mb->Get("g_364"));
  TGraphErrors* g_hm = static_cast<TGraphErrors*>(fdat_hm->Get("g_364"));

  TFile *fvar[20];
  TGraphErrors *gvar[20];
  for (int iV{0}; iV < 20; ++iV) {
    fvar[iV] = TFile::Open(Form("/home/mciacco/out_sys_test_%d_finalBinning_k2k1_n.root", iV));
    gvar[iV] = (TGraphErrors*)fvar[iV]->Get("g_364");
  }

  TH2D hChi2VsMult("hChi2VsMult", ";Mult;#tilde{v}_{2}", 8, 0, 8, 20, -0.005, 0.195);
  for (int im{0}; im < 7; ++im) {
    double v_dat = g_mb->GetPointY(im);
    double e_dat = g_mb->GetErrorY(im);
    for (int iv{0}; iv < 20; ++iv) {
      double v_var = gvar[iv]->GetPointY(im + 1);
      hChi2VsMult.SetBinContent(im + 2, iv + 1, std::pow( (v_dat - v_var) / e_dat, 2. ));
    }
  }

  double v_dat = g_hm->GetPointY(0);
  double e_dat = g_hm->GetErrorY(0);
  for (int iv{0}; iv < 20; ++iv) {
    double v_var = gvar[iv]->GetPointY(0);
    hChi2VsMult.SetBinContent(1, iv + 1, std::pow( (v_dat - v_var) / e_dat, 2. ));
  }

  TFile *fo = TFile::Open("tuning.root", "recreate");
  fo->cd();
  hChi2VsMult.Write();

  for (int im{0}; im < 8; ++im) {
    TH1D* hp = static_cast<TH1D*>(hChi2VsMult.ProjectionY(Form("hChi2_%d", im), im + 1, im + 1));
    double xmin = hp->GetBinCenter(hp->GetMinimumBin());
    hp->Fit("pol2", "QRM+", "", xmin - 0.04, xmin + 0.04);
    auto fitf = hp->GetFunction("pol2");
    double min = -0.5 * fitf->GetParameter(1) / fitf->GetParameter(2);
    std::cout << "v2 = " << min << std::endl;
    fo->cd();

    hp->Write();
  }

  fo->Close();
}


