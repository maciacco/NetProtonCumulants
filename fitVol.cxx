const char* modName[]{/*"26_", */"27_", "28_", "29_", "30_"};
const Color_t colors[]{kOrange + 7, kOrange, kGreen + 3, kBlue, kViolet};
const double vol[]{/*2.6, */2.7, 2.8, 2.9, 3.0};
const double fac = 0;

Double_t chi2(TGraphErrors *gstat, /*TGraphErrors *gsyst,*/ TGraphErrors *gmodel){
  double res = 0.;
  for (int i{0}; i < gstat->GetN(); ++i){
    if (gstat->GetPointX(i) < 6.) continue;
    double o = gstat->GetPointY(i);
    double ostat = gstat->GetErrorY(i);
    double osyst = 0.; //gsyst->GetErrorY(i);
    double model = gmodel->GetPointY(i);
    double modelstat = gmodel->GetErrorY(i);

    double tmp = 0.;
    tmp = ( o + fac * osyst - model ) * ( o + fac * osyst - model );
    tmp = tmp / ( std::pow(ostat / o * (o + fac * osyst), 2) /* + osyst * osyst *//*  + modelstat * modelstat */);
    res += tmp;
  }
  return res;
}

Double_t chi2err(TGraphErrors *gstat, /*TGraphErrors *gsyst,*/ TGraphErrors *gmodel){
  double res = 0.;
  for (int i{0}; i < gstat->GetN(); ++i){
    if (gstat->GetPointX(i) < 6.) continue;
    double o = gstat->GetPointY(i);
    double ostat = gstat->GetErrorY(i);
    double osyst = 0.; //gsyst->GetErrorY(i);
    double model = gmodel->GetPointY(i);
    double modelstat = gmodel->GetErrorY(i);

    double tmp = 0.;
    //tmp = 2 * (o - model) * sqrt( ostat * ostat + osyst * osyst);
    tmp = tmp + 2 * (o + fac * osyst - model) * sqrt( modelstat * modelstat );
    tmp = tmp / ( std::pow(ostat / o * (o + fac * osyst), 2 ) /* + osyst * osyst *//*  + modelstat * modelstat  */);
    res += ( tmp * tmp);
  }
  return sqrt(res);
}

void fitVol(){
  TFile *fo = TFile::Open("chi2Fit_hm.root", "recreate");
  TCanvas c("c", "c", 600, 600);
  TCanvas cChi2("c", "c", 600, 600);
  c.SetLeftMargin(0.15);
  c.SetTopMargin(0.03);
  c.SetRightMargin(0.03);
  cChi2.SetLeftMargin(0.15);
  cChi2.SetTopMargin(0.03);
  cChi2.SetRightMargin(0.03);

  TFile *data = TFile::Open("../ResultsNetP/final_plots_08_hadPID/out_sys_MB_18_08_finalBinning_k2k2sk.root");
  TGraphErrors* gr_dat_ = (TGraphErrors*)data->Get("g_364");
  TFile *data_hm = TFile::Open("../ResultsNetP/final_plots_08_hadPID/out_sys_HM_18_08_finalBinning_k2k2sk.root");
  TGraphErrors* gr_dat_hm = (TGraphErrors*)data_hm->Get("g_364");
  TGraphErrors* gr_dat = new TGraphErrors();
  gr_dat->AddPoint(gr_dat_hm->GetPointX(0), gr_dat_hm->GetPointY(0));
  gr_dat->SetPointError(0, 0., gr_dat_hm->GetErrorY(0));
  for (int i{0}; i < 7; ++i) {
    gr_dat->AddPoint(gr_dat_->GetPointX(i), gr_dat_->GetPointY(i));
    gr_dat->SetPointError(gr_dat->GetN() - 1, 0., gr_dat_->GetErrorY(i));
  }
  gr_dat->SetMarkerColor(kRed);
  gr_dat->SetLineColor(kRed);

  std::cout << gr_dat_hm->GetErrorY(0) << std::endl;
  c.cd();
  gr_dat->SetMarkerStyle(20);
  gr_dat->SetMarkerSize(1.2);
  gr_dat->GetYaxis()->SetRangeUser(0.89, 0.93);
  gr_dat->GetXaxis()->SetTitle("#LTd#it{N}/d#eta#GT");
  gr_dat->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{2,Sk}");
  gr_dat->Draw("ape");
  TFile *models[5];
  TGraphErrors* gr[5];
  TGraphErrors gr_chi2;
  TF1 fchi2("fchi2", "pol2(0)+pol4(3)", 0, 4);
  for (int iM{0}; iM < 4; ++iM) {
    models[iM] = TFile::Open(Form("/home/mciacco/Code/out_sys_18_%snewBW__finalBinning_k2k2sk.root", modName[iM]));
    gr[iM] = (TGraphErrors*)models[iM]->Get("g_364");
    gr[iM]->SetLineColor(colors[iM]);
    c.cd();
    gr[iM]->Draw("lesame");
    double chi2_d = chi2(gr_dat, gr[iM]);
    double chi2_d_err = chi2err(gr_dat, gr[iM]);
    gr_chi2.AddPoint(vol[iM], chi2_d);
    gr_chi2.SetPointError(iM, 0., chi2_d_err);
  }
  gr_chi2.Fit("fchi2");
  gr_chi2.GetXaxis()->SetTitle("#it{V}_{c} (dV/dy)");
  gr_chi2.GetYaxis()->SetTitle("#chi^{2}");
  gr_chi2.GetYaxis()->SetRangeUser(0., 5800.);
  cChi2.cd();
  gr_chi2.Draw();
  std::cout << fchi2.GetMinimumX(2.5, 3.2) << std::endl;
  fo->cd();
  gr_chi2.Write();
  c.Write();
  c.Print("cvol_hm.pdf");
  cChi2.Print("cchi2_hm.pdf");
  fo->Close();
}
