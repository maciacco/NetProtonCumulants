const char* obs_ax[]{"#kappa_{6}/#kappa_{2}", "#kappa_{4}/#kappa_{2}", "#kappa_{2}/#kappa_{2,Sk}"};
const char* obs_sr[]{"k6k2", "k4k2", "k2k2sk"};
const char* fname_ = "18_08";
const char* fname_m = "18_08_VF";
const double yax_lim[][2] = {{0.1, 1.1}, {0.8, 1.}, {0.89, 1.}};

void SetGraphStyle(TGraph* g, Color_t const color = kRed){
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  g->SetMarkerColor(color);
  g->SetLineWidth(2);
  g->SetLineColor(color);
};

void SetGraphStyleModel(TGraph* g, Color_t const color = kBlue){
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.);
  g->SetMarkerColor(color);
  g->SetLineWidth(2);
  g->SetLineColor(color);
};

void PlotCumulantsVsAcc(const int obs = 2){
  gStyle->SetOptStat(0);
  TFile *f_mb = TFile::Open(Form("../ResultsNetP/final_plots/out_sys_MB_%s_finalBinning_%s.root", fname_, obs_sr[obs]));
  TFile *f_hm = TFile::Open(Form("../ResultsNetP/final_plots/out_sys_HM_%s_finalBinning_%s.root", fname_, obs_sr[obs]));
//  TFile *f_mb_m = TFile::Open(Form("../test_models_thermalFistDowngrade/out_sys_MB_%s_finalBinning_%s_model.root", fname_m, obs_sr[obs]));
//  TFile *f_hm_m = TFile::Open(Form("../test_models_thermalFistDowngrade/out_sys_HM_%s_finalBinning_%s_model.root", fname_m, obs_sr[obs]));
  TGraphErrors gSys_mb, gSys_hm;
  TCanvas c("c", "c", 600, 600);
  c.SetTopMargin(0.03);
  c.SetRightMargin(0.03);
  c.SetLeftMargin(0.14);
  c.SetBottomMargin(0.16);
  TH2D hFrame("hFrame", Form(";#Delta#eta;%s", obs_ax[obs]), 1, 0, 1.8, 1, yax_lim[obs][0], yax_lim[obs][1]);
  hFrame.GetXaxis()->SetTitleFont(45);
  hFrame.GetXaxis()->SetTitleSize(25);
  hFrame.GetXaxis()->SetTitleOffset(1.5);
  hFrame.GetYaxis()->SetTitleFont(45);
  hFrame.GetYaxis()->SetTitleSize(25);
  TGraphErrors* g_mb = (TGraphErrors*)f_mb->Get("g_364");
  TGraphErrors* g_hm = (TGraphErrors*)f_hm->Get("g_364");
  TGraphErrors* g_mb_s = new TGraphErrors();
  TGraphErrors* g_hm_s = new TGraphErrors();
//  TGraphErrors* g_mb_m = (TGraphErrors*)f_mb_m->Get("g_gen_0");
//  TGraphErrors* g_hm_m = (TGraphErrors*)f_hm_m->Get("g_gen_0");
  for (int i{0}; i < 7; ++i) {
    auto hSys = (TH1D*)f_mb->Get(Form("hSys_%d", i));
    g_mb_s->AddPoint(g_mb->GetPointX(i), g_mb->GetPointY(i));
    g_mb_s->SetPointError(i, 0.8, hSys->GetStdDev());
  }
  for (int i{0}; i < 1; ++i) {
    auto hSys = (TH1D*)f_hm->Get(Form("hSys_%d", i));
    g_hm_s->AddPoint(g_hm->GetPointX(i), g_hm->GetPointY(i));
    g_hm_s->SetPointError(i, 0.8, hSys->GetStdDev());
  }
  SetGraphStyle(g_mb);
  SetGraphStyle(g_hm);
  SetGraphStyle(g_mb_s);
  SetGraphStyle(g_hm_s);
//  SetGraphStyleModel(g_mb_m);
//  SetGraphStyleModel(g_hm_m);
  c.cd();
  hFrame.Draw();
  g_mb_s->Draw("samee5");
  g_hm_s->Draw("samee5");
  g_mb->Draw("samepe");
  g_hm->Draw("samepe");
//  g_mb_m->Draw("samele");
//  g_hm_m->Draw("samele");
  TLatex txt;
  txt.SetNDC();
  txt.SetTextFont(45);
  txt.SetTextSize(20);
  txt.DrawLatex(0.5, 0.9, "ALICE Preliminary");
  txt.DrawLatex(0.5, 0.85, "pp, #sqrt{s} = 13 TeV");
  txt.DrawLatex(0.5, 0.8, "|#eta| < 0.8, 0.5 < #it{p}_{T} < 1.5 GeV/#it{c}");
  c.Print(Form("c%s.pdf", obs_sr[obs]));
  TFile *fout = TFile::Open("plot_out.root", "recreate");
  fout->cd();
  c.Write();
  fout->Close();
}
