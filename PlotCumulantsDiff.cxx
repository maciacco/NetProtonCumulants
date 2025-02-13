const char* obs_ax[]{"#kappa_{6}/#kappa_{2}", "#kappa_{4}/#kappa_{2}", "#kappa_{2}/#kappa_{2,Sk}", "#kappa_{2}^{#bar{p}}/#kappa_{1}^{#bar{p}}", "#kappa_{3}^{#bar{p}}/#kappa_{1}^{#bar{p}}"};
const char* obs_sr[]{"k6k2", "k4k2", "k2k2sk", "k2k1", "k3k1"};
const char* fname_ = "18";
const char* fname_m[]{"18_28_volF_23", "18_28_volF_23", "18_28", "18_28", "18_28"};
const char* fname_dat[]{"", "", "", "_singleParticleHighOrder", "_singleParticleHighOrder"};
const double yax_lim[][2] = {{0.24, 1.2}, {0.845, 1.03}, {0.9, .93}, {0.95, 1.02}, {0.86, 1.05}};
const double yax_lim_diff[][2] = {{-0.17, 0.17}, {-0.011, 0.011}, {-0.2, .2}, {0.95, 1.02}, {0.86, 1.05}};
const char* volf_str[]{" + vol. fluct.", " + vol. fluct.", "", "", ""};
//const double xerr[]{0.16, 0.12, 0.13, 0.11, 0.085, 0.064, 0.036};
//const double xerr_hm{0.26};

void SetGraphStyle(TGraph* g, Color_t const color = kRed){
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  g->SetMarkerColor(color);
  g->SetLineWidth(2);
  g->SetLineColor(color);
  g->SetFillStyle(0);
};

void SetGraphStyleModel(TGraph* g, Color_t const color = kBlue){
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.);
  g->SetMarkerColor(color);
  g->SetLineWidth(2);
  g->SetLineColor(color);
};

void PlotCumulantsDiff(const int obs = 0){
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile *f_mb = TFile::Open(Form("../ResultsNetP/final_plots/out_sys_MB_%s_08_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm = TFile::Open(Form("../ResultsNetP/final_plots/out_sys_HM_%s_08_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
/*  TFile *f_mb_6 = TFile::Open(Form("../ResultsNetP/final_plots_06/out_sys_MB_%s_06_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm_6 = TFile::Open(Form("../ResultsNetP/final_plots_06/out_sys_HM_%s_06_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_mb_4 = TFile::Open(Form("../ResultsNetP/final_plots_04/out_sys_MB_%s_04_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm_4 = TFile::Open(Form("../ResultsNetP/final_plots_04/out_sys_HM_%s_04_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
*/
  TFile *f_mb_m = TFile::Open(Form("out_sys_%s_finalBinning_%s.root", fname_m[obs], obs_sr[obs]));
//  TFile *f_hm_m = TFile::Open(Form("../test_models_thermalFistDowngrade/out_sys_HM_%s_finalBinning_%s_model.root", fname_m, obs_sr[obs]));
  TGraphErrors gSys_mb, gSys_hm;
  TCanvas c("c", "c", 600, 700);

  TPad p1("p1", "p1", 0., .35, 1., 1., 0);
  p1.SetTopMargin(0.03);
  p1.SetRightMargin(0.03);
  p1.SetLeftMargin(0.14);
  p1.SetBottomMargin(0.);
  c.cd();
  p1.Draw();

  TPad p2("p2", "p2", 0., 0., 1., .35, 0);
  p2.SetTopMargin(0.0);
  p2.SetRightMargin(0.03);
  p2.SetLeftMargin(0.14);
  p2.SetBottomMargin(0.32);
  c.cd();
  p2.Draw();

  TH2D hFrame("hFrame", Form(";#LTd#it{N}/d#eta#GT_{|#eta|<0.5};%s", obs_ax[obs]), 1, 0, 34., 100, yax_lim[obs][0], yax_lim[obs][1]);
  hFrame.GetXaxis()->SetTitleFont(45);
  hFrame.GetXaxis()->SetTitleSize(30);
  hFrame.GetXaxis()->SetTitleOffset(1.1);
  hFrame.GetYaxis()->SetTitleFont(45);
  hFrame.GetYaxis()->SetTitleSize(30);
  hFrame.GetYaxis()->SetTitleOffset(1.35);
  hFrame.GetXaxis()->SetLabelFont(45);
  hFrame.GetXaxis()->SetLabelSize(17);
  hFrame.GetYaxis()->SetLabelFont(45);
  hFrame.GetYaxis()->SetLabelSize(17);
  TGraphErrors* g_mb = (TGraphErrors*)f_mb->Get("g_364");
  TGraphErrors* g_hm = (TGraphErrors*)f_hm->Get("g_364");
/*  TGraphErrors* g_mb_6 = (TGraphErrors*)f_mb_6->Get("g_364");
  TGraphErrors* g_hm_6 = (TGraphErrors*)f_hm_6->Get("g_364");
  TGraphErrors* g_mb_4 = (TGraphErrors*)f_mb_4->Get("g_364");
  TGraphErrors* g_hm_4 = (TGraphErrors*)f_hm_4->Get("g_364");
*/
  TGraphErrors* g_mb_s = new TGraphErrors();
  TGraphErrors* g_hm_s = new TGraphErrors();
/*  TGraphErrors* g_mb_s_6 = new TGraphErrors();
  TGraphErrors* g_hm_s_6 = new TGraphErrors();
  TGraphErrors* g_mb_s_4 = new TGraphErrors();
  TGraphErrors* g_hm_s_4 = new TGraphErrors();
*/
  TGraphErrors* g_mb_m = (TGraphErrors*)f_mb_m->Get("g_364");
  g_mb_m->SetFillColor(kAzure + 1);
  g_mb_m->SetLineColor(kAzure + 1);
  g_mb_m->SetLineWidth(2);
  g_mb_m->SetFillStyle(3154);
  g_mb_m->SetMarkerStyle(20);
  g_mb_m->SetMarkerSize(0);
//  TGraphErrors* g_hm_m = (TGraphErrors*)f_hm_m->Get("g_gen_0");
  for (int i{0}; i < 7; ++i) {
    auto hSys = (TH1D*)f_mb->Get(Form("hSys_%d", i));
    g_mb_s->AddPoint(g_mb->GetPointX(i), g_mb->GetPointY(i));
    double min = hSys->GetBinCenter(hSys->FindFirstBinAbove(0));
    double max = hSys->GetBinCenter(hSys->FindLastBinAbove(0));
    g_mb_s->SetPointError(i, 0.5, /*xerr[i],*/ /*0.5 * (max - min)*/ hSys->GetStdDev());
  /*  hSys = (TH1D*)f_mb_6->Get(Form("hSys_%d", i));
    g_mb_s_6->AddPoint(g_mb_6->GetPointX(i), g_mb_6->GetPointY(i));
    g_mb_s_6->SetPointError(i, 0.5, hSys->GetStdDev());
    hSys = (TH1D*)f_mb_6->Get(Form("hSys_%d", i));
    g_mb_s_4->AddPoint(g_mb_4->GetPointX(i), g_mb_4->GetPointY(i));
    g_mb_s_4->SetPointError(i, 0.5, hSys->GetStdDev());
*/
  }
  for (int i{0}; i < 1; ++i) {
    auto hSys = (TH1D*)f_hm->Get(Form("hSys_%d", i));
    g_hm_s->AddPoint(g_hm->GetPointX(i), g_hm->GetPointY(i));
    g_hm_s->SetPointError(i, 0.5, /*xerr_hm,*/ hSys->GetStdDev());
  /*  hSys = (TH1D*)f_hm->Get(Form("hSys_%d", i));
    g_hm_s_6->AddPoint(g_hm_6->GetPointX(i), g_hm_6->GetPointY(i));
    g_hm_s_6->SetPointError(i, 0.5, hSys->GetStdDev());
    hSys = (TH1D*)f_hm->Get(Form("hSys_%d", i));
    g_hm_s_4->AddPoint(g_hm_4->GetPointX(i), g_hm_4->GetPointY(i));
    g_hm_s_4->SetPointError(i, 0.5, hSys->GetStdDev());
*/
  }
  SetGraphStyle(g_mb);
  SetGraphStyle(g_hm);
  SetGraphStyle(g_mb_s);
  SetGraphStyle(g_hm_s);
/*
  SetGraphStyle(g_mb_6, kBlue);
  SetGraphStyle(g_hm_6, kBlue);
  SetGraphStyle(g_mb_s_6, kBlue);
  SetGraphStyle(g_hm_s_6, kBlue);

  SetGraphStyle(g_mb_4, kGreen + 2);
  SetGraphStyle(g_hm_4, kGreen + 2);
  SetGraphStyle(g_mb_s_4, kGreen + 2);
  SetGraphStyle(g_hm_s_4, kGreen + 2);
*/
//  SetGraphStyleModel(g_mb_m);
//  SetGraphStyleModel(g_hm_m);
  p1.cd();
  hFrame.Draw();
  g_mb_m->Draw("samee3l");
  g_mb_s->Draw("samee5");
  g_hm_s->Draw("samee5");
  g_mb->Draw("samepez");
  g_hm->Draw("samepez");

//  g_mb_s_6->Draw("samee5");
//  g_hm_s_6->Draw("samee5");
//  g_mb_6->Draw("samepe");
//  g_hm_6->Draw("samepe");

//  g_mb_s_4->Draw("samee5");
//  g_hm_s_4->Draw("samee5");
//  g_mb_4->Draw("samepe");
//  g_hm_4->Draw("samepe");

//  g_hm_m->Draw("samele");
  TGraph g_ph;
  g_ph.SetLineWidth(0);
  g_ph.SetMarkerStyle(20);
  g_ph.SetMarkerSize(0);
  TLegend leg(0.18, 0.58, 0.4, 0.72);
  leg.SetTextFont(45);
  leg.SetTextSize(17);
  leg.SetBorderSize(0);
  leg.AddEntry(g_mb, "Data", "pe");
  leg.AddEntry(g_mb_m, Form("Thermal-FIST ev. gen. + BW%s, #it{V}_{c} = 2.8 d#it{V}/d#it{y}", volf_str[obs]), "f");
  leg.AddEntry(&g_ph, "#it{T}_{chem}, d#it{V}/d#it{y}, and #gamma_{s} from Phys. Rev. C 100 (2019) 054906");
  leg.Draw("same");

  TLatex txt;
  txt.SetNDC();
  txt.SetTextFont(45);
  txt.SetTextSize(25);
  txt.DrawLatex(0.18, 0.9, "ALICE Preliminary");
  txt.DrawLatex(0.18, 0.825, "pp, #sqrt{s} = 13 TeV");
  txt.DrawLatex(0.18, 0.75, "|#eta| < 0.8, 0.5 < #it{p}_{T} < 1.5 GeV/#it{c}");

  p2.cd();
  TH2D hFrame2("hFrame2", ";#LTd#it{N}/d#eta#GT_{|#eta|<0.5};Data - Model", 1, 0, 34., 100, yax_lim_diff[obs][0], yax_lim_diff[obs][1]);
  hFrame2.GetXaxis()->SetTitleFont(45);
  hFrame2.GetXaxis()->SetTitleSize(30);
  hFrame2.GetXaxis()->SetTitleOffset(1.);
  hFrame2.GetYaxis()->SetTitleFont(45);
  hFrame2.GetYaxis()->SetTitleSize(30);
  hFrame2.GetYaxis()->SetTitleOffset(1.35);
  hFrame2.GetXaxis()->SetLabelFont(45);
  hFrame2.GetXaxis()->SetLabelSize(17);
  hFrame2.GetYaxis()->SetLabelFont(45);
  hFrame2.GetYaxis()->SetLabelSize(17);
  hFrame2.GetYaxis()->SetNdivisions(405);

  hFrame2.Draw();
  TGraphErrors* diff_mb = new TGraphErrors(*g_mb);
  TGraphErrors* diff_hm = new TGraphErrors(*g_hm);
  TGraphErrors* diff_mb_s = new TGraphErrors(*g_mb_s);
  TGraphErrors* diff_hm_s = new TGraphErrors(*g_hm_s);
  TGraphErrors* diff_mb_m = new TGraphErrors(*g_mb_m);
  for (int i{0}; i < 7; ++i) {
    diff_mb->SetPointY(i, diff_mb->GetPointY(i) - g_mb_m->GetPointY(i + 1));
    diff_mb_s->SetPointY(i, diff_mb_s->GetPointY(i) - g_mb_m->GetPointY(i + 1));
    diff_mb_m->SetPointY(i + 1, 0.);
  }
  diff_hm->SetPointY(0, diff_hm->GetPointY(0) - g_mb_m->GetPointY(0));
  diff_hm_s->SetPointY(0, diff_hm_s->GetPointY(0) - g_mb_m->GetPointY(0));
  diff_mb_m->SetPointY(0, 0.);

  diff_mb_m->Draw("samee3l");
  diff_mb_s->Draw("samee5");
  diff_mb->Draw("samepez");
  diff_hm_s->Draw("samee5");
  diff_hm->Draw("samepez");

  c.Print(Form("c%s.pdf", obs_sr[obs]));

  TFile *fout = TFile::Open("plot_out.root", "recreate");
  fout->cd();
  c.Write();
  fout->Close();
}
