const char* obs_ax[]{"#kappa_{6}/#kappa_{2}", "#kappa_{4}/#kappa_{2}", "#kappa_{2}/#kappa_{2,Sk}", "#kappa_{2}^{#bar{p}}/#kappa_{1}^{#bar{p}}", "#kappa_{3}^{#bar{p}}/#kappa_{1}^{#bar{p}}"};
const char* obs_sr[]{"k6k2", "k4k2", "k2k2sk", "k2k1", "k3k1"};
const char* fname_ = "18";
const char* fname_m[]{"18_28_volF_23", "18_28_volF_23", "18_28", "18_28", "18_28"};
const char* fname_dat[]{"", "", "", "_singleParticleHighOrder", "_singleParticleHighOrder"};
const double yax_lim[][2] = {{0.25, 1.2}, {0.84, 1.03}, {0.89, 1.015}, {0.95, 1.02}, {0.86, 1.05}};
const char* volf_str[]{" + vol. fluct.", " + vol. fluct.", "", "", ""};

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

void PlotCumulantsVsAcc(const int obs = 2){
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile *f_mb = TFile::Open(Form("../ResultsNetP/final_plots/out_sys_MB_%s_08_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm = TFile::Open(Form("../ResultsNetP/final_plots/out_sys_HM_%s_08_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_mb_6 = TFile::Open(Form("../ResultsNetP/final_plots_06/out_sys_MB_%s_06_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm_6 = TFile::Open(Form("../ResultsNetP/final_plots_06/out_sys_HM_%s_06_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_mb_4 = TFile::Open(Form("../ResultsNetP/final_plots_04/out_sys_MB_%s_04_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm_4 = TFile::Open(Form("../ResultsNetP/final_plots_04/out_sys_HM_%s_04_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_mb_2 = TFile::Open(Form("../ResultsNetP/final_plots_02/out_sys_MB_%s_02_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm_2 = TFile::Open(Form("../ResultsNetP/final_plots_02/out_sys_HM_%s_02_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));

//  TFile *f_mb_m = TFile::Open(Form("out_sys_%s_finalBinning_%s.root", fname_m[obs], obs_sr[obs]));
//  TFile *f_hm_m = TFile::Open(Form("../test_models_thermalFistDowngrade/out_sys_HM_%s_finalBinning_%s_model.root", fname_m, obs_sr[obs]));
  TGraphErrors gSys_mb, gSys_hm;
  TCanvas c("c", "c", 600, 600);
  c.SetTopMargin(0.03);
  c.SetRightMargin(0.03);
  c.SetLeftMargin(0.14);
  c.SetBottomMargin(0.14);
  TH2D hFrame("hFrame", Form(";#LTd#it{N}/d#eta#GT_{|#eta|<0.5};%s", obs_ax[obs]), 1, 0, 34., 100, yax_lim[obs][0], yax_lim[obs][1]);
  hFrame.GetXaxis()->SetTitleFont(45);
  hFrame.GetXaxis()->SetTitleSize(30);
  hFrame.GetXaxis()->SetTitleOffset(1.1);
  hFrame.GetYaxis()->SetTitleFont(45);
  hFrame.GetYaxis()->SetTitleSize(30);
  hFrame.GetYaxis()->SetTitleOffset(1.3);
  TGraphErrors* g_mb = (TGraphErrors*)f_mb->Get("g_364");
  TGraphErrors* g_hm = (TGraphErrors*)f_hm->Get("g_364");
  TGraphErrors* g_mb_6 = (TGraphErrors*)f_mb_6->Get("g_364");
  TGraphErrors* g_hm_6 = (TGraphErrors*)f_hm_6->Get("g_364");
  TGraphErrors* g_mb_4 = (TGraphErrors*)f_mb_4->Get("g_364");
  TGraphErrors* g_hm_4 = (TGraphErrors*)f_hm_4->Get("g_364");
  TGraphErrors* g_mb_2 = (TGraphErrors*)f_mb_2->Get("g_364");
  TGraphErrors* g_hm_2 = (TGraphErrors*)f_hm_2->Get("g_364");

  TGraphErrors* g_mb_s = new TGraphErrors();
  TGraphErrors* g_hm_s = new TGraphErrors();
  TGraphErrors* g_mb_s_6 = new TGraphErrors();
  TGraphErrors* g_hm_s_6 = new TGraphErrors();
  TGraphErrors* g_mb_s_4 = new TGraphErrors();
  TGraphErrors* g_hm_s_4 = new TGraphErrors();
  TGraphErrors* g_mb_s_2 = new TGraphErrors();
  TGraphErrors* g_hm_s_2 = new TGraphErrors();

//  TGraphErrors* g_mb_m = (TGraphErrors*)f_mb_m->Get("g_364");
//  g_mb_m->SetFillColor(kAzure + 1);
//  g_mb_m->SetLineColor(kAzure + 1);
//  g_mb_m->SetLineWidth(2);
//  g_mb_m->SetFillStyle(3154);
//  g_mb_m->SetMarkerStyle(20);
//  g_mb_m->SetMarkerSize(0);
//  TGraphErrors* g_hm_m = (TGraphErrors*)f_hm_m->Get("g_gen_0");
  for (int i{0}; i < 7; ++i) {
    auto hSys = (TH1D*)f_mb->Get(Form("hSys_%d", i));
    g_mb_s->AddPoint(g_mb->GetPointX(i), g_mb->GetPointY(i));
    double min = hSys->GetBinCenter(hSys->FindFirstBinAbove(0));
    double max = hSys->GetBinCenter(hSys->FindLastBinAbove(0));
    g_mb_s->SetPointError(i, 0.5, /*0.5 * (max - min)*/ hSys->GetStdDev());
    hSys = (TH1D*)f_mb_6->Get(Form("hSys_%d", i));
    g_mb_s_6->AddPoint(g_mb_6->GetPointX(i), g_mb_6->GetPointY(i));
    g_mb_s_6->SetPointError(i, 0.5, hSys->GetStdDev());
    hSys = (TH1D*)f_mb_4->Get(Form("hSys_%d", i));
    g_mb_s_4->AddPoint(g_mb_4->GetPointX(i), g_mb_4->GetPointY(i));
    g_mb_s_4->SetPointError(i, 0.5, hSys->GetStdDev());
    hSys = (TH1D*)f_mb_2->Get(Form("hSys_%d", i));
    g_mb_s_2->AddPoint(g_mb_2->GetPointX(i), g_mb_2->GetPointY(i));
    g_mb_s_2->SetPointError(i, 0.5, hSys->GetStdDev());
  }
  for (int i{0}; i < 1; ++i) {
    auto hSys = (TH1D*)f_hm->Get(Form("hSys_%d", i));
    g_hm_s->AddPoint(g_hm->GetPointX(i), g_hm->GetPointY(i));
    g_hm_s->SetPointError(i, 0.5, hSys->GetStdDev());
    hSys = (TH1D*)f_hm_6->Get(Form("hSys_%d", i));
    g_hm_s_6->AddPoint(g_hm_6->GetPointX(i), g_hm_6->GetPointY(i));
    g_hm_s_6->SetPointError(i, 0.5, hSys->GetStdDev());
    hSys = (TH1D*)f_hm_4->Get(Form("hSys_%d", i));
    g_hm_s_4->AddPoint(g_hm_4->GetPointX(i), g_hm_4->GetPointY(i));
    g_hm_s_4->SetPointError(i, 0.5, hSys->GetStdDev());
    hSys = (TH1D*)f_hm_2->Get(Form("hSys_%d", i));
    g_hm_s_2->AddPoint(g_hm_2->GetPointX(i), g_hm_2->GetPointY(i));
    g_hm_s_2->SetPointError(i, 0.5, hSys->GetStdDev());
  }
  SetGraphStyle(g_mb);
  SetGraphStyle(g_hm);
  SetGraphStyle(g_mb_s);
  SetGraphStyle(g_hm_s);

  SetGraphStyle(g_mb_6, kOrange);
  SetGraphStyle(g_hm_6, kOrange);
  SetGraphStyle(g_mb_s_6, kOrange);
  SetGraphStyle(g_hm_s_6, kOrange);

  SetGraphStyle(g_mb_4, kGreen + 2);
  SetGraphStyle(g_hm_4, kGreen + 2);
  SetGraphStyle(g_mb_s_4, kGreen + 2);
  SetGraphStyle(g_hm_s_4, kGreen + 2);

  SetGraphStyle(g_mb_2, kBlue);
  SetGraphStyle(g_hm_2, kBlue);
  SetGraphStyle(g_mb_s_2, kBlue);
  SetGraphStyle(g_hm_s_2, kBlue);

//  SetGraphStyleModel(g_mb_m);
//  SetGraphStyleModel(g_hm_m);
  c.cd();
  hFrame.Draw();
//  g_mb_m->Draw("samee3l");
  g_mb_s->Draw("samee5");
  g_hm_s->Draw("samee5");
  g_mb->Draw("samepez");
  g_hm->Draw("samepez");

  g_mb_s_6->Draw("samee5");
  g_hm_s_6->Draw("samee5");
  g_mb_6->Draw("samepez");
  g_hm_6->Draw("samepez");

  g_mb_s_4->Draw("samee5");
  g_hm_s_4->Draw("samee5");
  g_mb_4->Draw("samepez");
  g_hm_4->Draw("samepez");

  g_mb_s_2->Draw("samee5");
  g_hm_s_2->Draw("samee5");
  g_mb_2->Draw("samepez");
  g_hm_2->Draw("samepez");

//  g_hm_m->Draw("samele");
  TGraph g_ph;
  g_ph.SetLineWidth(0);
  g_ph.SetMarkerStyle(20);
  g_ph.SetMarkerSize(0);
  TLegend leg(0.7, 0.75, 0.8, 0.95);
  leg.SetTextFont(45);
  leg.SetTextSize(25);
  leg.SetBorderSize(0);
//  leg.SetNColumns(2);
  leg.AddEntry(g_mb_2, "|#eta| < 0.2", "pe");
  leg.AddEntry(g_mb_4, "|#eta| < 0.4", "pe");
  leg.AddEntry(g_mb_6, "|#eta| < 0.6", "pe");
  leg.AddEntry(g_mb, "|#eta| < 0.8", "pe");

//  leg.AddEntry(g_mb_m, Form("Thermal-FIST ev. gen. + BW%s, #it{V}_{c} = 2.8 d#it{V}/d#it{y}", volf_str[obs]), "f");
//  leg.AddEntry(&g_ph, "#it{T}_{chem}, d#it{V}/d#it{y}, and #gamma_{s} from Phys. Rev. C 100 (2019) 054906");
  leg.Draw("same");

  TLatex txt;
  txt.SetNDC();
  txt.SetTextFont(45);
  txt.SetTextSize(25);
  txt.DrawLatex(0.18, 0.9, "ALICE Preliminary");
  txt.DrawLatex(0.18, 0.85, "pp, #sqrt{s} = 13 TeV");
  txt.DrawLatex(0.18, 0.8, "0.5 < #it{p}_{T} < 1.5 GeV/#it{c}");
  c.Print(Form("c%s_deta.pdf", obs_sr[obs]));

  TFile *fout = TFile::Open("plot_out.root", "recreate");
  fout->cd();
  c.Write();
  fout->Close();
}
