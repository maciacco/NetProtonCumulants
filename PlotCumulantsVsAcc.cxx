int colors[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"),
TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#595959"), TColor::GetColor("#0033a1")};

const char* obs_ax[]{"#kappa_{6}/#kappa_{2}", "#kappa_{4}/#kappa_{2}", "#kappa_{2}/#kappa_{2,Sk}", "#kappa_{2}^{#bar{p}}/#kappa_{1}^{#bar{p}}",
"#kappa_{3}^{#bar{p}}/#kappa_{1}^{#bar{p}}", "(#kappa_{2}/#kappa_{2,Sk} - 1)/#kappa_{2,Sk}"};
const char* obs_ax_2[]{"#kappa_{6}/#kappa_{2}", "#kappa_{4}/#kappa_{2}", "#kappa_{2}/#kappa_{2,Sk}", "#kappa_{2}^{#bar{p}}/#kappa_{1}^{#bar{p}}",
"#kappa_{3}^{#bar{p}}/#kappa_{1}^{#bar{p}}", "(#kappa_{2}/#kappa_{2,Sk} - 1)/#kappa_{2,Sk}"};

const char* obs_sr[]{"k6k2", "k4k2", "k2k2sk", "k2k1", "k3k1", "k2k2sk-1"};
const char* fname_ = "18";
const char* fname_m[]{"18_28_volF_23", "18_28_volF_23", "18_28", "18_28", "18_28", "18_28"};
const char* fname_dat[]{"", "", "", "_singleParticleHighOrder", "_singleParticleHighOrder", ""};
const double yax_lim[][2] = {{0.15, 1.2}, {0.805, 1.05}, {0.89, 1.015}, {0.95, 1.02}, {0.86, 1.05}, {-1.2, 0.5}};
const double yax_lim_2[][2] = {{0.15, 1.2}, {0.805, 1.05}, {0.89, 1.015}, {0.95, 1.02}, {0.86, 1.05}, {-1.4, 0.3}};
const char* volf_str[]{" + vol. fluct.", " + vol. fluct.", "", "", "", ""};

const char* classes[]{"0-10", "10-20", "20-30", "30-40", "40-50", "50-70", "70-100", "0-0.1"};
const double classes_mult[]{18.68, 12.90, 10.03, 7.95, 6.32, 4.49, 2.54, 31.53};
const double classes_mult_err[]{0.16, 0.12, 0.13, 0.11, 0.08, 0.06, 0.04, 0.36};
// const double scale[]{1., 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 1.};

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

void PlotCumulantsVsAcc(const int obs = 5){
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile *f_mb = TFile::Open(Form("../ResultsNetP/final_plots_08_hadPID/out_sys_MB_%s_08_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm = TFile::Open(Form("../ResultsNetP/final_plots_08_hadPID/out_sys_HM_%s_08_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_mb_6 = TFile::Open(Form("../ResultsNetP/final_plots_06_hadPID/out_sys_MB_%s_06_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm_6 = TFile::Open(Form("../ResultsNetP/final_plots_06_hadPID/out_sys_HM_%s_06_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_mb_4 = TFile::Open(Form("../ResultsNetP/final_plots_04_hadPID/out_sys_MB_%s_04_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm_4 = TFile::Open(Form("../ResultsNetP/final_plots_04_hadPID/out_sys_HM_%s_04_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_mb_2 = TFile::Open(Form("../ResultsNetP/final_plots_02_hadPID/out_sys_MB_%s_02_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm_2 = TFile::Open(Form("../ResultsNetP/final_plots_02_hadPID/out_sys_HM_%s_02_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));

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
  TGraphErrors* gKappaVsEta[8];
  TGraphErrors* gKappaVsEtaSys[8];
  for (int i{0}; i < 1; ++i) {
    gKappaVsEta[7] = new TGraphErrors();
    gKappaVsEta[7]->SetName(Form("gKE_%d", i));
    gKappaVsEtaSys[7] = new TGraphErrors();
    gKappaVsEtaSys[7]->SetName(Form("gKE_%d_sys", i));

    auto hSys = (TH1D*)f_hm->Get(Form("hSys_%d", i));
    g_hm_s->AddPoint(g_hm->GetPointX(i), g_hm->GetPointY(i));
    g_hm_s->SetPointError(i, 0.5, hSys->GetStdDev());
    auto hSys_6 = (TH1D*)f_hm_6->Get(Form("hSys_%d", i));
    g_hm_s_6->AddPoint(g_hm_6->GetPointX(i), g_hm_6->GetPointY(i));
    g_hm_s_6->SetPointError(i, 0.5, hSys_6->GetStdDev());
    auto hSys_4 = (TH1D*)f_hm_4->Get(Form("hSys_%d", i));
    g_hm_s_4->AddPoint(g_hm_4->GetPointX(i), g_hm_4->GetPointY(i));
    g_hm_s_4->SetPointError(i, 0.5, hSys_4->GetStdDev());
    auto hSys_2 = (TH1D*)f_hm_2->Get(Form("hSys_%d", i));
    g_hm_s_2->AddPoint(g_hm_2->GetPointX(i), g_hm_2->GetPointY(i));
    g_hm_s_2->SetPointError(i, 0.5, hSys_2->GetStdDev());

    gKappaVsEta[7]->AddPoint(0.4, g_hm_2->GetPointY(i));
    gKappaVsEta[7]->SetPointError(0, 0.0, g_hm_2->GetErrorY(i));
    gKappaVsEta[7]->AddPoint(0.8, g_hm_4->GetPointY(i));
    gKappaVsEta[7]->SetPointError(1, 0.0, g_hm_4->GetErrorY(i));
    gKappaVsEta[7]->AddPoint(1.2, g_hm_6->GetPointY(i));
    gKappaVsEta[7]->SetPointError(2, 0.0, g_hm_6->GetErrorY(i));
    gKappaVsEta[7]->AddPoint(1.6, g_hm->GetPointY(i));
    gKappaVsEta[7]->SetPointError(3, 0.0, g_hm->GetErrorY(i));

    gKappaVsEtaSys[7]->AddPoint(0.4, g_hm_2->GetPointY(i));
    gKappaVsEtaSys[7]->SetPointError(0, 0.02, hSys_2->GetStdDev());
    gKappaVsEtaSys[7]->AddPoint(0.8, g_hm_4->GetPointY(i));
    gKappaVsEtaSys[7]->SetPointError(1, 0.02, hSys_4->GetStdDev());
    gKappaVsEtaSys[7]->AddPoint(1.2, g_hm_6->GetPointY(i));
    gKappaVsEtaSys[7]->SetPointError(2, 0.02, hSys_6->GetStdDev());
    gKappaVsEtaSys[7]->AddPoint(1.6, g_hm->GetPointY(i));
    gKappaVsEtaSys[7]->SetPointError(3, 0.02, hSys->GetStdDev());

    SetGraphStyle(gKappaVsEta[7], colors[0]);
    SetGraphStyle(gKappaVsEtaSys[7], colors[0]);
  }

  for (int i{0}; i < 7; ++i) {
    gKappaVsEta[i] = new TGraphErrors();
    gKappaVsEta[i]->SetName(Form("gKE_%d", i + 1));
    gKappaVsEtaSys[i] = new TGraphErrors();
    gKappaVsEtaSys[i]->SetName(Form("gKE_%d_sys", i + 1));

    auto hSys = (TH1D*)f_mb->Get(Form("hSys_%d", i));
    g_mb_s->AddPoint(g_mb->GetPointX(i), g_mb->GetPointY(i));
    double min = hSys->GetBinCenter(hSys->FindFirstBinAbove(0));
    double max = hSys->GetBinCenter(hSys->FindLastBinAbove(0));
    g_mb_s->SetPointError(i, 0.5, /*0.5 * (max - min)*/ hSys->GetStdDev());
    auto hSys_6 = (TH1D*)f_mb_6->Get(Form("hSys_%d", i));
    g_mb_s_6->AddPoint(g_mb_6->GetPointX(i), g_mb_6->GetPointY(i));
    g_mb_s_6->SetPointError(i, 0.5, hSys_6->GetStdDev());
    auto hSys_4 = (TH1D*)f_mb_4->Get(Form("hSys_%d", i));
    g_mb_s_4->AddPoint(g_mb_4->GetPointX(i), g_mb_4->GetPointY(i));
    g_mb_s_4->SetPointError(i, 0.5, hSys_4->GetStdDev());
    auto hSys_2 = (TH1D*)f_mb_2->Get(Form("hSys_%d", i));
    g_mb_s_2->AddPoint(g_mb_2->GetPointX(i), g_mb_2->GetPointY(i));
    g_mb_s_2->SetPointError(i, 0.5, hSys_2->GetStdDev());

    double scale = 1.; // g_hm_2->GetPointY(0) / g_mb_2->GetPointY(i);

    gKappaVsEta[i]->AddPoint(0.4, scale * g_mb_2->GetPointY(i));
    gKappaVsEta[i]->SetPointError(0, 0.0, scale * g_mb_2->GetErrorY(i));
    gKappaVsEta[i]->AddPoint(0.8, scale * g_mb_4->GetPointY(i));
    gKappaVsEta[i]->SetPointError(1, 0.0, scale * g_mb_4->GetErrorY(i));
    gKappaVsEta[i]->AddPoint(1.2, scale * g_mb_6->GetPointY(i));
    gKappaVsEta[i]->SetPointError(2, 0.0, scale * g_mb_6->GetErrorY(i));
    gKappaVsEta[i]->AddPoint(1.6, scale * g_mb->GetPointY(i));
    gKappaVsEta[i]->SetPointError(3, 0.0, scale * g_mb->GetErrorY(i));

    gKappaVsEtaSys[i]->AddPoint(0.4, scale * g_mb_2->GetPointY(i));
    gKappaVsEtaSys[i]->SetPointError(0, 0.02, scale * hSys_2->GetStdDev());
    gKappaVsEtaSys[i]->AddPoint(0.8, scale * g_mb_4->GetPointY(i));
    gKappaVsEtaSys[i]->SetPointError(1, 0.02, scale * hSys_4->GetStdDev());
    gKappaVsEtaSys[i]->AddPoint(1.2, scale * g_mb_6->GetPointY(i));
    gKappaVsEtaSys[i]->SetPointError(2, 0.02, scale * hSys_6->GetStdDev());
    gKappaVsEtaSys[i]->AddPoint(1.6, scale * g_mb->GetPointY(i));
    gKappaVsEtaSys[i]->SetPointError(3, 0.02, scale * hSys->GetStdDev());

    SetGraphStyle(gKappaVsEta[i], colors[i + 1]);
    SetGraphStyle(gKappaVsEtaSys[i], colors[i + 1]);
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
  txt.DrawLatex(0.18, 0.85, "pp, #sqrt{#it{s}} = 13 TeV INEL > 0");
  txt.DrawLatex(0.18, 0.8, "0.5 < #it{p}_{T} < 1.5 GeV/#it{c}");
  c.Print(Form("c%s_deta.pdf", obs_sr[obs]));

  TFile *fout = TFile::Open("plot_out.root", "recreate");
  fout->cd();
  c.Write();

  // ---------------------------------------------- //
  TCanvas c2("c", "c", 600, 600);
  c2.SetTopMargin(0.03);
  c2.SetRightMargin(0.03);
  c2.SetLeftMargin(0.14);
  c2.SetBottomMargin(0.14);
  TH2D hFrame2("hFrame", Form(";#Delta#eta;%s", obs_ax_2[obs]), 1, 0, 2.6, 100, yax_lim_2[obs][0], yax_lim_2[obs][1]);
  hFrame2.GetXaxis()->SetTitleFont(45);
  hFrame2.GetXaxis()->SetTitleSize(30);
  hFrame2.GetXaxis()->SetTitleOffset(1.1);
  hFrame2.GetYaxis()->SetTitleFont(45);
  hFrame2.GetYaxis()->SetTitleSize(30);
  hFrame2.GetYaxis()->SetTitleOffset(1.3);
  c2.cd();
  hFrame2.Draw();

  TLegend leg_2(0.73913, 0.293913, 0.837793, 0.793043);
  leg_2.SetTextFont(45);
  leg_2.SetTextSize(25);
  leg_2.SetBorderSize(0);
  leg_2.SetHeader("V0M class");

  gKappaVsEta[7]->Write();
  gKappaVsEtaSys[7]->Write();

  gKappaVsEtaSys[7]->Draw("pe5same");
  gKappaVsEta[7]->Draw("pezsame");
  leg_2.AddEntry(gKappaVsEta[7], Form("%s%%", classes[7]), "pe");

  for (int i{0}; i < 7; ++i) {
    gKappaVsEta[i]->Write();
    gKappaVsEtaSys[i]->Write();

    gKappaVsEtaSys[i]->Draw("pe5same");
    gKappaVsEta[i]->Draw("pezsame");
    leg_2.AddEntry(gKappaVsEta[i], Form("%s%%", classes[i]), "pe");
  }

  txt.DrawLatex(0.18, 0.9, "ALICE Preliminary");
  txt.DrawLatex(0.18, 0.85, "pp, #sqrt{#it{s}} = 13 TeV INEL > 0");
  txt.DrawLatex(0.18, 0.19, "0.5 < #it{p}_{T} < 1.5 GeV/#it{c}");

  leg_2.Draw("same");

  c2.Write();
  c2.Print(Form("c%s_deta_vsacc.pdf", obs_sr[obs]));

  // ---------------------------------------------- //
  TCanvas c3("c3", "c3", 800, 600);
  c3.SetTopMargin(0.);
  c3.SetRightMargin(0.);
  c3.SetLeftMargin(0.);
  c3.SetBottomMargin(0.);

  c3.Divide(3, 3);

  c3.cd(1);
  auto pad = c3.cd(1);
  pad->SetTopMargin(0.01);
  pad->SetRightMargin(0.01);
  pad->SetLeftMargin(0.19);
  pad->SetBottomMargin(0.17);

  gKappaVsEta[7]->SetTitle(Form(";#Delta#eta;%s", obs_ax_2[obs]));
  gKappaVsEta[7]->GetXaxis()->SetTitleFont(45);
  gKappaVsEta[7]->GetXaxis()->SetTitleSize(16);
  gKappaVsEta[7]->GetXaxis()->SetLabelFont(45);
  gKappaVsEta[7]->GetXaxis()->SetLabelSize(9);
  gKappaVsEta[7]->GetXaxis()->SetTitleOffset(.8);
  gKappaVsEta[7]->GetYaxis()->SetTitleFont(45);
  gKappaVsEta[7]->GetYaxis()->SetTitleSize(16);
  gKappaVsEta[7]->GetYaxis()->SetLabelFont(45);
  gKappaVsEta[7]->GetYaxis()->SetLabelSize(9);
  gKappaVsEta[7]->GetYaxis()->SetTitleOffset(1.5);
  gKappaVsEta[7]->Draw("apez");
  gKappaVsEtaSys[7]->Draw("pe5same");

  txt.SetTextSize(18);
  txt.SetTextAlign(obs > 4 ? 11 : 31);
  txt.DrawLatex(obs > 4 ? 0.25 : 0.95, 0.86, Form("V0M %s%%", classes[7]));
  txt.SetTextSize(12);
  txt.DrawLatex(obs > 4 ? 0.25 : 0.95, 0.76, Form("#LTd#it{N}_{ch}/d#eta#GT_{|#eta|<0.5} = %.2f #pm %.2f", classes_mult[7], classes_mult_err[7]));

  for (int i{0}; i < 7; ++i) {
    auto pad = c3.cd(i + 2);
    pad->SetTopMargin(0.01);
    pad->SetRightMargin(0.01);
    pad->SetLeftMargin(0.19);
    pad->SetBottomMargin(0.17);

    gKappaVsEta[i]->SetTitle(Form(";#Delta#eta;%s", obs_ax_2[obs]));
    gKappaVsEta[i]->GetXaxis()->SetTitleFont(45);
    gKappaVsEta[i]->GetXaxis()->SetTitleSize(16);
    gKappaVsEta[i]->GetXaxis()->SetTitleOffset(.8);
    gKappaVsEta[i]->GetXaxis()->SetLabelFont(45);
    gKappaVsEta[i]->GetXaxis()->SetLabelSize(9);
    gKappaVsEta[i]->GetYaxis()->SetTitleFont(45);
    gKappaVsEta[i]->GetYaxis()->SetTitleSize(16);
    gKappaVsEta[i]->GetYaxis()->SetLabelFont(45);
    gKappaVsEta[i]->GetYaxis()->SetLabelSize(9);
    gKappaVsEta[i]->GetYaxis()->SetTitleOffset(1.5);
    gKappaVsEta[i]->Draw("apez");
    gKappaVsEtaSys[i]->Draw("pe5same");

    txt.SetTextSize(18);
    txt.DrawLatex(obs > 4 ? 0.25 : 0.95, 0.86, Form("V0M %s%%", classes[i]));
    txt.SetTextSize(12);
    txt.DrawLatex(obs > 4 ? 0.25 : 0.95, 0.76, Form("#LTd#it{N}_{ch}/d#eta#GT_{|#eta|<0.5} = %.2f #pm %.2f", classes_mult[i], classes_mult_err[i]));
  }
  c3.cd(9);
  txt.SetTextAlign(11);
  txt.SetTextSize(18);
  txt.DrawLatex(0.15, 0.8, "ALICE Preliminary");
  txt.DrawLatex(0.15, 0.6, "pp, #sqrt{#it{s}} = 13 TeV INEL > 0");
  txt.DrawLatex(0.15, 0.4, "0.5 < #it{p}_{T} < 1.5 GeV/#it{c}");

  c3.Write();
  c3.Print(Form("c%s_deta_divide.pdf", obs_sr[obs]));

  fout->Close();
}
