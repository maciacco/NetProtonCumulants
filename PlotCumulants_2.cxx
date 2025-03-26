const char* obs_ax[]{"#it{#kappa}_{6}/#it{#kappa}_{2}", "#it{#kappa}_{4}/#it{#kappa}_{2}", "#it{#kappa}_{2}/#it{#kappa}_{2,Sk}", "#it{#kappa}_{2}^{#bar{p}}/#it{#kappa}_{1}^{#bar{p}}", "#it{#kappa}_{3}^{#bar{p}}/#it{#kappa}_{1}^{#bar{p}}"};
const char* obs_sr[]{"k6k2", "k4k2", "k2k2sk", "k2k1", "k3k1"};
const char* obs_sr_m[]{"k6k2", "k4k2", "k2k2sk", "k2k1_n", "k3k1_n"};
const char* fname_ = "18";
const char* fname_m[]{"18_278", "18_278", "18_278", "18_278", "18_278"};
const char* fname_dat[]{"", "", "", "_singleParticleHighOrder", "_singleParticleHighOrder"};
const double yax_lim[][2] = {{0.2, 1.05}, {0.815, 1.01}, {0.7, 1.05}, {0.948, 1.016}, {0.84, 1.05}};
const char* volf_str[]{" + vol. fluct.", " + vol. fluct.", "", "", ""};
//const double xerr[]{0.16, 0.12, 0.13, 0.11, 0.085, 0.064, 0.036};
//const double xerr_hm{0.26};
const double pythia8_monash[][7]{{0.803634, 0.994558, 0.956312, 0.967531, 0.956223, 0.952725, 1.11756}, {0.950518, 0.991494, 0.980016, 0.981019, 0.979718, 0.977441, 1.00332}, {0.833258, 0.802092, 0.785658, 0.771072, 0.756551, 0.742023, 0.728226/*, 0.721047*/}};
const double pythia8_monash_err[][7]{{0.00347729, 0.00917395, 0.0125335, 0.0207804, 0.0169664, 0.0169349, 0.0255032}, {0.000502528, 0.00119946, 0.0016136, 0.00134277, 0.00159809, 0.00111598, 0.0016718}, {0.000128186, 0.000121104, 0.000301248, 0.000184409, 0.000225784, 0.000243093, 0.000224489/*, 0.00212855*/}};
const double pythia8_monash_mult[]{2.70857, 4.74809, 6.93366, 8.87174, 11.4479, 14.8722, 20.4656};
const double pythia8_monash_mult_err[]{0., 0., 0., 0., 0., 0., 0.};

const double mult_5[]{21.18, 16.18, 13.78, 12.01, 10.6742, 9.45727, 8.39124, 7.4559, 6.63397, 5.91062, 5.27308, 4.71037, 4.21301, 3.77281, 3.38268, 3.03647, 2.72883, 2.45514, 2.21133, 1.9939};
const double mult[]{18.68, 12.90, 10.03, 7.95, 6.32, 4.49, 2.54};

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

void PlotCumulants_2(const int obs = 1){
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile *f_mb = TFile::Open(Form("../ResultsNetP/final_plots_08_hadPID/out_sys_MB_%s_08_pct_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_mb_5 = TFile::Open(Form("../ResultsNetP/final_plots_08_hadPID/out_sys_MB_%s_08_5pct_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));

  TGraphErrors gSys_mb;
  TCanvas c("c", "c", 600, 600);
  c.SetTopMargin(0.03);
  c.SetRightMargin(0.03);
  c.SetLeftMargin(0.14);
  c.SetBottomMargin(0.14);
  TH2D hFrame("hFrame", Form(";#LTd#it{N}/d#it{#eta}#GT_{|#it{#eta}|<0.5};%s", obs_ax[obs]), 1, 0, 24., 100, yax_lim[obs][0], yax_lim[obs][1]);
  hFrame.GetXaxis()->SetTitleFont(45);
  hFrame.GetXaxis()->SetTitleSize(30);
  hFrame.GetXaxis()->SetTitleOffset(1.1);
  hFrame.GetYaxis()->SetTitleFont(45);
  hFrame.GetYaxis()->SetTitleSize(30);
  hFrame.GetYaxis()->SetTitleOffset(1.3);
  TGraphErrors* g_mb = (TGraphErrors*)f_mb->Get("g_364");
  TGraphErrors* g_mb_s = new TGraphErrors();
  TGraphErrors* g_mb_5 = (TGraphErrors*)f_mb_5->Get("g_364");
  TGraphErrors* g_mb_s_5 = new TGraphErrors();

  for (int i{0}; i < 7; ++i) {
    auto hSys = (TH1D*)f_mb->Get(Form("hSys_%d", i));
    g_mb->SetPointX(i, mult[i]);
    g_mb_s->AddPoint(g_mb->GetPointX(i), g_mb->GetPointY(i));
    double min = hSys->GetBinCenter(hSys->FindFirstBinAbove(0));
    double max = hSys->GetBinCenter(hSys->FindLastBinAbove(0));
    g_mb_s->SetPointError(i, 0.15, /*xerr[i],*/ /*0.5 * (max - min)*/ hSys->GetStdDev());
  }

  for (int i{0}; i < 20; ++i) {
    auto hSys = (TH1D*)f_mb_5->Get(Form("hSys_%d", i));
    g_mb_5->SetPointX(i, mult_5[i]);
    g_mb_s_5->AddPoint(g_mb_5->GetPointX(i), g_mb_5->GetPointY(i));
    double min = hSys->GetBinCenter(hSys->FindFirstBinAbove(0));
    double max = hSys->GetBinCenter(hSys->FindLastBinAbove(0));
    g_mb_s_5->SetPointError(i, 0.15, /*xerr[i],*/ /*0.5 * (max - min)*/ hSys->GetStdDev());
  }

  SetGraphStyle(g_mb);
  SetGraphStyle(g_mb_s);

  SetGraphStyle(g_mb_5, kBlue);
  SetGraphStyle(g_mb_s_5, kBlue);

  c.cd();
  hFrame.Draw();
//  g_mb_m->Draw("samee3l");
  g_mb_s->Draw("samee5");
  g_mb->Draw("samepez");

  g_mb_s_5->Draw("samee5");
  g_mb_5->Draw("samepez");

  TLegend leg(0.18, 0.72, 0.4, 0.84);
  leg.SetTextFont(45);
  leg.SetTextSize(17);
  leg.SetBorderSize(0);
  leg.AddEntry(g_mb, "Data, [0-10]%, [10-20]%,..., [50-70]%, [70-100]% V0M bins", "pe");
  leg.AddEntry(g_mb_5, "Data, 5% V0M bins", "pe");
  leg.Draw("same");

  TLatex txt;
  txt.SetNDC();
  txt.SetTextFont(45);
  txt.SetTextSize(25);
  txt.DrawLatex(0.18, 0.91, "ALICE Preliminary");
  txt.DrawLatex(0.18, 0.86, "pp, #sqrt{#it{s}} = 13 TeV, INEL > 0");
  txt.DrawLatex(0.18, 0.18, "|#it{#eta}| < 0.8, 0.5 < #it{p}_{T} < 1.5 GeV/#it{c}");
  c.Print(Form("c%s_pct.pdf", obs_sr[obs]));

  TFile *fout = TFile::Open("plot_out.root", "recreate");
  fout->cd();
  c.Write();
  fout->Close();
}
