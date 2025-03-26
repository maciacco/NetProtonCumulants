const char* obs_ax[]{"#it{#kappa}_{6}/#it{#kappa}_{2}", "#it{#kappa}_{4}/#it{#kappa}_{2}", "#it{#kappa}_{2}/#it{#kappa}_{2,Sk}", "#it{#kappa}_{2}^{#bar{p}}/#it{#kappa}_{1}^{#bar{p}}", "#it{#kappa}_{3}^{#bar{p}}/#it{#kappa}_{1}^{#bar{p}}"};
const char* obs_sr[]{"k6k2", "k4k2", "k2k2sk", "k2k1", "k3k1"};
const char* fname_ = "18";
const char* fname_m[]{"18_278", "18_278", "18_278", "18_278", "18_278"};
const char* fname_dat[]{"", "", "", "_singleParticleHighOrder", "_singleParticleHighOrder"};
const double yax_lim[][2] = {{/*0.15,*/-0.65, 1.7}, {/*0.805*/0.58, 1.13}, {0.845/*0.805*/, 0.94/*1.06*/}, {0.88, 1.02}, {0.86, 1.05}};
const double yax_lim_diff[][2] = {{-0.17, 0.17}, {-0.014, 0.014}, {-0.007, .007}, {0.95, 1.02}, {0.86, 1.05}};
const char* volf_str[]{" + vol. fluct.", " + vol. fluct.", "", "", ""};
//const double xerr[]{0.16, 0.12, 0.13, 0.11, 0.085, 0.064, 0.036};
//const double xerr_hm{0.26};

//const double pythia8_monash[][7]{{0.803634, 0.994558, 0.956312, 0.967531, 0.956223, 0.952725, 1.11756}, {0.950518, 0.991494, 0.980016, 0.981019, 0.979718, 0.977441, 1.00332}, {0.833258, 0.802092, 0.785658, 0.771072, 0.756551, 0.742023, 0.728226/*, 0.721047*/}};
//const double pythia8_monash_err[][7]{{0.00347729, 0.00917395, 0.0125335, 0.0207804, 0.0169664, 0.0169349, 0.0255032}, {0.000502528, 0.00119946, 0.0016136, 0.00134277, 0.00159809, 0.00111598, 0.0016718}, {0.000128186, 0.000121104, 0.000301248, 0.000184409, 0.000225784, 0.000243093, 0.000224489/*, 0.00212855*/}};

const double pythia8_monash[][7]{{1.12506, 1.16207, 1.19859, 1.15557, 0.947055, 1.01417, 0.987916}, {1.01832, 1.03148, 1.03109, 1.02433, 0.995888, 0.982448, 0.951488}, {0.909469, 0.902476, 0.897163, 0.89186, 0.885977, 0.87954, 0.873931}};
const double pythia8_monash_err[][7]{{0.0208973, 0.0290274, 0.057634, 0.0743559, 0.0854517, 0.0950019, 0.138353}, {0.00267924, 0.002859, 0.00721343, 0.0061776, 0.00708425, 0.00687774, 0.0107421}, {0.000412089, 0.000648351, 0.000700927, 0.000700781, 0.000677897, 0.000665452, 0.00060533}};
const double pythia8_monash_mult[]{3.05523, 4.85737, 6.49253, 8.3376, 10.7163, 13.8762, 18.7087};
const double pythia8_monash_mult_err[]{0., 0., 0., 0., 0., 0., 0.};

const double epos[][7]{{0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0.}, {0.953341, 0.959568, 0.980787, 0.965528, 0.963337, 0.971922, 0.967258}};
const double epos_err[][7]{{0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0.}, {0.00550079, 0.00699218, 0.00512526, 0.0117841, 0.0075855, 0.00811657, 0.00879356}};
const double epos_mult[]{3.02264, 5.01595, 7.30462, 8.38705, 10.3089, 12.6717, 18.2672};

const double lqcd[][2]{{-2.9967e-02, 1.0936e-01}, {8.1287e-02, 1.0936e-01}, {0., 0.}, {0., 0.}, {0., 0.}};
const double lqcd_err[][2]{{1.8776e-02, 2.7321e-04}, {2.8665e-03, 2.7321e-04}, {0., 0.}, {0., 0.}, {0., 0.}};

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

void PlotCumulantsDiff_2(const int obs = 2){
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile *f_mb = TFile::Open(Form("../ResultsNetP/final_plots_08_hadPID/out_sys_MB_%s_08_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_hm = TFile::Open(Form("../ResultsNetP/final_plots_08_hadPID/out_sys_HM_%s_08_finalBinning%s_%s.root", fname_, fname_dat[obs], obs_sr[obs]));
  TFile *f_mb_m_fixv = TFile::Open(Form("out_sys_%s_finalBinning_%s.root", fname_m[obs], obs_sr[obs]));
//  TFile *f_mb_m_fixv = TFile::Open(Form("out_sys_%s_vfix_finalBinning_%s.root", fname_m[obs], obs_sr[obs]));
//  TFile *f_mb_m_volf = TFile::Open(Form("out_sys_test-100000_finalBinning_%s.root", obs_sr[obs]));

//  TFile *f_hm_m = TFile::Open(Form("../test_models_thermalFistDowngrade/out_sys_HM_%s_finalBinning_%s_model.root", fname_m, obs_sr[obs]));
  TGraphErrors gSys_mb, gSys_hm;

/*  double v_lqcd = lqcd[obs][0] / lqcd[obs][1];
  double e_lqcd = std::abs(v_lqcd) * std::hypot(lqcd_err[obs][0] / lqcd[obs][0], lqcd_err[obs][1] / lqcd[obs][1]);
  std::cout << v_lqcd << " +/- " << e_lqcd << std::endl;
  TGraphErrors gLQCD;
  gLQCD.SetName("gLQCD");
  gLQCD.AddPoint(31.2, v_lqcd);
  gLQCD.SetPointError(0, 0.5, e_lqcd);
  gLQCD.SetFillColor(kGray + 2);
  gLQCD.SetLineColor(kGray + 2);
  gLQCD.SetLineWidth(2);
  gLQCD.SetFillStyle(3154);
  gLQCD.SetMarkerStyle(20);
  gLQCD.SetMarkerSize(0);
*/
  TGraphErrors gPythia8Monash(7, pythia8_monash_mult, pythia8_monash[obs], pythia8_monash_mult_err, pythia8_monash_err[obs]);
  gPythia8Monash.SetFillColor(kOrange + 1);
  gPythia8Monash.SetLineColor(kOrange + 1);
  gPythia8Monash.SetLineWidth(2);
  gPythia8Monash.SetFillStyle(3154);
  gPythia8Monash.SetMarkerStyle(20);
  gPythia8Monash.SetMarkerSize(0);

  TGraphErrors gEPOS4(7, epos_mult, epos[obs], pythia8_monash_mult_err, epos_err[obs]);
  gEPOS4.SetFillColor(kGreen + 2);
  gEPOS4.SetLineColor(kGreen + 2);
  gEPOS4.SetLineWidth(2);
  gEPOS4.SetFillStyle(3154);
  gEPOS4.SetMarkerStyle(20);
  gEPOS4.SetMarkerSize(0);


/*  TGraphErrors* g_mb_m_volf = (TGraphErrors*)f_mb_m_volf->Get("g_364");
  g_mb_m_volf->SetLineColor(kBlack);
  g_mb_m_volf->SetLineWidth(2);
  g_mb_m_volf->SetLineStyle(kDashed);
  g_mb_m_volf->SetMarkerStyle(20);
  g_mb_m_volf->SetMarkerSize(0);
*/
  TCanvas c("c", "c", 600, 700);

  TPad p1("p1", "p1", 0., .25, 1., 1., 0);
  p1.SetTopMargin(0.03);
  p1.SetRightMargin(0.03);
  p1.SetLeftMargin(0.14);
  p1.SetBottomMargin(0.);
  c.cd();
  p1.Draw();

  TPad p2("p2", "p2", 0., 0., 1., .25, 0);
  p2.SetTopMargin(0.0);
  p2.SetRightMargin(0.03);
  p2.SetLeftMargin(0.14);
  p2.SetBottomMargin(0.4);
  c.cd();
  p2.Draw();

  TH2D hFrame("hFrame", Form(";#LTd#it{N}/d#it{#eta}#GT_{|#it{#eta}|<0.5};%s", obs_ax[obs]), 1, 0, 34., 100, yax_lim[obs][0], yax_lim[obs][1]);
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

  TGraphErrors* g_mb_s = new TGraphErrors();
  TGraphErrors* g_hm_s = new TGraphErrors();

/*  TGraphErrors* g_mb_m = (TGraphErrors*)f_mb_m->Get("g_364");
  g_mb_m->SetFillColor(kAzure + 1);
  g_mb_m->SetLineColor(kAzure + 1);
  g_mb_m->SetLineWidth(2);
  g_mb_m->SetFillStyle(3154);
  g_mb_m->SetMarkerStyle(20);
  g_mb_m->SetMarkerSize(0);
*/
  TGraphErrors* g_mb_m_fixv = (TGraphErrors*)f_mb_m_fixv->Get("g_364");
  g_mb_m_fixv->SetFillColor(kAzure - 1);
  g_mb_m_fixv->SetLineColor(kAzure - 1);
  g_mb_m_fixv->SetLineWidth(2);
  g_mb_m_fixv->SetFillStyle(3154);
  g_mb_m_fixv->SetMarkerStyle(20);
  g_mb_m_fixv->SetMarkerSize(0);

  g_mb_m_fixv->SetPointX(0, 31.53);

//  TGraphErrors* g_hm_m = (TGraphErrors*)f_hm_m->Get("g_gen_0");
  for (int i{0}; i < 7; ++i) {
    auto hSys = (TH1D*)f_mb->Get(Form("hSys_%d", i));
    g_mb_s->AddPoint(g_mb->GetPointX(i), g_mb->GetPointY(i));
    double min = hSys->GetBinCenter(hSys->FindFirstBinAbove(0));
    double max = hSys->GetBinCenter(hSys->FindLastBinAbove(0));
    g_mb_s->SetPointError(i, 0.3, /*xerr[i],*/ /*0.5 * (max - min)*/ hSys->GetStdDev());
  }
  for (int i{0}; i < 1; ++i) {
    auto hSys = (TH1D*)f_hm->Get(Form("hSys_%d", i));
    g_hm_s->AddPoint(g_hm->GetPointX(i), g_hm->GetPointY(i));
    g_hm_s->SetPointError(i, 0.3, /*xerr_hm,*/ hSys->GetStdDev());
  }

  g_hm->SetPointX(0, 31.53);
  g_hm_s->SetPointX(0, 31.53);

  SetGraphStyle(g_mb);
  SetGraphStyle(g_hm);
  SetGraphStyle(g_mb_s);
  SetGraphStyle(g_hm_s);

  TLine lGC(0., 1., 34., 1.);
  lGC.SetLineStyle(kDashed);

//  SetGraphStyleModel(g_mb_m);
//  SetGraphStyleModel(g_hm_m);
  p1.cd();
  hFrame.Draw();
//  g_mb_m->Draw("samee3l");
  gPythia8Monash.Draw("samee3l");
//  gEPOS4.Draw("samee3l");
  g_mb_m_fixv->Draw("samee3l");
//  gLQCD.Draw("samee5l");
//  g_mb_m_volf->Draw("samel");
  g_mb_s->Draw("samee5");
  g_hm_s->Draw("samee5");
  g_mb->Draw("samepez");
  g_hm->Draw("samepez");
  lGC.Draw("same");

//  g_hm_m->Draw("samele");
  TGraph g_ph;
  g_ph.SetLineWidth(0);
  g_ph.SetMarkerStyle(20);
  g_ph.SetMarkerSize(0);
  TLegend leg(0.18, 0.03, 0.35, 0.2983);
  leg.SetTextFont(45);
  leg.SetTextSize(17);
  leg.SetBorderSize(0);

  leg.AddEntry(&g_ph, "");
  leg.AddEntry(g_mb, "Data", "pe");

//  leg.AddEntry(&gLQCD, "Lattice QCD, #it{T} = 155 MeV, PRD 101 (2020) 074502", "f");

//  leg.AddEntry(&g_ph, "");
  leg.AddEntry(&gPythia8Monash, "Pythia 8.313, rope hadronization + QCD CR", "f");
//  leg.AddEntry(&gEPOS4, "EPOS 4.0.0", "f");
  leg.AddEntry(&g_ph);

  leg.AddEntry(&g_ph, ""); // "Thermal-FIST+ Blast Wave, #it{V}_{c} = 2.78 d#it{V}/d#it{y}", "f");
  leg.AddEntry(&g_ph, ""); //"#it{T}_{chem}, d#it{V}/d#it{y}, and #it{#gamma}_{s} from PRC 100 (2019) 054906");
  leg.AddEntry(g_mb_m_fixv, "CE SHM ev. gen.", "f");
//  leg.AddEntry(g_mb_m, "Thermal-FIST + vol. fluct. (PRC 88 (2013) 034911)", "f");
//  leg.AddEntry(g_mb_m_volf, "Thermal-FIST + vol. fluct. (Gaussian sampling)");
  leg.Draw("same");

  TLatex txt;
  txt.SetNDC();
  txt.SetTextFont(45);
  txt.SetTextSize(17);
  txt.DrawLatex(0.18, 0.12, "Thermal-FIST+ Blast Wave, #it{V}_{c} = 2.78 d#it{V}/d#it{y}");
  txt.DrawLatex(0.18, 0.085, "#it{T}_{chem}, d#it{V}/d#it{y}, and #it{#gamma}_{s} from PRC 100 (2019) 054906");

//  TLatex txt;
//  txt.SetNDC();
  txt.SetTextFont(45);
  txt.SetTextSize(25);
  txt.DrawLatex(0.18, 0.9, "ALICE Preliminary");
  txt.DrawLatex(0.18, 0.84, "pp, #sqrt{#it{s}} = 13 TeV, INEL > 0");
  txt.DrawLatex(0.18, 0.78, "|#it{#eta}| < 0.8, 0.5 < #it{p}_{T} < 1.5 GeV/#it{c}");

  p2.cd();
  TH2D hFrame2("hFrame2", ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};Data - FIST", 1, 0, 34., 100, yax_lim_diff[obs][0], yax_lim_diff[obs][1]);
  hFrame2.GetXaxis()->SetTitleFont(45);
  hFrame2.GetXaxis()->SetTitleSize(30);
  hFrame2.GetXaxis()->SetTitleOffset(.85);
  hFrame2.GetYaxis()->SetTitleFont(45);
  hFrame2.GetYaxis()->SetTitleSize(20);
  hFrame2.GetYaxis()->SetTitleOffset(1.8);
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
  TGraphErrors* diff_mb_m = new TGraphErrors(*g_mb_m_fixv);
  for (int i{0}; i < 7; ++i) {
    diff_mb->SetPointY(i, diff_mb->GetPointY(i) - g_mb_m_fixv->GetPointY(i + 1));
    diff_mb_s->SetPointY(i, diff_mb_s->GetPointY(i) - g_mb_m_fixv->GetPointY(i + 1));
    diff_mb_m->SetPointY(i + 1, 0.);
  }
  diff_hm->SetPointY(0, diff_hm->GetPointY(0) - g_mb_m_fixv->GetPointY(0));
  diff_hm_s->SetPointY(0, diff_hm_s->GetPointY(0) - g_mb_m_fixv->GetPointY(0));
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
