/*
  several function used for PbPb combined spectra
  Blast Wave is also implemented here
  further documentation will come

  author: Roberto Preghenella
  email : preghenella@bo.infn.it
*/


TH1 *
ReturnExtremeHighHisto(TH1 *hin, TH1 *herr = NULL)
{
  TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehigh", hin->GetName()));
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    if (herr) err = herr->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val + err);
  }
  return hout;
}

TH1 *
ReturnExtremeLowHisto(TH1 *hin, TH1 *herr = NULL)
{
  TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremelow", hin->GetName()));
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    if (herr) err = herr->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val - err);
  }
  return hout;
}

TH1 *
ReturnExtremeHisto(TH1 *hin, TH1 *herr = NULL, Float_t sign = 1.)
{
  Double_t ptlow, pthigh;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    ptlow = hin->GetBinLowEdge(ibin + 1);
    break;
  }
  for (Int_t ibin = hin->GetNbinsX(); ibin >= 0; ibin--) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    pthigh = hin->GetBinLowEdge(ibin + 2);
    break;
  }

  Double_t mean = hin->GetMean();
  Double_t maxdiff = 0.;
  TH1 *hmax = NULL;
  for (Int_t inode = 0; inode < hin->GetNbinsX(); inode++) {

    Double_t ptnode = hin->GetBinCenter(inode + 1);
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehard", hin->GetName()));

    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      Double_t val = hin->GetBinContent(ibin + 1);
      Double_t err = hin->GetBinError(ibin + 1);
      if (herr) err = herr->GetBinError(ibin + 1);
      Double_t cen = hin->GetBinCenter(ibin + 1);
      //    err *= -1. - (cen - ptlow) * (-1. - 1.) / (pthigh - ptlow);
      if (cen < ptnode)
        err *= -1. + (cen - ptlow) / (ptnode - ptlow);
      else
        err *= (cen - ptnode) / (pthigh - ptnode);

      hout->SetBinContent(ibin + 1, val + sign * err);
    }

    Double_t diff = TMath::Abs(mean - hout->GetMean());
    if (diff > maxdiff) {
      //      printf("found max at %f\n", ptnode);
      if (hmax) delete hmax;
      hmax = (TH1 *)hout->Clone("hmax");
      maxdiff = diff;
    }
    delete hout;
  }
  return hmax;
}

TH1 *
ReturnExtremeSoftHisto(TH1 *hin, TH1 *herr = NULL)
{
  return ReturnExtremeHisto(hin, herr, -1.);
}

TH1 *
ReturnExtremeHardHisto(TH1 *hin, TH1 *herr = NULL)
{
  return ReturnExtremeHisto(hin, herr, 1.);
}

TGraphErrors *
ReturnExtremeSoftGraph(TGraphErrors *hin, TGraphErrors *herr = NULL)
{
  TGraphErrors *hout = (TGraphErrors *)hin->Clone(Form("%s_extremesoft", hin->GetName()));
  Double_t ptlow, pthigh;
  ptlow = hin->GetX()[0];
  pthigh = hin->GetX()[hin->GetN() - 1];

  for (Int_t ibin = 0; ibin < hin->GetN(); ibin++) {
    Double_t val = hin->GetY()[ibin];
    Double_t err = hin->GetEY()[ibin];
    if (herr) err = herr->GetEY()[ibin];
    Double_t cen = hin->GetX()[ibin];
    err *= 1. + (cen - ptlow) * (-1. - 1.) / (pthigh - ptlow);
    hout->SetPoint(ibin, cen, val + err);
  }
  return hout;
}


TGraphErrors *
ReturnExtremeHardGraph(TGraphErrors *hin, TGraphErrors *herr = NULL)
{
  TGraphErrors *hout = (TGraphErrors *)hin->Clone(Form("%s_extremehard", hin->GetName()));
  Double_t ptlow, pthigh;
  ptlow = hin->GetX()[0];
  pthigh = hin->GetX()[hin->GetN() - 1];

  for (Int_t ibin = 0; ibin < hin->GetN(); ibin++) {
    Double_t val = hin->GetY()[ibin];
    Double_t err = hin->GetEY()[ibin];
    if (herr) err = herr->GetEY()[ibin];
    Double_t cen = hin->GetX()[ibin];
    err *= -1. - (cen - ptlow) * (-1. - 1.) / (pthigh - ptlow);
    hout->SetPoint(ibin, cen, val + err);
  }
  return hout;
}

/*****************************************************************/
// BOLTZMANN
/*****************************************************************/

Double_t
Boltzmann_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm * mt * TMath::Exp(-mt / T);
}

TF1 *
Boltzmann(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{

  TF1 *fBoltzmann = new TF1(name, Boltzmann_Func, 0., 10., 3);
  fBoltzmann->SetParameters(mass, T, norm);
  fBoltzmann->SetParNames("mass", "T", "norm");
  fBoltzmann->FixParameter(0, mass);
  return fBoltzmann;
}

/*****************************************************************/
/* LEVY-TSALLIS */
/*****************************************************************/

Double_t
LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 *
LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{

  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 0., 1.e5);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}

/*****************************************************************/
/* BOLTZMANN-GIBBS BLAST-WAVE */
/*****************************************************************/

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;
Double_t
BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{

  /*
     x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
  */

  Double_t r = x[0];
  Double_t mt = p[0];
  Double_t pt = p[1];
  Double_t beta_max = p[2];
  Double_t temp_1 = 1. / p[3];
  Double_t n = p[4];

  Double_t beta = beta_max * TMath::Power(r, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;
  Double_t rho = TMath::ATanH(beta);
  Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
  if (argI0 > 700.) argI0 = 700.;
  Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
  //  if (argI0 > 100 || argI0 < -100)
  //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
  return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);

}

Double_t
BGBlastWave_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t norm = p[4];

  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., 1.e-6);
  return norm * pt * integral;
}

Double_t
BGBlastWaveRatio_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max_num = p[1];
  Double_t temp_num = p[2];
  Double_t n_num = p[3];
  Double_t norm_num = p[4];
  Double_t beta_max_den = p[5];
  Double_t temp_den = p[6];
  Double_t n_den = p[7];
  Double_t norm_den = p[8];

  if (!fBGBlastWave_Integrand_num)
    fBGBlastWave_Integrand_num = new TF1("fBGBlastWave_Integrand_num", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_num->SetParameters(mt, pt, beta_max_num, temp_num, n_num);
  Double_t integral_num = fBGBlastWave_Integrand_num->Integral(0., 1.);

  if (!fBGBlastWave_Integrand_den)
    fBGBlastWave_Integrand_den = new TF1("fBGBlastWave_Integrand_den", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_den->SetParameters(mt, pt, beta_max_den, temp_den, n_den);
  Double_t integral_den = fBGBlastWave_Integrand_den->Integral(0., 1.);

  return (norm_num / norm_den) * (integral_num / integral_den);
}

Double_t
BGBlastWaveParticleRatio_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass_num = p[0];
  Double_t mass_den = p[1];
  Double_t mt_num = TMath::Sqrt(pt * pt + mass_num * mass_num);
  Double_t mt_den = TMath::Sqrt(pt * pt + mass_den * mass_den);
  Double_t beta_max = p[2];
  Double_t temp = p[3];
  Double_t n = p[4];
  Double_t norm_num = p[5];
  Double_t norm_den = p[6];

  if (!fBGBlastWave_Integrand_num)
    fBGBlastWave_Integrand_num = new TF1("fBGBlastWave_Integrand_num", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_num->SetParameters(mt_num, pt, beta_max, temp, n);
  Double_t integral_num = fBGBlastWave_Integrand_num->Integral(0., 1.);

  if (!fBGBlastWave_Integrand_den)
    fBGBlastWave_Integrand_den = new TF1("fBGBlastWave_Integrand_den", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_den->SetParameters(mt_den, pt, beta_max, temp, n);
  Double_t integral_den = fBGBlastWave_Integrand_den->Integral(0., 1.);

  return (norm_num / norm_den) * (integral_num / integral_den);
}

Double_t
BGBlastWave_Func_OneOverPt(const Double_t *x, const Double_t *p)
{
  /* 1/pt dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t norm = p[4];

  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., 1.e-3);

  return norm * integral;
}

TF1 *
BGBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{

  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 50.);
  return fBGBlastWave;
}

TF1 *
BGBlastWaveRatio(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{

  TF1 *fBGBlastWave = new TF1(name, BGBlastWaveRatio_Func, 0., 10., 9);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max_num", "T_num", "n_num", "norm_num", "beta_max_den", "T_den", "n_den", "norm_den");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 10.);
  fBGBlastWave->SetParLimits(5, 0.01, 0.99);
  fBGBlastWave->SetParLimits(6, 0.01, 1.);
  fBGBlastWave->SetParLimits(7, 0.01, 10.);
  return fBGBlastWave;
}

TF1 *
BGBlastWaveParticleRatio(const Char_t *name, Double_t mass_num, Double_t mass_den, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm_num = 1.e6, Double_t norm_den = 1.e6)
{

  TF1 *fBGBlastWave = new TF1(name, BGBlastWaveParticleRatio_Func, 0., 10., 7);
  fBGBlastWave->SetParameters(mass_num, mass_den, beta_max, temp, n, norm_num, norm_den);
  fBGBlastWave->SetParNames("mass_num", "mass_den", "beta_max", "T", "n", "norm_num", "norm_den");
  fBGBlastWave->FixParameter(0, mass_num);
  fBGBlastWave->FixParameter(1, mass_den);
  fBGBlastWave->SetParLimits(2, 0.01, 0.99);
  fBGBlastWave->SetParLimits(3, 0.01, 1.);
  fBGBlastWave->SetParLimits(4, 0.01, 10.);
  return fBGBlastWave;
}

TF1 *BGBlastWave_OneOverPT(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{

  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func_OneOverPt, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 50.);
  return fBGBlastWave;
}


TF1 *
BGBlastWave_SingleFit(TH1 *h, Double_t mass, Option_t *opt = "")
{

  TF1 *f = BGBlastWave(Form("fBGBW_%s", h->GetName()), mass);
  h->Fit(f);
  h->Fit(f);
  h->Fit(f, opt);
  return f;

}

Int_t nBW;
TF1 *fBGBW[1000];
TGraphErrors *gBW[1000];

void
BGBlastWave_FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  /* beta -> beta_max */
  Double_t beta = par[nBW+0];
  Double_t T = par[nBW+1];
  Double_t n = par[nBW+2];
  Double_t beta_max = 0.5 * (2. + n) * beta;
#if 1
  /* check beta_max */
  if (beta_max >= 1. || beta_max <= 0.) {
    f = kMaxInt;
    return;
  }
  /* check T */
  if (T <= 0.) {
    f = kMaxInt;
    return;
  }
#endif

  Double_t pt, pte, val, vale, func, pull, chi = 0;
  /* loop over all the data */
  for (Int_t iBW = 0; iBW < nBW; iBW++) {
    /* set BGBW parameters */
    fBGBW[iBW]->SetParameter(4, par[iBW]);
    fBGBW[iBW]->SetParameter(1, beta_max);
    fBGBW[iBW]->SetParameter(2, T);
    fBGBW[iBW]->SetParameter(3, n);
    /* loop over all the points */
    for (Int_t ipt = 0; ipt < gBW[iBW]->GetN(); ipt++) {
      pt = gBW[iBW]->GetX()[ipt];
      // cout << "PT = " << pt << endl;
      pte = gBW[iBW]->GetEX()[ipt];
      // cout << "PTE = " << pte << endl;
      val = gBW[iBW]->GetY()[ipt];
      // cout << "VAL = " << val << endl;
      vale = gBW[iBW]->GetEY()[ipt];
      // cout << "VALE = " << vale << endl;
      func = fBGBW[iBW]->Eval(pt);
      //      func = fBGBW[iBW]->Integral(pt - pte, pt + pte);
      pull = (val - func) / vale;
      chi += pull * pull;
    }
  }

  f = chi;
}

TObjArray *
BGBlastWave_GlobalFit(TObjArray *data, Double_t *mass, Double_t profile = 0.5, Bool_t computeCont = kFALSE, Bool_t fixProfile = kFALSE)
{

  /* get data */
  Int_t ndf = 0;
  nBW = data->GetEntries();
  for (Int_t idata = 0; idata < nBW; idata++) {
    gBW[idata] = (TGraphErrors *)data->At(idata);
    gBW[idata]->SetName(Form("gBW%d", idata));
    ndf += gBW[idata]->GetN();
  }

  /* init BG blast-wave functions */
  for (Int_t idata = 0; idata < nBW; idata++) {
    printf("init BG-BlastWave function #%d: mass = %f\n", idata, mass[idata]);
    fBGBW[idata] = BGBlastWave(Form("fBGBW%d", idata), mass[idata]);
  }

  if (computeCont)
    printf("-> compute contours requested\n");

  /* display data */
  TCanvas *cBW = new TCanvas("cBW");
  cBW->Divide(nBW, 1);
  for (Int_t idata = 0; idata < nBW; idata++) {
    cBW->cd(idata + 1);
    gBW[idata]->Draw("ap*");
  }
  cBW->Update();

  /* init minuit: nBW normalizations + 3 (beta, T, n) BG-BlastWave params */
  const Int_t nbwpars = 3;
  const Int_t nfitpars = nBW + nbwpars;
  TMinuit *minuit = new TMinuit(nfitpars);
  minuit->SetFCN(BGBlastWave_FCN);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg);
  for (Int_t idata = 0; idata < nBW; idata++) {
    minuit->mnparm(idata, Form("norm%d", idata), 1.e6, 1., 0., 0., ierflg);
    ndf--;
  }
  //  minuit->mnparm(nBW + 0, "<beta>", 0.55, 0.01, 0., 1., ierflg);
  //  minuit->mnparm(nBW + 1, "T", 0.14, 0.01, 0., 1., ierflg);
  //  minuit->mnparm(nBW + 2, "n", profile, 0.1, 0., 10., ierflg);

  minuit->mnparm(nBW + 0, "<beta>", 0.7, 0.01, 0.2, 0.7, ierflg);
  minuit->mnparm(nBW + 1, "T", 0.07, 0.001, 0.07, 0.2, ierflg);
  minuit->mnparm(nBW + 2, "n", profile, 0.1, 0.6, 13., ierflg);

  ndf -= 3;
  if (fixProfile) {
    minuit->FixParameter(nBW + 2);
    ndf++;
  }

  /* set strategy */
  arglist[0] = 1;
  minuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);

  /* start MIGRAD minimization */
  arglist[0] = 1;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  /* set strategy */
  arglist[0] = 2;
  minuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);

  /* start MIGRAD minimization */
  arglist[0] = 500000;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  /* start IMPROVE minimization */
  arglist[0] = 500000;
  minuit->mnexcm("IMPROVE", arglist, 1, ierflg);

  /* start MINOS */
  arglist[0] = 500000;
  arglist[1] = nBW + 1;
  arglist[2] = nBW + 2;
  arglist[3] = nBW + 3;
  minuit->mnexcm("MINOS", arglist, 4, ierflg);

  /* print results */
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  minuit->mnprin(4, amin);

  /* get parameters */
  Double_t beta, betae, betaeplus, betaeminus, betagcc, temp, tempe, tempeplus, tempeminus, tempgcc, prof, profe, profeplus, profeminus, profgcc;
  minuit->GetParameter(nBW + 0, beta, betae);
  minuit->mnerrs(nBW + 0, betaeplus, betaeminus, betae, betagcc);
  minuit->GetParameter(nBW + 1, temp, tempe);
  minuit->mnerrs(nBW + 1, tempeplus, tempeminus, tempe, tempgcc);
  minuit->GetParameter(nBW + 2, prof, profe);
  minuit->mnerrs(nBW + 2, profeplus, profeminus, profe, profgcc);
  Double_t beta_max = 0.5 * (2. + prof) * beta;
  Double_t norm[1000], norme[1000];
  for (Int_t idata = 0; idata < nBW; idata++)
    minuit->GetParameter(idata, norm[idata], norme[idata]);

  /* printout */
  printf("[x] *********************************\n");
  printf("[x] beta_max = %f\n", beta_max);
  printf("[x] <beta>   = %f +- %f (e+ = %f, e- = %f)\n", beta, betae, betaeplus, betaeminus);
  printf("[x] T        = %f +- %f (e+ = %f, e- = %f)\n", temp, tempe, tempeplus, tempeminus);
  printf("[x] n        = %f +- %f (e+ = %f, e- = %f)\n", prof, profe, profeplus, profeminus);
  printf("[x] chi2     = %f\n", amin);
  printf("[x] ndf      = %d\n", ndf);

  /* 1-sigma contour */
  minuit->SetErrorDef(1);
  TGraph *gCont1 = NULL;
  if (computeCont) gCont1 = (TGraph *) minuit->Contour(50, nBW + 0, nBW + 1);
  if (gCont1) gCont1->SetName("gCont1");

  /* 2-sigma contour */
  minuit->SetErrorDef(4);
  TGraph *gCont2 = NULL;
  //  if (computeCont) gCont2 = (TGraph *) minuit->Contour(50, nBW + 0, nBW + 1);
  if (gCont2) gCont2->SetName("gCont2");

  /* display fits */
  for (Int_t idata = 0; idata < nBW; idata++) {
    cBW->cd(idata + 1);
    fBGBW[idata]->SetParameter(4, norm[idata]);
    fBGBW[idata]->SetParameter(1, beta_max);
    fBGBW[idata]->SetParameter(2, temp);
    fBGBW[idata]->SetParameter(3, prof);
    fBGBW[idata]->Draw("same");
  }
  cBW->Update();

  /* histo params */
  TH1D *hBW = new TH1D("hBW", "", 4, 0., 4.);
  hBW->SetBinContent(1, beta);
  hBW->SetBinError(1, betae);
  hBW->SetBinContent(2, temp);
  hBW->SetBinError(2, tempe);
  hBW->SetBinContent(3, prof);
  hBW->SetBinError(3, profe);
  hBW->SetBinContent(4, amin/ndf);

  /* BW graph */
  TGraphAsymmErrors *gBetaT = new TGraphAsymmErrors();
  gBetaT->SetName("gBetaT");
  gBetaT->SetPoint(0, beta, temp);
  gBetaT->SetPointEXlow(0, TMath::Abs(betaeminus));
  gBetaT->SetPointEXhigh(0, TMath::Abs(betaeplus));
  gBetaT->SetPointEYlow(0, TMath::Abs(tempeminus));
  gBetaT->SetPointEYhigh(0, TMath::Abs(tempeplus));

  /* prepare output array */
  TObjArray *outoa = new TObjArray();
  for (Int_t idata = 0; idata < nBW; idata++) {
    outoa->Add(gBW[idata]);
    outoa->Add(fBGBW[idata]);
  }
  outoa->Add(cBW);
  outoa->Add(hBW);
  outoa->Add(gBetaT);
  if (gCont1) outoa->Add(gCont1);
  if (gCont2) outoa->Add(gCont2);

  return outoa;

}

/*****************************************************************/

Double_t
y2eta(Double_t pt, Double_t mass, Double_t y){
  Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
  return TMath::ASinH(mt / pt * TMath::SinH(y));
}
Double_t
eta2y(Double_t pt, Double_t mass, Double_t eta){
  Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

/*****************************************************************/
