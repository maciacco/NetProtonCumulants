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
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., 1.e-4);
  return norm * pt * integral;
}

constexpr double betaAvg[] = {0.431924, 0.363985, 0.3252060, 0.2875364, 0.2524384, 0.2028491, 0.1156451, 0.527852};
constexpr double Tkin[] = {0.177070, 0.184004, 0.1837213, 0.1843517, 0.1838040, 0.1812372, 0.1737870, 0.180095};
constexpr double n[] = {1.8086846, 2.407019, 2.893228, 3.483689, 4.209985, 5.714837, 11.61721, 1.199265};

void fitBW(){
  TFile *fHM = TFile::Open("/home/mciacco/Downloads/HEPData-ins1928822-v1-(Anti)proton_spectrum_in_HM_V0M_multiplicity_class.root");
  TFile *fMB = TFile::Open("/home/mciacco/Downloads/HEPData-ins1784041-v1-Table_5(2).root");
  TFile *fout = TFile::Open("foutBW.root", "recreate");
  TGraphErrors *g[10];
  TF1 *fBW[10];
  for (int iM{0}; iM < 8; ++iM) {
    if (iM == 0) {
      g[iM] = static_cast<TGraphErrors*>(fMB->Get("Table 5/Graph1D_y1"));
      TGraphErrors* g1 = static_cast<TGraphErrors*>(fMB->Get("Table 5/Graph1D_y2"));
      TGraphErrors* g2 = static_cast<TGraphErrors*>(fMB->Get("Table 5/Graph1D_y3"));
      for (int iP{0}; iP < g[iM]->GetN(); ++iP) {
        g[0]->SetPoint(iP, g[iM]->GetPointX(iP), g[iM]->GetPointY(iP) + g1->GetPointY(iP) + g2->GetPointY(iP));
      }
    }
    else if (iM == 1) {
      g[iM] = static_cast<TGraphErrors*>(fMB->Get("Table 5/Graph1D_y2"));
      TGraphErrors* g1 = static_cast<TGraphErrors*>(fMB->Get("Table 5/Graph1D_y4"));
      for (int iP{0}; iP < g[iM]->GetN(); ++iP) {
        g[0]->SetPoint(iP, g[iM]->GetPointX(iP), g[iM]->GetPointY(iP) + g1->GetPointY(iP));
      }
    }

    else if (iM > 1 && iM < 7) g[iM] = static_cast<TGraphErrors*>(fMB->Get(Form("Table 5/Graph1D_y%d", iM + 3)));
    else g[iM] = static_cast<TGraphErrors*>(fHM->Get("(Anti)proton spectrum in HM V0M multiplicity class/Graph1D_y1"));
    g[iM]->SetName(Form("graph_%d", iM));
    fBW[iM] = new TF1(Form("fBW_%d", iM), BGBlastWave_Func, 0., 20., 5);
    fBW[iM]->FixParameter(0, 0.93827208943);
    fBW[iM]->SetParameter(1, 0.5);
    fBW[iM]->SetParLimits(1, 0.2, 1.);
    fBW[iM]->SetParameter(2, 0.18);
    fBW[iM]->SetParLimits(2, 0., .3);
    fBW[iM]->SetParLimits(3, .5, 15.);
    if (iM > -1) {
      fBW[iM]->FixParameter(1, (2 + n[iM])/2. * betaAvg[iM]);
      fBW[iM]->FixParameter(2, Tkin[iM]);
      fBW[iM]->FixParameter(3, n[iM]);
    }
    //fBW[iM]->FixParameter(4, 1.);
    g[iM]->Fit(Form("fBW_%d", iM), "QNRM+", "", 0.7, 3.);
    g[iM]->Fit(Form("fBW_%d", iM), "QNRM+", "", 0.7, 3.);

    fBW[iM]->FixParameter(1, fBW[iM]->GetParameter(1));
    fBW[iM]->FixParameter(2, fBW[iM]->GetParameter(2));
    fBW[iM]->FixParameter(3, fBW[iM]->GetParameter(3));
    g[iM]->Fit(Form("fBW_%d", iM), "QSRM+", "", 0.3, 3.);

    std::cout << "f_pt = " << g[iM]->GetFunction(Form("fBW_%d", iM))->Integral(0.5, 1.5, 1.e-5) / g[iM]->GetFunction(Form("fBW_%d", iM))->Integral(0., 20., 1.e-5)  << std::endl;

    fout->cd();
    g[iM]->Write();
    fBW[iM]->Write();
  }
  fout->Close();
}
