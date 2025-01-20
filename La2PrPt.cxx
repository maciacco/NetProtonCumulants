#include <TF1.h>
#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TH1D.h>
#include <TLorentzVector.h> // ugly, but TGenPhaseSpace uses this.
#include <TMath.h>
#include <TRandom.h>

constexpr uint64_t kNtrials{10000000};
constexpr double kLaMass = 1.115683;
constexpr double kSigMass = 1.118937;
constexpr double kPrMass = 0.938;
constexpr double kPiMass = 0.139;
constexpr double kPiZeroMass = 0.140;

double beta_max(const double beta_avg, const double n){
  return beta_avg*0.5*(2+n);
}

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

constexpr int kNMultCharged = 7;
constexpr double kMultCharged[] = {/*31.25, */18.68, 12.90, 10.03, 7.95, 6.32, 4.49, 2.54};
constexpr double PrToPi[] = {0.05695641, 0.05770467, 0.05808717, 0.05805706, 0.0581584, 0.05661002, 0.05106217}; // TODO: interpolation {0.05658572}
constexpr double LaToPi[] = {0.04045519, 0.04027806, 0.03957963, 0.03885696, 0.03753549, 0.03515464, 0.02568511}; // {0.04055572}

constexpr int kNMultChargedHM = 1;
constexpr double kMultChargedHM[] = {31.25};
constexpr double PrToPiHM[] = {0.05658572};
constexpr double LaToPiHM[] = {0.04055572};

void La2PrPt(const bool MB = true) {

  std::unique_ptr<TFile> myFile(TFile::Open(MB ? "data/la2prPt.root" : "data/la2prPtHM.root", "RECREATE"));
  for (int iM{0}; iM < (MB ? kNMultCharged : kNMultChargedHM); ++iM) {
    gRandom->SetSeed(1234);
    double mult_tmp = MB ? kMultCharged[iM] : kMultChargedHM[iM];
    //define exponential function to describe distribution of pT
    double Tkin = /* T_kin[kMultClass]; */ 0.1823 - mult_tmp * 0.0006;
    double betaAvg = /* beta_avg[kMultClass]; */ 0.69 - 7.4 / (mult_tmp + 10.4);
    double nBW = /* n[kMultClass]; */ exp(1.28 - 0.035 * mult_tmp) + exp(3.21 - 0.46 * mult_tmp);


    TF1 blastWave("blastWave", BGBlastWave_Func, 0, 10, 5);
    blastWave.SetParameter(0, kLaMass);
    blastWave.SetParameter(1, beta_max(betaAvg, nBW));
    blastWave.SetParameter(2, Tkin);
    blastWave.SetParameter(3, nBW);
    blastWave.SetParameter(4, 10.);

    TF1 blastWaveSig("blastWaveSig", BGBlastWave_Func, 0, 10, 5);
    blastWaveSig.SetParameter(0, kSigMass);
    blastWaveSig.SetParameter(1, beta_max(betaAvg, nBW));
    blastWaveSig.SetParameter(2, Tkin);
    blastWaveSig.SetParameter(3, nBW);
    blastWaveSig.SetParameter(4, 10.);

    TF1 blastWavePr("blastWavePr", BGBlastWave_Func, 0, 10, 5);
    blastWavePr.SetParameter(0, kPrMass);
    blastWavePr.SetParameter(1, beta_max(betaAvg, nBW));
    blastWavePr.SetParameter(2, Tkin);
    blastWavePr.SetParameter(3, nBW);
    blastWavePr.SetParameter(4, 10.);

    // lorentz vectors to save particles
    TLorentzVector mother, motherSig, pip, pim;
    // GenPhaseSpace to generate decay
    TGenPhaseSpace genPiPr;
    //masses of resulting particles
    const double massesDau[2]{kPrMass, kPiMass};
    const double massesDauSig[2]{kPrMass, kPiZeroMass};

    TH1D hPrPt("hPrPtTmp", ";#it{p}_{T} (GeV/#it{c});Entries", 200, 0.f, 10.f);
    TH1D hPrPtSig("hPrPtSigTmp", ";#it{p}_{T} (GeV/#it{c});Entries", 200, 0.f, 10.f);
    TH1D hPrPtPrim("hPrPtPrimTmp", ";#it{p}_{T} (GeV/#it{c});Entries", 200, 0.f, 10.f);

    // simulate decay
    // loop over number of trials

    for (uint64_t i=0; i < kNtrials; ++i) {
      // initialise variables
      float pT_pr = blastWavePr.GetRandom();
      float pT_original_sig = blastWaveSig.GetRandom();
      float pT_original = blastWave.GetRandom();

      float eta = gRandom->Uniform(-0.8, 0.8);
      float phi = gRandom->Uniform(0, TMath::TwoPi());

      // generate mother and decay
      mother.SetPtEtaPhiM(pT_original, eta, phi, kLaMass);
      motherSig.SetPtEtaPhiM(pT_original_sig, eta, phi, kSigMass);

      genPiPr.SetDecay(mother, 2, massesDau);
      genPiPr.Generate();
      float pT_pos = genPiPr.GetDecay(0)->Pt();

      genPiPr.SetDecay(motherSig, 2, massesDauSig);
      genPiPr.Generate();
      float pT_pos_sig = genPiPr.GetDecay(0)->Pt();

      hPrPt.Fill(pT_pos);
      hPrPtSig.Fill(pT_pos_sig);
      hPrPtPrim.Fill(pT_pr);
    }

    // create directon for this pT value
    myFile->cd();
    std::cout << hPrPt.Integral(1, hPrPt.GetNbinsX()) / hPrPt.GetXaxis()->GetBinWidth(1) << std::endl;
    // hPrPt.Scale(3.754e-2 * hPrPt.GetXaxis()->GetBinWidth(1)/hPrPt.Integral(1, hPrPt.GetNbinsX()));
    // hPrPtPrim.Scale(6.825e-2 * hPrPtPrim.GetXaxis()->GetBinWidth(1)/hPrPtPrim.Integral(1, hPrPt.GetNbinsX()));
    hPrPt.Scale((MB ? LaToPi[iM] : LaToPiHM[iM]) * 0.641 * hPrPt.GetXaxis()->GetBinWidth(1) / hPrPt.Integral(1, hPrPt.GetNbinsX()));
    hPrPtSig.Scale((MB ? LaToPi[iM] * 0.272 : LaToPiHM[iM] * 0.275) * 0.5147 * hPrPtSig.GetXaxis()->GetBinWidth(1) / hPrPtSig.Integral(1, hPrPtSig.GetNbinsX()));
    hPrPtPrim.Scale((MB ? PrToPi[iM] : PrToPiHM[iM]) * hPrPtPrim.GetXaxis()->GetBinWidth(1) / hPrPtPrim.Integral(1, hPrPtPrim.GetNbinsX()));

    // save only pt interval of interest
    const int kNBinsPt = 20;
    float kMinPt = 0.5f;
    float kDeltaPt = 0.05f;
    float ptBins[kNBinsPt + 1];
    for (int iB = 0; iB < kNBinsPt + 1; ++iB){
      ptBins[iB] = kMinPt + kDeltaPt * iB;
    }

    TH1D hPrPt2Save(Form("hPrPt_%d", iM), ";#it{p}_{T} (GeV/#it{c});Entries", kNBinsPt, ptBins);
    TH1D hPrPtSig2Save(Form("hPrPtSig_%d", iM), ";#it{p}_{T} (GeV/#it{c});Entries", kNBinsPt, ptBins);
    TH1D hPrPtPrim2Save(Form("hPrPtPrim_%d", iM), ";#it{p}_{T} (GeV/#it{c});Entries", kNBinsPt, ptBins);
    for (int iB{1}; iB < hPrPt2Save.GetNbinsX() + 1; ++iB) {
      hPrPt2Save.SetBinContent(iB, hPrPt.GetBinContent(hPrPt.FindBin(hPrPt2Save.GetBinCenter(iB))));
      hPrPtSig2Save.SetBinContent(iB, hPrPtSig.GetBinContent(hPrPtSig.FindBin(hPrPt2Save.GetBinCenter(iB))));
      hPrPtPrim2Save.SetBinContent(iB, hPrPtPrim.GetBinContent(hPrPt.FindBin(hPrPt2Save.GetBinCenter(iB))));
      hPrPt2Save.SetBinError(iB, hPrPt.GetBinError(hPrPt.FindBin(hPrPt2Save.GetBinCenter(iB))));
      hPrPtSig2Save.SetBinError(iB, hPrPtSig.GetBinError(hPrPtSig.FindBin(hPrPtSig2Save.GetBinCenter(iB))));
      hPrPtPrim2Save.SetBinError(iB, hPrPtPrim.GetBinError(hPrPt.FindBin(hPrPt2Save.GetBinCenter(iB))));
    }
    blastWavePr.Write();
    hPrPt.Write();
    hPrPtSig.Write();
    hPrPtPrim.Write();
    hPrPt2Save.Write();
    hPrPtSig2Save.Write();
    hPrPtPrim2Save.Write();
  }

  myFile->Close();
}
