#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <limits.h>
#include <memory>
#include <vector>
#include <string>

#include <TColor.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TList.h>
#include <TObject.h>
#include <TString.h>

#include <ROOT/RDataFrame.hxx>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooArgList.h>
#include <RooArgSet.h>

using namespace std;
using ROOT::RDataFrame;

constexpr double fac = 0.;

namespace utils
{

  class TH2TF { // https://root-forum.cern.ch/t/create-a-tf1-from-a-th1f/41429/7
    public:
      void SetInputHist (TH1D* h_in) {h_input = (TH1D*)h_in;}
      double Eval (double *x, double *p) {
        double xx = x[0];
        if (xx < h_input->GetBinLowEdge(1) || xx > h_input->GetBinLowEdge(h_input->GetNbinsX())+ h_input->GetBinWidth(1))
          return 0.;
        return p[0]*h_input->GetBinContent(h_input->FindBin(xx)); // add scaling parameter
      }
    private:
      TH1D* h_input;
  };

  bool replace(std::string &str, const std::string &from, const std::string &to)
  {
    size_t start_pos = str.find(from);
    if (start_pos == std::string::npos)
      return false;
    str.replace(start_pos, from.length(), to);
    return true;
  }

  template <class T>
  void Requires(T *obj, string msg = "")
  {
    if (!obj)
    {
      std::cout << "Missing object: " << msg.data() << "." << std::endl;
      abort();
    }
  }

  class TTList : public TList
  {
  public:
    TObject *Get(std::string name)
    {
      TObject *obj = this->FindObject(name.data());
      Requires(obj, name);
      return obj;
    }
  };

  /// Older version of ROOT 6 do not provide TMath::Sq
  template <typename T>
  constexpr T Sq(T x)
  {
    return ((x) * (x));
  }

  constexpr bool ValueInInterval(const double val, const double min, const double max)
  {
    return (val >= min && val < max);
  }

  double QuadratureSum(std::vector<double> &args)
  {
    double sum = 0.;
    for (double &val : args)
      sum += Sq(val);
    return std::sqrt(sum);
  }

  /// Workaround if C++14 is not available
  template <typename T, typename... Args>
  std::unique_ptr<T> make_unique(Args &&...args)
  {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  double Ztest(const double mu0, const double sig0, const double mu1, const double sig1, const double corr = 0.)
  {
    const double sigma = std::sqrt(sig0 * sig0 + sig1 * sig1 + 2 * corr);
    if (sigma < FLT_MIN * 10.f)
      return FLT_MAX;
    else
      return (mu0 - mu1) / sigma;
  }

  void Divide(TH1 *h, TGraph *gr)
  {
    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
      if (h->GetBinContent(i) == 0 || std::abs(gr->Eval(h->GetBinCenter(i))) < FLT_MIN * 10.f)
      {
        continue;
      }
      h->SetBinContent(i, h->GetBinContent(i) / gr->Eval(h->GetBinCenter(i)));
      h->SetBinError(i, h->GetBinError(i) / gr->Eval(h->GetBinCenter(i)));
    }
  }

  TH1 *ComputeEfficiency(TH1 *tof, TH1 *tot)
  {
    TH1 *efftof = (TH1 *)tof->Clone();
    for (int iBin = 1; iBin <= efftof->GetNbinsX(); ++iBin)
    {
      const double den = tot->GetBinContent(iBin);
      if (std::abs(den) < FLT_MIN * 10.f)
        continue;
      double eff = tof->GetBinContent(iBin) / den;
      efftof->SetBinContent(iBin, eff);
      if (eff>1.e-9)
        efftof->SetBinError(iBin, std::sqrt(eff * (1. - eff) / den));
      else
        efftof->SetBinError(iBin, 0.);
    }
    return efftof;
  }

  void MeanRebin(TH1F *&hOriginal, int n_bins, const double *bin_limits, float Max_val = 100)
  {
    int iBinCounter = 1;
    std::vector<float> content_tmp(n_bins, 0.);
    std::vector<float> err2_tmp(n_bins, 0.);
    std::vector<int> count_tmp(n_bins, 0);
    TH1F *hTmp = (TH1F *)hOriginal->Clone("hTmp");
    const char *original_name = hOriginal->GetName();
    printf("original name: %s\n", original_name);
    const char *original_title = hOriginal->GetTitle();
    printf("original title: %s\n", original_title);
    const char *original_X = hOriginal->GetXaxis()->GetTitle();
    printf("original X: %s\n", original_X);
    const char *original_Y = hOriginal->GetYaxis()->GetTitle();
    printf("original Y: %s\n", original_Y);
    delete hOriginal;
    for (int iB = 1; iB <= hTmp->GetNbinsX(); iB++)
    {
      if (hTmp->GetBinCenter(iB) < bin_limits[0])
        continue;
      if (hTmp->GetBinCenter(iB) > Max_val)
        break;
      if (hTmp->GetBinCenter(iB) < bin_limits[iBinCounter])
      {
        content_tmp[iBinCounter - 1] += hTmp->GetBinContent(iB);
        err2_tmp[iBinCounter - 1] += Sq(hTmp->GetBinError(iB));
        count_tmp[iBinCounter - 1]++;
      }
      else
      {
        iBinCounter++;
        content_tmp[iBinCounter - 1] += hTmp->GetBinContent(iB);
        err2_tmp[iBinCounter - 1] += Sq(hTmp->GetBinError(iB));
        count_tmp[iBinCounter - 1]++;
      }
    }
    hOriginal = new TH1F(original_name, Form("%s;%s;%s", original_title, original_X, original_Y), n_bins, bin_limits);
    hOriginal->SetDirectory(0);
    for (int i = 0; i < n_bins; i++)
    {
      if (hOriginal->GetBinCenter(i + 1) > Max_val)
        continue;
      hOriginal->SetBinContent(i + 1, content_tmp[i] / count_tmp[i]);
      hOriginal->SetBinError(i + 1, TMath::Sqrt(err2_tmp[i]) / count_tmp[i]);
    }
  }

  float zTest(const float mu0, const float sig0, const float mu1, const float sig1)
  {
    const float sigma = sqrt(sig0 * sig0 + sig1 * sig1);
    if (sigma < FLT_MIN * 10.f)
      return FLT_MAX;
    else
      return (mu0 - mu1) / std::abs(sig1 - sig0);
  }

  void SmoothInRange(TH1F *histo, float low_limit, float up_limit, int n_times = 1)
  {
    histo->GetXaxis()->SetRange(histo->FindBin(low_limit), histo->FindBin(up_limit));
    histo->Smooth(n_times, "R");
  }

  Double_t Eff(TH1D *hist_num, TH1D *hist_den, Double_t bin_cent)
  {
    Double_t num = hist_num->GetBinContent(hist_num->GetXaxis()->FindBin(bin_cent));
    Double_t den = hist_den->GetBinContent(hist_den->GetXaxis()->FindBin(bin_cent));
    return num / den;
  }

  Double_t EffErr(TH1D *hist_eff, TH1D *hist_gen, Double_t bin_cent)
  {
    Double_t gen = hist_gen->GetBinContent(hist_gen->GetXaxis()->FindBin(bin_cent));
    Double_t eff = hist_eff->GetBinContent(hist_eff->GetXaxis()->FindBin(bin_cent));
    Double_t q = 1.00000000000 - eff;
    return TMath::Sqrt(gen * eff * q) / gen;
  }

  void MergeGraphPoints(TGraphErrors *gmodel, int bMin, int bMax){
    double x1 = gmodel->GetPointX(bMin);
    double x2 = gmodel->GetPointX(bMax);
    double y1 = gmodel->GetPointY(bMin);
    double y2 = gmodel->GetPointY(bMax);
    double erry1 = gmodel->GetErrorY(bMin);
    double erry2 = gmodel->GetErrorY(bMax);

    double xavg = 0.5 * (x1 + x2);
    double yavg = (y1 / erry1 / erry1 + y2 / erry2 / erry2) / (1. / erry1 / erry1 + 1. / erry2 / erry2);
    double eavg = 1. / sqrt(1. / erry1 / erry1 + 1. / erry2 / erry2);

    gmodel->RemovePoint(bMax);
    gmodel->SetPointX(bMin, xavg);
    gmodel->SetPointY(bMax, yavg);
    gmodel->SetPointError(bMin, 0., eavg);
  }

  Double_t chi2(TGraphErrors *gstat, TGraphErrors *gsyst, TGraphErrors *gmodel){
    double res = 0.;
    for (int i{0}; i < gstat->GetN(); ++i){
      double o = gstat->GetPointY(i);
      double ostat = gstat->GetErrorY(i);
      double osyst = gsyst->GetErrorY(i);
      double model = gmodel->GetPointY(i);
      double modelstat = gmodel->GetErrorY(i);

      double tmp = 0.;
      tmp = ( o + fac * osyst - model ) * ( o + fac * osyst - model );
      tmp = tmp / ( std::pow(ostat / o * (o + fac * osyst), 2) /* + osyst * osyst *//*  + modelstat * modelstat */);
      res += tmp;
    }
    return res;
  }

  Double_t chi2interp(TGraphErrors *gstat, TGraphErrors *gsyst, TGraphErrors *gmodel){
    double res = 0.;
    for (int i{0}; i < gstat->GetN(); ++i){
      double o = gstat->GetPointY(i);
      double ostat = gstat->GetErrorY(i);
      double osyst = gsyst->GetErrorY(i);
      int pL = -999, pR = -999;
      double x = gstat->GetPointX(i);
      for (int iP{0}; iP < gmodel->GetN(); ++iP){
        // std::cout << x << "\t" << gmodel->GetPointX(iP) << "\n";
        if (x > gmodel->GetPointX(iP)) continue;
        pR = iP;
        pL = iP - 1;
        break;
      }
      if (pR == pL) continue;
      TF1 f("fit", "pol1", gmodel->GetPointX(pL) - 0.5, gmodel->GetPointX(pR) + 0.5);
      auto r = gmodel->Fit("fit", "NRSQ", "", gmodel->GetPointX(pL) -0.5, gmodel->GetPointX(pR) + 0.5);
      double xx[1] = { x };
      double err[1];  // error on the function at point x0
      r->GetConfidenceIntervals(1, 1, 1, xx, err, 0.683, false);
      double model = f.Eval(x);
      double model_err = err[0];
      double tmp = 0.;
      std::cout << o << "\t" << model << "\n";
      tmp = ( o - model ) * ( o - model );
      tmp = tmp / ( ostat * ostat + osyst * osyst + model_err * model_err);
      res += tmp;
    }
    return res;
  }

  Double_t chi2err(TGraphErrors *gstat, TGraphErrors *gsyst, TGraphErrors *gmodel){
    double res = 0.;
    for (int i{0}; i < gstat->GetN(); ++i){
      double o = gstat->GetPointY(i);
      double ostat = gstat->GetErrorY(i);
      double osyst = gsyst->GetErrorY(i);
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

}

#endif
