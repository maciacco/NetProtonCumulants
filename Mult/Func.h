#ifndef __FUNC_H__
#define __FUNC_H__

#include <tuple>
#include <TH1D.h>

const std::string amnames[2] = {"A", "M"};

template<class TH>
std::tuple<TH1D*, TH1D*> projectAM(TH const& h, int const xl, int const xu, int const zl = 0, int const zu = -1){
  TH1D* proj[2];
  for (int iM{0}; iM < 2; ++iM){
    proj[iM] = new TH1D(Form("%s_%s", h->GetName(), amnames[iM].data()), ";#it{p}_{T} (GeV/#it{c});Entries", h->GetNbinsY(), 0, h->GetYaxis()->GetBinUpEdge(h->GetNbinsY()));
  }
  int bzero = h->GetYaxis()->FindBin(0.);
  TH1D* proj_tmp = nullptr;
  if (zu > 0 && zl < zu) {
    proj_tmp = (TH1D*)h->ProjectionY("proj", xl, xu, zl, zu);
  }
  else {
    proj_tmp = (TH1D*)h->ProjectionY("proj", xl, xu, zl, zu);
  }
  for (int ib{bzero}; ib < h->GetNbinsY(); ++ib){
    proj[1]->SetBinContent(ib - bzero + 1, proj_tmp->GetBinContent(ib));
    proj[1]->SetBinError(ib - bzero + 1, proj_tmp->GetBinError(ib));
    proj[0]->SetBinContent(ib - bzero + 1, proj_tmp->GetBinContent(2 * bzero - ib - 1));
    proj[0]->SetBinError(ib - bzero + 1, proj_tmp->GetBinError(2 * bzero - ib - 1));
  }
  return {proj[0], proj[1]};
}

template<class TH>
std::tuple<TH2D*, TH2D*> projectAM2D(TH const& h, int const xl, int const xu){
  TH2D* proj[2];
  for (int iM{0}; iM < 2; ++iM){
    proj[iM] = new TH2D(Form("%s_%s", h->GetName(), amnames[iM].data()), ";#it{p}_{T} (GeV/#it{c});Entries", h->GetNbinsY(), 0, h->GetYaxis()->GetBinUpEdge(h->GetNbinsY()), h->GetNbinsZ(), h->GetZaxis()->GetBinLowEdge(1), h->GetZaxis()->GetBinUpEdge(h->GetNbinsZ()));
  }
  int bzero = h->GetYaxis()->FindBin(0.);
  TH2D* proj_tmp = (TH2D*)h->Project3D("zy");
  for (int iby{bzero}; iby < h->GetNbinsY(); ++iby){
    for (int ibz{1}; ibz < h->GetNbinsZ(); ++ ibz){
      proj[1]->SetBinContent(iby - bzero + 1, ibz, proj_tmp->GetBinContent(iby, ibz));
      proj[1]->SetBinError(iby - bzero + 1, ibz, proj_tmp->GetBinError(iby, ibz));
      proj[0]->SetBinContent(iby - bzero + 1, ibz, proj_tmp->GetBinContent(2 * bzero - iby - 1, ibz));
      proj[0]->SetBinError(iby - bzero + 1, ibz, proj_tmp->GetBinError(2 * bzero - iby - 1, ibz));
    }
  }
  return {proj[0], proj[1]};
}

#endif // __FUNC_H__