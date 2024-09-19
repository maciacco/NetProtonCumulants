#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <Riostream.h>
#include "utils.h"
#include <TRandom3.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TStopwatch.h>

#define FILL_MC

void ReadTree(const char* fname = "newTree", const char* ofname = "LHC18", const int iVarMin = 0, const int iVarMax = 3)
{
  TStopwatch w;
  w.Start();

  TFile *fEffPr = TFile::Open(Form("%s/%s.root", kDataDir, kEffPrFile));
  TFile f(Form("%s/%s.root", kDataDir, fname));

  const int nF = iVarMax - iVarMin;
  TFile *o[N_SAMPLE][nF];
  for (int iS{0}; iS < N_SAMPLE; ++iS){
    for (int iVar{iVarMin}; iVar < iVarMax; ++iVar){
      o[iS][iVar - iVarMin] = new TFile(Form("%s%d_var_%d.root", ofname, iS, iVar), "recreate");
    }
  }

  // track table (w/ MC info)
  float fPt = 0.f;
  uint8_t fEtaMask = 0u;
  int fSelMask = 0;
  float fOuterPID = 0.f;
  float fGenPt = -999.f;
  uint8_t fGenEtaMask = 100u;
  bool fIsReco = 0;
  int32_t fIndexMiniCollTables = 0;

  // collision table
  int64_t icoll = 0;
  uint8_t fZvtxMask = 0u;
  uint8_t fTriggerMask = 0u;
  uint8_t fNtracklets = 0u;
  uint8_t fV0Multiplicity = 0u;

  TTree *t = (TTree*)f.Get("newtree");

  TClonesArray *tracks = new TClonesArray("miniTrack");
  TClonesArray &tr = *tracks;
  t->SetBranchAddress("fZvtxMask", &fZvtxMask);
  t->SetBranchAddress("fTriggerMask", &fTriggerMask);
  t->SetBranchAddress("fNtracklets", &fNtracklets);
  t->SetBranchAddress("fV0Multiplicity", &fV0Multiplicity);
  t->GetBranch("fTracks")->SetAutoDelete(kFALSE);
  t->SetBranchAddress("fTracks", &tracks);

  // Init histos and tuple
  TH1D *hCent[N_SAMPLE];
  #ifdef FILL_MC
    TH1D *hProtonQ1_Gen[2][N_SAMPLE][nF];
    TH1D *hProtonQ11_Gen[2][N_SAMPLE][nF];
    TH1D *hProtonQ1Sq_Gen[2][N_SAMPLE][nF];
  #endif // FILL_MC

  TH1D *hEffPr[2][kNCentBins][kNEtaBins][N_SAMPLE][nF];
  TNtuple *evtTuple[N_SAMPLE][nF];

  for (int iS = 0; iS < N_SAMPLE; ++iS){
    hCent[iS] = new TH1D(Form("hCent"), ";Centrality (%);Entries", kNCentBinsSmall, kCentBinsSmall);
    for (int iVar{iVarMin}; iVar < iVarMax; ++iVar)
    {
      evtTuple[iS][iVar - iVarMin] = new TNtuple(Form("evtTuple_%d", iVar), Form("evtTuple_%d", iVar), "cent:q1pP:q1pN:q2pP:q2pN");
      evtTuple[iS][iVar - iVarMin]->SetDirectory(o[iS][iVar - iVarMin]);
      for (int iC = 0; iC < 2; ++iC){
      #ifdef FILL_MC
        hProtonQ1_Gen[iC][iS][iVar - iVarMin] = new TH1D(Form("h%sProtonQ1_Gen_%d", kAntiMatterLabel[iC], iVar), ";Centrality (%);q_{1}^{K}", kNCentBinsSmall, kCentBinsSmall);
        hProtonQ11_Gen[iC][iS][iVar - iVarMin] = new TH1D(Form("h%sProtonQ11_Gen_%d", kAntiMatterLabel[iC], iVar), ";Centrality (%);q_{1}^{K}", kNCentBinsSmall, kCentBinsSmall);
        hProtonQ1Sq_Gen[iC][iS][iVar - iVarMin] = new TH1D(Form("h%sProtonQ1Sq_Gen_%d", kAntiMatterLabel[iC], iVar), ";Centrality (%);(q_{1}^{K})^{2}", kNCentBinsSmall, kCentBinsSmall);
      #endif // FILL_MC
        for (int iCent = 0; iCent < kNCentBins; ++iCent){
          for (int iEta = 0; iEta < kNEtaBins; ++iEta){
            if (fEffPr){
              hEffPr[iC][iCent][iEta][iS][iVar - iVarMin] = (TH1D*)fEffPr->Get(Form("subsample_%d_var_%d/h%sEff%s_%d_%d_%d", 1, iVar, kAntiMatterLabel[iC], kPartLabel[0], iCent, 0, iVar));
            }
          }
        }
      }
    }
  }

  #ifdef FILL_MC
    TH3F *hGenRecProton[2][nF];
    TH2D *hGenProton[2][nF];
    TH2D *hRecProton[2][nF];
    double ptBins[kNBinsPt + 1];
    for (int iB = 0; iB < kNBinsPt + 1; ++iB){
      ptBins[iB] = kMinPt + kDeltaPt * iB;
    }

    for (int iV{iVarMin}; iV < iVarMax; ++iV){
      for (int iC = 0; iC < 2; ++iC){
        hGenProton[iC][iV - iVarMin] = new TH2D(Form("h%sGenProton_%d", kAntiMatterLabel[iC], iV), ";Centrality (%);#it{p}_{T} (GeV/#it{c})", kNCentBins, kCentBins, kNBinsPt, ptBins);
        hRecProton[iC][iV - iVarMin] = new TH2D(Form("h%sRecProton_%d", kAntiMatterLabel[iC], iV), ";Centrality (%);#it{p}_{T} (GeV/#it{c})", kNCentBins, kCentBins, kNBinsPt, ptBins);
        hGenRecProton[iC][iV - iVarMin] = new TH3F(Form("h%sGenRecProton_%d", kAntiMatterLabel[iC], iV), ";Centrality (%);#it{N}_{gen};#it{N_rec}", kNCentBins, 0, 100, 200, 0, 200, 200, 0, 200);
      }
    }
  #endif // FILL_MC

  // loop variables
  Long64_t nEntries = kLimitSample ? kLimitedSample : t->GetEntries();
  TH1D hCentTmp("hCentTmp", "hCentTmp", kNCentBins, kCentBins);
  TH1D hCentSmallTmp("hCentSmallTmp", "hCentSmallTmp", kNCentBinsSmall, kCentBinsSmall);
  TH1D hEtaTmp("hEtaTmp", "hEtaTmp", kNBinsPt, kMinEta, kMinEta + kDeltaEta * kNBinsPt);

  #ifdef FILL_MC
    double qPr_1_Gen[kNCentBinsSmall][2][N_SAMPLE][nF];
    double qPr_11_Gen[kNCentBinsSmall][2][N_SAMPLE][nF];
    double qPr_1Sq_Gen[kNCentBinsSmall][2][N_SAMPLE][nF];
    double qPr_2_Gen[kNCentBinsSmall][2][N_SAMPLE][nF];
    for (int iV{iVarMin}; iV < iVarMax; ++iV){
      for (int i = 0; i < 2; ++i){
        for (int j = 0; j < kNCentBinsSmall; ++j){
          for (int is = 0; is < N_SAMPLE; ++is){
              qPr_1_Gen[j][i][is][iV - iVarMin] = 0;
              qPr_11_Gen[j][i][is][iV - iVarMin] = 0;
              qPr_1Sq_Gen[j][i][is][iV - iVarMin] = 0;\
          }
        }
      }
    }
  #endif // FILL_MC

  // Event loop
  gRandom->SetSeed(42);
  for (Long64_t i = 0; i < nEntries; ++i){
    const int iS = (int)(gRandom->Rndm() * N_SAMPLE);

    Long64_t e = i;
    if (!(i%1000000)) std::cout << "n_ev = " << i << std::endl;

    Long64_t tentry = t->LoadTree(e);
    t->GetEntry(tentry);

    float cent = fV0Multiplicity;
    if (cent > kMaxCent) continue;
    int ic = hCentTmp.FindBin(cent);
    int ic_sm = hCentSmallTmp.FindBin(cent);

    hCent[iS]->Fill(cent);

    for (int iVar{iVarMin}; iVar < iVarMax; ++iVar)
    {
      int iTPCcls = iVar % kNTPCcls;
      int iChi2TPC = (iVar / kNTPCcls) % kNChi2TPC;
      int iDCAxy = (iVar / kNTPCcls / kNChi2TPC) % kNDCAxy;
      int iDCAz = (iVar / kNTPCcls / kNChi2TPC / kNDCAxy) % kNDCAz;
      // int iITSPID = (iVar / kNTPCcls / kNChi2TPC / kNDCAxy / kNDCAz) % kNITSPID;
      int iTPCPID = (iVar / kNTPCcls / kNChi2TPC / kNDCAxy / kNDCAz / kNITSPID) % kNTPCPID;
      // std::cout << iTPCcls << "\t" << iChi2TPC << "\t" << iDCAxy << "\t" << iDCAz << "\t" << iTPCPID << std::endl;

      #ifdef FILL_MC
        double qPr_1_gen_tmp[] = {0, 0};
        double qPr_2_gen_tmp[] = {0, 0};
        Long64_t nPr_gen[] = {0, 0};
      #endif // FILL_MC
      double qPr_1_tmp[][2] = {{0, 0}, {0, 0}};
      double qPr_1_tmp_update[][2] = {{0, 0}, {0, 0}};
      double qPr_1_sq_tmp[][2] = {{0, 0}, {0, 0}};
      double qPr_2_tmp[][2] = {{0, 0}, {0, 0}};
      Long64_t nPr[] = {0, 0};

      for (int itrk = 0; itrk < tracks->GetEntries(); ++itrk) {

        miniTrack* trk_tmp = (miniTrack*)tracks->At(itrk);

        #ifdef FILL_MC
          if ( std::abs(trk_tmp->fGenPt) > kTOFptCut || std::abs(trk_tmp->fGenPt) < kPtLowLimitPr || trk_tmp->fGenEtaMask > kEtaCut ) continue;
            int im_MC = trk_tmp->fGenPt > 0 ? 1 : 0;
            // std::cout << im_MC << std::endl;
            qPr_1_gen_tmp[im_MC] += 1.;
            nPr_gen[im_MC] += 1;
            hGenProton[im_MC][iVar - iVarMin]->Fill(cent, std::abs(trk_tmp->fGenPt));
        #endif // FILL_MC
        if (
            ( ( (((trk_tmp->fSelMask & kCutTPCcls[iTPCcls]) == kCutTPCcls[iTPCcls]) || ((trk_tmp->fSelMask & kCutTPCcls2[iTPCcls]) == kCutTPCcls2[iTPCcls]) ) && kRequireTPCclsCut[iTPCcls] ) || !kRequireTPCclsCut[iTPCcls] ) &&
            ( ( (((trk_tmp->fSelMask & kCutChi2TPC[iChi2TPC]) == kCutChi2TPC[iChi2TPC]) || ((trk_tmp->fSelMask & kCutChi2TPC2[iChi2TPC]) == kCutChi2TPC2[iChi2TPC]) ) && kRequireChi2TPCCut[iChi2TPC] ) || !kRequireChi2TPCCut[iChi2TPC] ) &&
            ( ( ((trk_tmp->fSelMask & kCutDCAxy[iDCAxy]) == kCutDCAxy[iDCAxy] || (trk_tmp->fSelMask & kCutDCAxy2[iDCAxy]) == kCutDCAxy2[iDCAxy]) && kRequireDCAxyCut[iDCAxy] ) || !kRequireDCAxyCut[iDCAxy] ) &&
            ( ( ((trk_tmp->fSelMask & kCutDCAz[iDCAz]) == kCutDCAz[iDCAz] || (trk_tmp->fSelMask & kCutDCAz2[iDCAz]) == kCutDCAz2[iDCAz]) && kRequireDCAzCut[iDCAz] ) || !kRequireDCAzCut[iDCAz] ) &&
            std::abs(trk_tmp->fPt) > kPtLowLimitPr && std::abs(trk_tmp->fPt) < kTOFptCut &&
            (trk_tmp->fEtaMask < kEtaCut)
          )
        {
          #ifdef FILL_MC
            int im_ = trk_tmp->fPt > 0 ? 1 : 0;
            int ie_ = hEtaTmp.FindBin(trk_tmp->fGenEtaMask);
            double eff_ = fEffPr ? hEffPr[im_][ic - 1][ie_ - 1][iS][iVar - iVarMin]->GetBinContent(hEffPr[im_][ic - 1][ie_ - 1][iS][iVar - iVarMin]->FindBin(std::abs(trk_tmp->fGenPt))) : kDummyEffPr;
            if (!trk_tmp->fIsReco) continue;
            int im_tmp = trk_tmp->fPt > 0 ? 1 : 0;
            hRecProton[im_tmp][iVar - iVarMin]->Fill(cent, std::abs(trk_tmp->fPt), 1./eff_);
          #endif // FILL_MC
          int im = trk_tmp->fPt > 0 ? 1 : 0;
          int ie = hEtaTmp.FindBin(static_cast<float>(trk_tmp->fEtaMask) / 10.f);
          double eff = fEffPr ? hEffPr[im][ic - 1][ie - 1][iS][0]->GetBinContent(hEffPr[im][ic - 1][ie - 1][iS][iVar - iVarMin]->FindBin(std::abs(trk_tmp->fPt))) : kDummyEffPr;
          qPr_1_tmp[im][0] += 1.;
          qPr_1_tmp_update[im][0] += 1.;
          qPr_2_tmp[im][0] += 1.;
          qPr_1_tmp[im][1] += (1. / eff);
          qPr_1_tmp_update[im][1] += (1. / eff);
          qPr_2_tmp[im][1] += (1. / eff / eff);
          nPr[im] += 1;
        }
      }
      for (int iM = 0; iM < 2; ++iM){
      #ifdef FILL_MC
        qPr_1_Gen[ic_sm-1][iM][iS][iVar - iVarMin] += qPr_1_gen_tmp[iM];
        qPr_11_Gen[ic_sm-1][iM][iS][iVar - iVarMin] += (qPr_1_gen_tmp[iM] * qPr_1_gen_tmp[1 - iM]);
        qPr_1Sq_Gen[ic_sm-1][iM][iS][iVar - iVarMin] += (qPr_1_gen_tmp[iM] * qPr_1_gen_tmp[iM]);
        hGenRecProton[iM][iVar - iVarMin]->Fill(cent, nPr_gen[iM], nPr[iM]);
      #endif // FILL_MC
        for (int iCorr = 0; iCorr < 2; ++iCorr){
          double q1_sq = qPr_1_tmp[iM][iCorr] * qPr_1_tmp[iM][iCorr];
          qPr_1_sq_tmp[iM][iCorr] += q1_sq;
          qPr_1_tmp_update[iM][iCorr] = 0.;
        }
      }

      evtTuple[iS][iVar - iVarMin]->Fill(cent, qPr_1_tmp[1][1], qPr_1_tmp[0][1], qPr_2_tmp[1][1], qPr_2_tmp[0][1]);
    }
  }

  // Process output
  for (int iVar{iVarMin}; iVar < iVarMax; ++iVar)
  {
    for (int iS{0}; iS < N_SAMPLE; ++iS){
    #ifdef FILL_MC
      if (isMC || (!isMC && kUseIndex)){
        o[iS][iVar - iVarMin]->mkdir(Form("subsample_%s%d", kSubsampleFlag, iS + 1));
        o[iS][iVar - iVarMin]->cd(Form("subsample_%s%d", kSubsampleFlag, iS + 1));
      }
      else {
        o[iS][iVar - iVarMin]->mkdir(Form("subsample_%s", ofname));
        o[iS][iVar - iVarMin]->cd(Form("subsample_%s", ofname));
      }
      for (int iM = 0; iM < 2; ++iM){
        for (int iC = 0; iC < kNCentBinsSmall; ++iC){
          double ev = hCent[iS]->GetBinContent(iC+1);
          hProtonQ1_Gen[iM][iS][iVar - iVarMin]->SetBinContent(iC + 1, ev > 0 ? qPr_1_Gen[iC][iM][iS][iVar - iVarMin] / ev : 0);
          hProtonQ11_Gen[iM][iS][iVar - iVarMin]->SetBinContent(iC + 1, ev > 0 ? qPr_11_Gen[iC][iM][iS][iVar - iVarMin] / ev : 0);
          hProtonQ1Sq_Gen[iM][iS][iVar - iVarMin]->SetBinContent(iC + 1, ev > 0 ? qPr_1Sq_Gen[iC][iM][iS][iVar - iVarMin] / ev : 0);
        }
        hProtonQ1_Gen[iM][iS][iVar - iVarMin]->Write();
        hProtonQ11_Gen[iM][iS][iVar - iVarMin]->Write();
        hProtonQ1Sq_Gen[iM][iS][iVar - iVarMin]->Write();
      }
      o[iS][iVar - iVarMin]->cd();

      for (int iM = 0; iM < 2; ++iM){
        hGenProton[iM][iVar - iVarMin]->Write();
        hRecProton[iM][iVar - iVarMin]->Write();
        hGenRecProton[iM][iVar - iVarMin]->Write();
      }
    #endif // FILL_MC
      o[iS][iVar - iVarMin]->cd();
      hCent[iS]->Write();
      evtTuple[iS][iVar - iVarMin]->Write();
    }
  }

  for (int iVar{iVarMin}; iVar < iVarMax; ++iVar){
    for (int iS{0}; iS < N_SAMPLE; ++iS){
      o[iS][iVar - iVarMin]->Close();
    }
  }

  w.Stop();
  w.Print();

  f.Close();
}
