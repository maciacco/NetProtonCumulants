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

void ReadTree(const char* fname = "newTree_mc", const char* ofname = "LHC21d3", const int iVarMin = 364, const int iVarMax = 365)
{
  TStopwatch w;
  w.Start();

  TFile *fEffPr = TFile::Open(Form("%s/%s.root", kResDir, kEffPrFile));
  TFile f(Form("%s/%s.root", kResDir, fname));

  const int nF = iVarMax - iVarMin;
  TFile *o[N_SAMPLE][nF];
  for (int iS{0}; iS < N_SAMPLE; ++iS){
    for (int iVar{iVarMin}; iVar < iVarMax; ++iVar){
      o[iS][iVar - iVarMin] = new TFile(Form("%s/%s%d_var_%d.root", kResDir, ofname, iS, iVar), "recreate");
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
  TH1D *hCent[N_SAMPLE]{nullptr};
  TH1D* hNtrkl[N_SAMPLE]{nullptr};
  TH3D *outerPID[N_SAMPLE][2];
  #ifdef FILL_MC
    TNtupleD *evtTupleGen[N_SAMPLE][nF];
  #endif // FILL_MC

  TH1D *hEffPr[2][kNCentBins][kNEtaBins][N_SAMPLE][nF];
  TNtupleD *evtTuple[N_SAMPLE][nF];

  for (int iS = 0; iS < N_SAMPLE; ++iS){
    hCent[iS] = new TH1D(Form("hCent_%d", iS), ";Centrality (%);Entries", kNCentBinsSmall, kCentBinsSmall);
    hNtrkl[iS] = new TH1D(Form("hNtrkl_%d", iS), ";#it{N}_{tracklets}^{0.7 < |#eta| < 1.2};Entries", 100, 0, 100);

    for (int iVar{iVarMin}; iVar < iVarMax; ++iVar)
    {
      evtTuple[iS][iVar - iVarMin] = new TNtupleD(Form("evtTuple_%d", iVar), Form("evtTuple_%d", iVar), "cent:q1pP:q1pN:q2pP:q2pN:q3pP:q3pN:q4pP:q4pN:q5pP:q5pN:q6pP:q6pN");
      evtTuple[iS][iVar - iVarMin]->SetDirectory(o[iS][iVar - iVarMin]);
      for (int iC = 0; iC < 2; ++iC){
      #ifdef FILL_MC
        evtTupleGen[iS][iVar - iVarMin] = new TNtupleD(Form("evtTupleGen_%d", iVar), Form("evtTupleGen_%d", iVar), "cent:q1pP:q1pN");
        evtTupleGen[iS][iVar - iVarMin]->SetDirectory(o[iS][iVar - iVarMin]);
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

  // histos
  float ptBins[kNBinsPt + 1];
  for (int iB = 0; iB < kNBinsPt + 1; ++iB){
    ptBins[iB] = kMinPt + kDeltaPt * iB;
  }
  float pidBins[kNBinsPID + 1];
  for (int iB = 0; iB < kNBinsPID + 1; ++iB){
    pidBins[iB] = kMinPID + kDeltaPID * iB;
  }
  for (int iS = 0; iS < N_SAMPLE; ++iS){
    for (int iC = 0; iC < 2; ++iC){
      outerPID[iS][iC] = new TH3D(Form("h%sOuterPID_%d", kAntiMatterLabel[iC], iS), ";Centrality (%);#it{p}_{T} (GeV/#it{c});n#sigma (a.u.)", kNCentBins, kCentBins, kNBinsPt, ptBins, kNBinsPID, pidBins);
    }
  }

  #ifdef FILL_MC
    TH3F *hGenRecProton[2][nF];
    TH2D *hGenProton[2][nF];
    TH2D *hRecProton[2][nF];

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
    hNtrkl[iS]->Fill(fNtracklets);

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
        Long64_t nPr_gen[] = {0, 0};
      #endif // FILL_MC
      double qPr_1_tmp[2] = {0, 0};
      double qPr_2_tmp[2] = {0, 0};
      double qPr_3_tmp[2] = {0, 0};
      double qPr_4_tmp[2] = {0, 0};
      double qPr_5_tmp[2] = {0, 0};
      double qPr_6_tmp[2] = {0, 0};
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
            (trk_tmp->fEtaMask < kEtaCut) &&
            (std::abs(trk_tmp->fOuterPID) < 3.f)
          )
        {
          #ifdef FILL_MC
            int im_ = trk_tmp->fPt > 0 ? 1 : 0;
            int ie_ = hEtaTmp.FindBin(static_cast<float>(trk_tmp->fEtaMask) / 10.f);
            double eff_ = fEffPr ? hEffPr[im_][ic - 1][ie_ - 1][iS][iVar - iVarMin]->GetBinContent(hEffPr[im_][ic - 1][ie_ - 1][iS][iVar - iVarMin]->FindBin(std::abs(trk_tmp->fGenPt))) : kDummyEffPr;
            if (!trk_tmp->fIsReco || trk_tmp->fGenPt < -998.f) continue;
            int im_tmp = trk_tmp->fPt > 0 ? 1 : 0;
            hRecProton[im_tmp][iVar - iVarMin]->Fill(cent, std::abs(trk_tmp->fPt), 1./eff_);
          #endif // FILL_MC
          int im = trk_tmp->fPt > 0 ? 1 : 0;
          int ie = hEtaTmp.FindBin(static_cast<float>(trk_tmp->fEtaMask) / 10.f);
          double eff = fEffPr ? hEffPr[im][ic - 1][ie - 1][iS][0]->GetBinContent(hEffPr[im][ic - 1][ie - 1][iS][iVar - iVarMin]->FindBin(std::abs(trk_tmp->fPt))) : kDummyEffPr;
          qPr_1_tmp[im] += (1. / eff);
          qPr_2_tmp[im] += (1. / powI(eff, 2));
          qPr_3_tmp[im] += (1. / powI(eff, 3));
          qPr_4_tmp[im] += (1. / powI(eff, 4));
          qPr_5_tmp[im] += (1. / powI(eff, 5));
          qPr_6_tmp[im] += (1. / powI(eff, 6));
          nPr[im] += 1;

          // histos
          outerPID[iS][im]->Fill(cent, std::abs(trk_tmp->fPt), trk_tmp->fOuterPID);
        }
      }
      #ifdef FILL_MC
      for (int iM = 0; iM < 2; ++iM){
        hGenRecProton[iM][iVar - iVarMin]->Fill(cent, nPr_gen[iM], nPr[iM]);
        evtTupleGen[iS][iVar - iVarMin]->Fill(cent, qPr_1_gen_tmp[1], qPr_1_gen_tmp[0]);
      }
      #endif // FILL_MC

      evtTuple[iS][iVar - iVarMin]->Fill(cent, qPr_1_tmp[1], qPr_1_tmp[0], qPr_2_tmp[1], qPr_2_tmp[0], qPr_3_tmp[1], qPr_3_tmp[0], qPr_4_tmp[1], qPr_4_tmp[0], qPr_5_tmp[1], qPr_5_tmp[0], qPr_6_tmp[1], qPr_6_tmp[0]);
    }
  }

  // Process output
  for (int iVar{iVarMin}; iVar < iVarMax; ++iVar)
  {
    for (int iS{0}; iS < N_SAMPLE; ++iS){
      #ifdef FILL_MC
        o[iS][iVar - iVarMin]->cd();
        for (int iM = 0; iM < 2; ++iM){
          hGenProton[iM][iVar - iVarMin]->Write();
          hRecProton[iM][iVar - iVarMin]->Write();
          hGenRecProton[iM][iVar - iVarMin]->Write();
        }
        evtTupleGen[iS][iVar - iVarMin]->Write();
      #endif // FILL_MC

      o[iS][iVar - iVarMin]->cd();
      hCent[iS]->Write();
      hNtrkl[iS]->Write();
      for (int iC{0}; iC < 2; ++iC){
        outerPID[iS][iC]->Write();
      }
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
