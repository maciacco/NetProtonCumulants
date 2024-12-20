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
#include <TNtupleD.h>
#include <TStopwatch.h>
#include <TClonesArray.h>

void ReadTree(const char* fname = "newTree_noTOF", const char* ofname = "LHC18pp_noTOF", const int iVarMin = 526, const int iVarMax = 527, const bool isMC = false)
{
  TStopwatch w;
  w.Start();

  int nSample = isMC ? 1 : kNSample;

  TFile *fEffPr = isMC ? nullptr : TFile::Open(Form("%s/%s.root", kResDir, kEffPrFile));
  TFile f(Form("%s/%s.root", kResDir, fname));

  const int nF = iVarMax - iVarMin;
  TFile ***o = new TFile**[nSample];
  for (int i{0}; i < nSample; ++i) {
    o[i] = new TFile*[nF];
  }
  for (int iS{0}; iS < nSample; ++iS){
    for (int iVar{iVarMin}; iVar < iVarMax; ++iVar){
      o[iS][iVar - iVarMin] = new TFile(Form("%s/%s%d_var_%d.root", kResDir, ofname, iS, iVar), "recreate");
    }
  }

  // track table (w/ MC info)
  float fPt = 0.f;
  int8_t fEtaMask = 0;
  int fSelMask = 0;
  float fOuterPID = 0.f;
  float fGenPt = -999.f;
  int8_t fGenEtaMask = 100;
  bool fIsReco = 0;
  int32_t fIndexMiniCollTables = 0;

  // collision table
  int64_t icoll = 0;
  char fZvtxMask = 0;
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
  TH1D **hCent = new TH1D*[nSample];
  TH1D **hNtrkl = new TH1D*[nSample];

  TH3D ***outerPID = new TH3D**[nSample];
  TNtupleD ***evtTuple = new TNtupleD**[nSample];
  TNtupleD ***evtTupleGen = new TNtupleD**[nSample];
  for (int i{0}; i < nSample; ++i) {
    outerPID[i] = new TH3D*[2];
    evtTuple[i] = new TNtupleD*[nF];
    evtTupleGen[i] = new TNtupleD*[nF];
  }

  TH1D ******hEffPr = new TH1D*****[2];
  for (int i{0}; i < 2; ++i) {
    hEffPr[i] = new TH1D****[kNCentBins];
    for (int j{0}; j < kNCentBins; ++j) {
      hEffPr[i][j] = new TH1D***[kNEtaBins];
      for (int k{0}; k < kNEtaBins; ++k) {
        hEffPr[i][j][k] = new TH1D**[nSample];
        for (int l{0}; l < nSample; ++l) {
          hEffPr[i][j][k][l] = new TH1D*[nF];
        }
      }
    }
  }

  TH1D ******hEffPrTrkl = new TH1D*****[2];
  for (int i{0}; i < 2; ++i) {
    hEffPrTrkl[i] = new TH1D****[kNTrklBins];
    for (int j{0}; j < kNTrklBins; ++j) {
      hEffPrTrkl[i][j] = new TH1D***[kNEtaBins];
      for (int k{0}; k < kNEtaBins; ++k) {
        hEffPrTrkl[i][j][k] = new TH1D**[nSample];
        for (int l{0}; l < nSample; ++l) {
          hEffPrTrkl[i][j][k][l] = new TH1D*[nF];
        }
      }
    }
  }

  for (int iS = 0; iS < nSample; ++iS){
    hCent[iS] = new TH1D(Form("hCent_%d", iS), ";Centrality (%);Entries", kNCentBinsSmall, kCentBinsSmall);
    hNtrkl[iS] = new TH1D(Form("hNtrkl_%d", iS), ";#it{N}_{tracklets}^{0.7 < |#eta| < 1.2};Entries", kNTrklBinsSmall, kTrklBinsSmall);

    for (int iVar{iVarMin}; iVar < iVarMax; ++iVar)
    {
      evtTuple[iS][iVar - iVarMin] = new TNtupleD(Form("evtTuple_%d", iVar), Form("evtTuple_%d", iVar), "cent:q1pP:q1pN:q2pP:q2pN:q3pP:q3pN:q4pP:q4pN:q5pP:q5pN:q6pP:q6pN:ntrkl");
      evtTuple[iS][iVar - iVarMin]->SetDirectory(o[iS][iVar - iVarMin]);
      for (int iC = 0; iC < 2; ++iC){
        if (isMC) {
          evtTupleGen[iS][iVar - iVarMin] = new TNtupleD(Form("evtTupleGen_%d", iVar), Form("evtTupleGen_%d", iVar), "cent:q1pP:q1pN:ntrkl");
          evtTupleGen[iS][iVar - iVarMin]->SetDirectory(o[iS][iVar - iVarMin]);
        }
        for (int iEta = 0; iEta < kNEtaBins; ++iEta){
          if (fEffPr){
            for (int iCent = 0; iCent < kNCentBins; ++iCent){
              hEffPr[iC][iCent][iEta][iS][iVar - iVarMin] = (TH1D*)fEffPr->Get(Form("subsample_%d_var_%d/h%sEff%s_%d_%d_%d", 1, iVar, kAntiMatterLabel[iC], kPartLabel[0], iCent, iEta, iVar));
            }
            for (int iTrkl = 0; iTrkl < kNTrklBins; ++iTrkl){
              hEffPrTrkl[iC][iTrkl][iEta][iS][iVar - iVarMin] = (TH1D*)fEffPr->Get(Form("subsample_%d_var_%d/h%sEffTrkl%s_%d_%d_%d", 1, iVar, kAntiMatterLabel[iC], kPartLabel[0], iTrkl, iEta, iVar));
            }
          }
        }
      }
    }
  }

  // histos
  float etaBins[kNEtaBins + 1];
  for (int iB = 0; iB < kNEtaBins + 1; ++iB){
    etaBins[iB] = kMinEta + kDeltaEta * iB;
  }
  float ptBins[kNBinsPt + 1];
  for (int iB = 0; iB < kNBinsPt + 1; ++iB){
    ptBins[iB] = kMinPt + kDeltaPt * iB;
  }
  float pidBins[kNBinsPID + 1];
  for (int iB = 0; iB < kNBinsPID + 1; ++iB){
    pidBins[iB] = kMinPID + kDeltaPID * iB;
  }
  for (int iS = 0; iS < nSample; ++iS){
    for (int iC = 0; iC < 2; ++iC){
      outerPID[iS][iC] = new TH3D(Form("h%sOuterPID_%d", kAntiMatterLabel[iC], iS), ";Centrality (%);#it{p}_{T} (GeV/#it{c});n#sigma (a.u.)", kNCentBins, kCentBins, kNBinsPt, ptBins, kNBinsPID, pidBins);
    }
  }

  TH3F *hGenRecProton[2][nF];
  TH3D *hGenProton[2][nF];
  TH3D *hRecProton[2][nF];
  TH3D *hGenProtonTrkl[2][nF];
  TH3D *hRecProtonTrkl[2][nF];

  //if (isMC) {
    for (int iV{iVarMin}; iV < iVarMax; ++iV){
      for (int iC = 0; iC < 2; ++iC){
        hRecProton[iC][iV - iVarMin] = new TH3D(Form("h%sRecProton_%d", kAntiMatterLabel[iC], iV), ";Centrality (%);#it{p}_{T} (GeV/#it{c});#eta", kNCentBins, kCentBins, kNBinsPt, ptBins, kNEtaBins, etaBins);
        hRecProtonTrkl[iC][iV - iVarMin] = new TH3D(Form("h%sRecProtonTrkl_%d", kAntiMatterLabel[iC], iV), ";#it{N}_{trkl};#it{p}_{T} (GeV/#it{c});#eta", kNTrklBins, kTrklBins, kNBinsPt, ptBins, kNEtaBins, etaBins);
        if (isMC) {
          hGenProton[iC][iV - iVarMin] = new TH3D(Form("h%sGenProton_%d", kAntiMatterLabel[iC], iV), ";Centrality (%);#it{p}_{T} (GeV/#it{c});#eta", kNCentBins, kCentBins, kNBinsPt, ptBins, kNEtaBins, etaBins);
          hGenProtonTrkl[iC][iV - iVarMin] = new TH3D(Form("h%sGenProtonTrkl_%d", kAntiMatterLabel[iC], iV), ";#it{N}_{trkl};#it{p}_{T} (GeV/#it{c});#eta", kNTrklBins, kTrklBins, kNBinsPt, ptBins, kNEtaBins, etaBins);
          hGenRecProton[iC][iV - iVarMin] = new TH3F(Form("h%sGenRecProton_%d", kAntiMatterLabel[iC], iV), ";Centrality (%);#it{N}_{gen};#it{N_rec}", kNCentBins, 0, 100, 200, 0, 200, 200, 0, 200);
        }
      }
    }
  //}

  // loop variables
  Long64_t nEntries = kLimitSample ? kLimitedSample : t->GetEntries();
  TH1D hCentTmp("hCentTmp", "hCentTmp", kNCentBins, kCentBins);
  TH1D hCentSmallTmp("hCentSmallTmp", "hCentSmallTmp", kNCentBinsSmall, kCentBinsSmall);
  TH1D hEtaTmp("hEtaTmp", "hEtaTmp", kNEtaBins, etaBins);

  // Event loop
  gRandom->SetSeed(42);
  for (Long64_t i = 0; i < nEntries; ++i){
    const int iS = (int)(gRandom->Rndm() * nSample);

    Long64_t e = i;
    if (!(i%1000000)) std::cout << "n_ev = " << i << std::endl;

    Long64_t tentry = t->LoadTree(e);
    t->GetEntry(tentry);

    float cent = fV0Multiplicity;

    if ((fTriggerMask & 0x2) == 0x2) { // high granularity
      cent = cent / 100.f;
    }

    if (!((fTriggerMask & kTriggerSel) == kTriggerSel) && !isMC) {
      continue;
    }

    if (cent > kMaxCent) continue;
    int ic = hCentTmp.FindBin(cent);
    int ic_sm = hCentSmallTmp.FindBin(cent);
    if (std::abs(fZvtxMask) > kZvtxCut || std::abs(fZvtxMask) < kZvtxCutMin) continue;

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

      double qPr_1_gen_tmp[] = {0, 0};
      Long64_t nPr_gen[] = {0, 0};

      double qPr_1_tmp[2] = {0, 0};
      double qPr_2_tmp[2] = {0, 0};
      double qPr_3_tmp[2] = {0, 0};
      double qPr_4_tmp[2] = {0, 0};
      double qPr_5_tmp[2] = {0, 0};
      double qPr_6_tmp[2] = {0, 0};
      Long64_t nPr[] = {0, 0};

      for (int itrk = 0; itrk < tracks->GetEntries(); ++itrk) {

        miniTrack* trk_tmp = (miniTrack*)tracks->At(itrk);

        if (isMC) {
          if ( std::abs(trk_tmp->fGenPt) > (kTOFptCut + 1.) || std::abs(trk_tmp->fGenPt) < (kPtLowLimitPr - .1) || trk_tmp->fGenEtaMask > (kEtaCut + 10.) || trk_tmp->fGenEtaMask < -(kEtaCut + 10.) ) continue;
          if ( trk_tmp->fGenEtaMask > -kEtaCutMin && trk_tmp->fGenEtaMask < kEtaCutMin ) continue;
          float eta_MC = static_cast<float>(trk_tmp->fGenEtaMask) / 100.f;
          int ie_MC = hEtaTmp.FindBin(eta_MC);
          int im_MC = trk_tmp->fGenPt > 0 ? 1 : 0;
          hGenProton[im_MC][iVar - iVarMin]->Fill(cent, std::abs(trk_tmp->fGenPt), eta_MC);
          hGenProtonTrkl[im_MC][iVar - iVarMin]->Fill(fNtracklets, std::abs(trk_tmp->fGenPt), eta_MC);
          if ( std::abs(trk_tmp->fGenPt) < kTOFptCut && std::abs(trk_tmp->fGenPt) > kPtLowLimitPr ) {
            // std::cout << im_MC << std::endl;
            qPr_1_gen_tmp[im_MC] += 1.;
            nPr_gen[im_MC] += 1;
          }
        }
        if (
            ( ( (((trk_tmp->fSelMask & kCutTPCcls[iTPCcls]) == kCutTPCcls[iTPCcls]) || ((trk_tmp->fSelMask & kCutTPCcls2[iTPCcls]) == kCutTPCcls2[iTPCcls]) ) && kRequireTPCclsCut[iTPCcls] ) || !kRequireTPCclsCut[iTPCcls] ) &&
            ( ( (((trk_tmp->fSelMask & kCutChi2TPC[iChi2TPC]) == kCutChi2TPC[iChi2TPC]) || ((trk_tmp->fSelMask & kCutChi2TPC2[iChi2TPC]) == kCutChi2TPC2[iChi2TPC]) ) && kRequireChi2TPCCut[iChi2TPC] ) || !kRequireChi2TPCCut[iChi2TPC] ) &&
            ( ( ((trk_tmp->fSelMask & kCutDCAxy[iDCAxy]) == kCutDCAxy[iDCAxy] || (trk_tmp->fSelMask & kCutDCAxy2[iDCAxy]) == kCutDCAxy2[iDCAxy]) && kRequireDCAxyCut[iDCAxy] ) || !kRequireDCAxyCut[iDCAxy] ) &&
            ( ( ((trk_tmp->fSelMask & kCutDCAz[iDCAz]) == kCutDCAz[iDCAz] || (trk_tmp->fSelMask & kCutDCAz2[iDCAz]) == kCutDCAz2[iDCAz]) && kRequireDCAzCut[iDCAz] ) || !kRequireDCAzCut[iDCAz] ) &&
            std::abs(trk_tmp->fPt) > kPtLowLimitPr && std::abs(trk_tmp->fPt) < kTOFptCut &&
            (trk_tmp->fEtaMask < kEtaCut && trk_tmp->fEtaMask > -kEtaCut) &&
            (trk_tmp->fEtaMask < -kEtaCutMin || trk_tmp->fEtaMask > kEtaCutMin) &&
            (std::abs(trk_tmp->fOuterPID) < 2.f)
          )
        {
          if (isMC) {
            if (!trk_tmp->fIsReco || trk_tmp->fGenPt < -998.f) continue;
          }
          int im = trk_tmp->fPt > 0 ? 1 : 0;
          float eta_ = static_cast<float>(trk_tmp->fEtaMask) / 100.f;
          int ie = hEtaTmp.FindBin(eta_);
          double eff = fEffPr ? hEffPr[im][ic - 1][ie - 1][iS][iVar - iVarMin]->GetBinContent(hEffPr[im][ic - 1][ie - 1][iS][iVar - iVarMin]->FindBin(std::abs(trk_tmp->fPt))) : kDummyEffPr;
          if (fEffPr && eff < 1.e-12) continue;
          // if (im > 0.5 && !isMC) {eff *= 0.96;} // test systematic effect on k1
          hRecProton[im][iVar - iVarMin]->Fill(cent, std::abs(trk_tmp->fPt), eta_, 1./eff);
          hRecProtonTrkl[im][iVar - iVarMin]->Fill(fNtracklets, std::abs(trk_tmp->fPt), eta_, 1./eff);
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
      if (isMC) {
        for (int iM = 0; iM < 2; ++iM){
          hGenRecProton[iM][iVar - iVarMin]->Fill(cent, nPr_gen[iM], nPr[iM]);
          hGenRecProton[iM][iVar - iVarMin]->Fill(cent, nPr_gen[iM], nPr[iM]);
          evtTupleGen[iS][iVar - iVarMin]->Fill(cent, qPr_1_gen_tmp[1], qPr_1_gen_tmp[0], fNtracklets);
        }
      }

      evtTuple[iS][iVar - iVarMin]->Fill(cent, qPr_1_tmp[1], qPr_1_tmp[0], qPr_2_tmp[1], qPr_2_tmp[0], qPr_3_tmp[1], qPr_3_tmp[0], qPr_4_tmp[1], qPr_4_tmp[0], qPr_5_tmp[1], qPr_5_tmp[0], qPr_6_tmp[1], qPr_6_tmp[0], fNtracklets);
    }
  }

  // Process output
  for (int iVar{iVarMin}; iVar < iVarMax; ++iVar)
  {
    for (int iS{0}; iS < nSample; ++iS){
      //if (isMC) {
        o[iS][iVar - iVarMin]->cd();
        for (int iM = 0; iM < 2; ++iM){
          if (isMC) {
            hGenProton[iM][iVar - iVarMin]->Write();
            hGenProtonTrkl[iM][iVar - iVarMin]->Write();
          }
          hRecProton[iM][iVar - iVarMin]->Write();
          hRecProtonTrkl[iM][iVar - iVarMin]->Write();
          if (isMC) hGenRecProton[iM][iVar - iVarMin]->Write();
        }
        if (isMC) evtTupleGen[iS][iVar - iVarMin]->Write();
      //}

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
    for (int iS{0}; iS < nSample; ++iS){
      o[iS][iVar - iVarMin]->Close();
    }
  }

  w.Stop();
  w.Print();

  f.Close();
}
