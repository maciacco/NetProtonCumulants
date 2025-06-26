#include <TFile.h>
#include "utils.h"

void AOD2TTree(const char *inFile = "AO2D_20250618_smalleta_merged", const char *outFile = "newTree_20250618_smalleta", const bool isMC = false){
  TFile *fout = TFile::Open(Form("%s/%s.root", kResDir, outFile), "recreate");

  TTree *new_tree = new TTree("newtree", "new tree");
  TClonesArray *tracks = new TClonesArray("miniTrack");
  TClonesArray &tr = *tracks;
  char zvtxMask = 0;
  uint8_t triggerMask = 0u;
  uint8_t ntracklets = 0u;
  uint8_t v0Multiplicity = 0u;
  new_tree->Branch("fZvtxMask", &zvtxMask);
  new_tree->Branch("fTriggerMask", &triggerMask);
  new_tree->Branch("fNtracklets", &ntracklets);
  new_tree->Branch("fV0Multiplicity", &v0Multiplicity);
  new_tree->Branch("fTracks", &tracks);

  // track table (w/ MC info)
  float fPt = 0.f;
  char fEtaMask = 0;
  int fSelMask = 0;
  float fOuterPID = 0.f;
  float fGenPt = -999.f;
  char fGenEtaMask = 100;
  bool fIsReco = 0;
  int32_t fIndexMiniCollTables = 0;

  // collision table
  int64_t icoll = 0;
  char fZvtxMask = 0;
  uint8_t fTriggerMask = 0;
  uint8_t fNtracklets = 0u;
  uint8_t fV0Multiplicity = 0u;

  TFile *_file0 = TFile::Open(Form("%s/%s.root", kDataDir, inFile));
  auto _list_keys = _file0->GetListOfKeys();
  for (auto _key : *_list_keys) {
    const char* c_keyName = _key->GetName();
    string keyName(c_keyName);
    if (keyName.find("DF_") != 0)
      continue;

    TTree *_tree_coll = (TTree*)_file0->Get(Form("%s/O2minicolltable", c_keyName));
    TTree *_tree_trk = (TTree*)_file0->Get(Form("%s/O2%sminitrktable", c_keyName, isMC ? "mc" : ""));
    int64_t ncoll = _tree_coll->GetEntries();
    int64_t ntrk = _tree_trk->GetEntries();

    _tree_trk->SetBranchAddress("fIndexMiniCollTables", &fIndexMiniCollTables);
    _tree_trk->SetBranchAddress("fPt", &fPt);
    _tree_trk->SetBranchAddress("fEtaMask", &fEtaMask);
    _tree_trk->SetBranchAddress("fSelMask", &fSelMask);
    _tree_trk->SetBranchAddress("fOuterPID", &fOuterPID);
    if (isMC) {
      _tree_trk->SetBranchAddress("fGenPt", &fGenPt);
      _tree_trk->SetBranchAddress("fGenEtaMask", &fGenEtaMask);
      _tree_trk->SetBranchAddress("fIsReco", &fIsReco);
    }

    _tree_coll->SetBranchAddress("fZvtxMask", &fZvtxMask);
    _tree_coll->SetBranchAddress("fTriggerMask", &fTriggerMask);
    _tree_coll->SetBranchAddress("fNtracklets", &fNtracklets);
    _tree_coll->SetBranchAddress("fV0Multiplicity", &fV0Multiplicity);

    int64_t itrk = 0;
    for (int64_t ic = 0; ic < ncoll; ++ic) {
      if (!(ic % 1000)) std::cout << ic << std::endl;
      tracks->Clear();
      int it = 0;

      _tree_coll->GetEntry(ic);
      zvtxMask = fZvtxMask;
      triggerMask = fTriggerMask;
      ntracklets = fNtracklets;
      v0Multiplicity = fV0Multiplicity;

      // fill tracks array
      for (; itrk < ntrk; ++itrk) {
        _tree_trk->GetEntry(itrk);
        if (fIndexMiniCollTables == ic) {
          new (tr[it]) miniTrack();
          auto trk = dynamic_cast<miniTrack*>(tr[it++]);
          trk->fPt = fPt;
          trk->fEtaMask = fEtaMask;
          trk->fSelMask = fSelMask;
          trk->fOuterPID = fOuterPID;
          trk->fGenPt = fGenPt;
          trk->fGenEtaMask = fGenEtaMask;
          trk->fIsReco = fIsReco;
        }
        else if (fIndexMiniCollTables > ic) {
          itrk--;
          break;
        }
      }

      new_tree->Fill();
    }
  }

  fout->cd();
  new_tree->Write();
  fout->Close();
}
