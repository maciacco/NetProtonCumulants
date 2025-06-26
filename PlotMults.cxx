#include <TFile.h>
#include <ROOT/RDataFrame.hxx>

void PlotMults(const char* fname = "results/newTree_LHC18pp_20250611.root", const char* tname = "newtree"){
  ROOT::EnableImplicitMT(4);
  ROOT::RDataFrame df(tname, fname);
  TFile *out = TFile::Open("mulPlotsOut.root", "recreate");
  auto h2 = df.Filter("std::abs(fZvtxMask) < 100 && (fTriggerMask & 0x1) == 0x1 && fV0Multiplicity < 100")
              .Histo2D({"hTracksVsMult", ";Multiplicity (%);#it{N}_{Tracks}", 100, 0, 100, 200, 0, 200}, "fV0Multiplicity", "fNtracklets");
  out->cd();
  h2->Write();
  out->Close();
}