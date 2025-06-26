#ifndef __UTILS_H__
#define __UTILS_H__

#include <TObject.h>
#include <TColor.h>
#include <TClonesArray.h>

enum kSources{
  kPrim = 0,
  kWD = 1
};

struct miniTrack : public TObject{
  miniTrack() :
  fPt{-999.f},
  fEtaMask{0},
  fSelMask{0},
  fOuterPID{-999.f},
  fGenPt{-999.f},
  fGenEtaMask{0},
  fIsReco{false}
  {};
  miniTrack(const float pt, const int8_t etaMask, const int selMask, const float outerPID, const float genPt, const int8_t genEtaMask, const bool isReco) :
  fPt{pt},
  fEtaMask{etaMask},
  fSelMask{selMask},
  fOuterPID{outerPID},
  fGenPt{genPt},
  fGenEtaMask{genEtaMask},
  fIsReco{isReco}
  {};
  miniTrack(miniTrack const& other) :
  fPt{other.fPt},
  fEtaMask{other.fEtaMask},
  fSelMask{other.fSelMask},
  fOuterPID{other.fOuterPID},
  fGenPt{other.fGenPt},
  fGenEtaMask{other.fGenEtaMask},
  fIsReco{other.fIsReco}
  {};
  float fPt;
  int8_t fEtaMask;
  int fSelMask;
  float fOuterPID;
  float fGenPt;
  int8_t fGenEtaMask;
  bool fIsReco;
};

struct miniEvent : public TObject{
  miniEvent() :
  fZvtxMask{0},
  fTriggerMask{0u},
  fNtracklets{0u},
  fV0Multiplicity{0u}
  {};
  miniEvent(char const zvtxMask, uint8_t const triggerMask, uint8_t const ntracklets, uint8_t const v0Multiplicity) :
  fZvtxMask{zvtxMask},
  fTriggerMask{triggerMask},
  fNtracklets{ntracklets},
  fV0Multiplicity{v0Multiplicity}
  {};
  miniEvent(char const zvtxMask, uint8_t const triggerMask, uint8_t const ntracklets, uint8_t const v0Multiplicity, TClonesArray* const tracks) :
  fZvtxMask{zvtxMask},
  fTriggerMask{triggerMask},
  fNtracklets{ntracklets},
  fV0Multiplicity{v0Multiplicity}
  {
    for (int itrk = 0; itrk < tracks->GetEntries(); ++itrk) {
      miniTrack* trk_tmp = static_cast<miniTrack*>(tracks->At(itrk));
      this->fTracks.emplace_back(*trk_tmp);
    }
  };
  char fZvtxMask;
  uint8_t fTriggerMask;
  uint8_t fNtracklets;
  uint8_t fV0Multiplicity;
  std::vector<miniTrack> fTracks;
};

constexpr int maxDequeSize = 10000;
constexpr int kMixThrZvtx = 20;
constexpr int kMixThrV0M =  2;

double powI(double a, int p){
  if (p > 1)
    return a * powI(a, p - 1);
  return a;
}

constexpr const char* kAntiMatterLabel[2] = {"A", "M"};
constexpr const char* kAntiMatterLabelML[2] = {"antimatter", "matter"};
constexpr const char* kCorrLabel[3] = {"Uncorr", "Corr", "Gen"};
constexpr const char* kPartLabel[2] = {"Pr"};
constexpr const char* kPartLabelExtend[1] = {"Proton"};


int colors[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"),
TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#595959"), TColor::GetColor("#0033a1")};

const char *kEffPrFile = "prEff";
constexpr float kDummyEffPr = 1.;
constexpr const char* kSubsampleFlag = "_";
const bool kMultV0M = true; //false;

constexpr int kLimitSample = 1;
constexpr int kLimitedSample = 3000000;
const char *kDataDir = "data";
const char *kResDir = "/home/mciacco/Code/NetProtonCumulants/results";
//const char *kResDir = "results";
const char *kCalibDir = "calib";
constexpr bool kUseIndex = false;
// constexpr bool isMC = false;

constexpr int kNSample = 20;

// constexpr int kNCentBins = 2;
// constexpr double kCentBins[kNCentBins + 1] = {0.f, 0.1f, 100.f};
constexpr int kNCentBins = 7;
constexpr double kCentBins[kNCentBins + 1] = {0, 10, 20, 30, 40, 50, 70, 100};

//constexpr int kNCentBins = 100;
//constexpr double kCentBins[kNCentBins + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
//26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
// 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92,
//93, 94, 95, 96, 97, 98, 99, 100};

constexpr int kNVtxBins = 8;
constexpr double kVtxBins[kNVtxBins + 1]{-100., -58., -35, -16., 0., 16., 35., 58., 100.};

constexpr int kNTrklBins = 11;
constexpr double kTrklBins[kNTrklBins + 1] = /* {0, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f, 24.f, 25.f, 26.f, 27.f, 28.f, 29.f, 30.f, 35.f, 40.f}; */ {0, 3, 4, 5, 7, 9, 11, 15, 19, 27, 50, 200};

constexpr uint8_t kTriggerSel = 0x1; //0x1;
constexpr int kNEtaBins = 1;
constexpr int kNBinsPt = 20;
constexpr int kNBinsPID = 50;
// constexpr int kNVtxBins = 8;
constexpr int kMultBins = 30;
constexpr double kMinEta = -0.8f;
constexpr double kMinMult = 0.f;
constexpr double kDeltaEta = 1.6f;
constexpr double kMinPt = 0.5f;
// constexpr double kMinVtx = -100.f;
constexpr double kDeltaPt = 0.05f;
constexpr double kMinPID = -5.f;
constexpr double kDeltaPID = 0.175f;
// constexpr double kDeltaVtx = 25.f;
constexpr double kDeltaMult = 4.f;

constexpr int kNTPCcls = 3;
constexpr int kNChi2TPC = 3;
constexpr int kNDCAxy = 3;
constexpr int kNDCAz = 3;
constexpr int kNITSPID = 3;
constexpr int kNTPCPID = 3;

constexpr double kMaxCent = 100;
constexpr double kMinCent = 0;

constexpr int8_t kEtaCut = 80;
constexpr int8_t kEtaCutMin = 0;
constexpr float kPtLowLimitPr = 0.5f;
constexpr float kTOFptCut = 1.5f;
constexpr int8_t kZvtxCut = 100;
constexpr int8_t kZvtxCutMin = 0;

// variations
// tpc cls
constexpr int kCutTPCcls[] = {BIT(0), BIT(0), BIT(0)};
constexpr int kCutTPCcls2[] = {BIT(0), BIT(1), BIT(1)};
constexpr bool kRequireTPCclsCut[] = {true, true, false};

// chi2 tpc
constexpr int kCutChi2TPC[] = {BIT(2), BIT(2), BIT(2)};
constexpr int kCutChi2TPC2[] = {BIT(2), BIT(3), BIT(3)};
constexpr bool kRequireChi2TPCCut[] = {true, true, false};

// dcaxy
constexpr int kCutDCAxy[] = {BIT(4), BIT(4), BIT(4)};
constexpr int kCutDCAxy2[] = {BIT(4), BIT(5), BIT(5)};
constexpr bool kRequireDCAxyCut[] = {true, true, false};

// dcaz
constexpr int kCutDCAz[] = {BIT(6), BIT(6), BIT(6)};
constexpr int kCutDCAz2[] = {BIT(6), BIT(7), BIT(7)};
constexpr bool kRequireDCAzCut[] = {true, true, false};

// its pid
constexpr int kCutITSPID[] = {BIT(8), BIT(8), BIT(8)};
constexpr int kCutITSPID2[] = {BIT(8), BIT(9), BIT(9)};
constexpr bool kRequireITSPIDCut[] = {true, true, false};

// tpc pid
constexpr int kCutTPCPID[] = {BIT(10), BIT(10), BIT(10)};
constexpr int kCutTPCPID2[] = {BIT(10), BIT(11), BIT(11)};
constexpr bool kRequireTPCPIDCut[] = {true, true, false};

// constexpr int kNCentBinsSmall = 32;
// constexpr float kCentBinsSmall[kNCentBinsSmall + 1] = {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f, 24.f, 25.f, 26.f, 27.f, 28.f, 29.f, 30.f, 35.f, 40.f};// {0.f, 3.f, 4.f, 5.f, 7.f, 9.f, 11.f, 15.f, 19.f, 27.f, 100.f};

//constexpr int kNCentBinsSmall = 2;
//constexpr double kCentBinsSmall[kNCentBinsSmall + 1] = {0.f, 0.1f, 100.f};
constexpr int kNCentBinsSmall = 7;
constexpr double kCentBinsSmall[kNCentBinsSmall + 1] = {0, 10, 20, 30, 40, 50, 70, 100};

//constexpr int kNCentBinsSmall = 100;
//constexpr double kCentBinsSmall[kNCentBinsSmall + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
//26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
// 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92,
//93, 94, 95, 96, 97, 98, 99, 100};

constexpr int kNCentBinsSmallTmp = 100;
constexpr double kCentBinsSmallTmp[kNCentBinsSmallTmp + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92,
93, 94, 95, 96, 97, 98, 99, 100};

constexpr int kNCentBinsMidTmp = 20;
constexpr double kCentBinsMidTmp[kNCentBinsMidTmp + 1] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.};

//constexpr int kNCentBinsSmall = 50;
//constexpr float kCentBinsSmall[kNCentBinsSmall + 1] = {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f, 24.f, 25.f, 26.f, 27.f, 28.f, 29.f, 30.f, 31.f, 32.f, 33.f, 34.f, 35.f, 36.f, 37.f, 38.f, 39.f, 40.f, 41.f, 42.f, 43.f, 44.f, 45.f, 46.f, 47.f, 48.f, 49.f, 50.f};
constexpr int kNTrklBinsSmall = 51;
constexpr double kTrklBinsSmall[kNTrklBinsSmall + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 200};

// systematics vars
constexpr int kNVar = 42;
constexpr int kVar[kNVar + 1] = {364, 81, 87, 93, 99, 105, 111, 117, 123, 129, 135, 141, 147, 153, 159, 327, 333, 339, 345, 351, 357, 363, 369, 375, 381, 387, 393, 399, 567, 573, 579, 585, 591, 597, 603, 609, 615, 621, 627, 633, 639, 645};

#endif
