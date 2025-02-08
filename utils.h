#ifndef __UTILS_H__
#define __UTILS_H__

#include <TObject.h>
#include <TColor.h>

struct miniTrack : public TObject{
  float fPt;
  int8_t fEtaMask;
  int fSelMask;
  float fOuterPID;
  float fGenPt;
  int8_t fGenEtaMask;
  bool fIsReco;
};

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


int colors[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#595959"), TColor::GetColor("#0033a1")};

const char *kEffPrFile = "prEff";
constexpr float kDummyEffPr = 1.;
constexpr const char* kSubsampleFlag = "_";
const bool kMultV0M = true; //false;

constexpr int kLimitSample = 0;
constexpr int kLimitedSample = 30000000;
const char *kDataDir = "data";
const char *kResDir = "/home/mciacco/Code/NetProtonCumulants/results";
//const char *kResDir = "results";
const char *kCalibDir = "calib";
constexpr bool kUseIndex = false;
// constexpr bool isMC = false;

constexpr int kNSample = 20;

//constexpr int kNCentBins = 2;
//constexpr float kCentBins[kNCentBins + 1] = {0.f, 0.1f, 100.f};
constexpr int kNCentBins = 7;
constexpr float kCentBins[kNCentBins + 1] = {0.f, 10.f, 20.f, 30.f, 40.f, 50.f, 70.f, 100.f};

constexpr int kNTrklBins = 10;
constexpr float kTrklBins[kNTrklBins + 1] = /* {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f, 24.f, 25.f, 26.f, 27.f, 28.f, 29.f, 30.f, 35.f, 40.f}; */ {0.f, 3.f, 4.f, 5.f, 7.f, 9.f, 11.f, 15.f, 19.f, 27.f, 50.f};

constexpr uint8_t kTriggerSel = 0x1; //0x1;
constexpr int kNEtaBins = 1;
constexpr int kNBinsPt = 20;
constexpr int kNBinsPID = 50;
constexpr float kMinEta = -0.8f;
constexpr float kDeltaEta = 1.6f;
constexpr float kMinPt = 0.5f;
constexpr float kDeltaPt = 0.05f;
constexpr float kMinPID = -5.f;
constexpr float kDeltaPID = 0.175f;

constexpr int kNTPCcls = 3;
constexpr int kNChi2TPC = 3;
constexpr int kNDCAxy = 3;
constexpr int kNDCAz = 3;
constexpr int kNITSPID = 3;
constexpr int kNTPCPID = 3;

constexpr int kMaxCent = 100;

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
//constexpr float kCentBinsSmall[kNCentBinsSmall + 1] = {0.f, 0.1f, 100.f};
constexpr int kNCentBinsSmall = 7;
constexpr float kCentBinsSmall[kNCentBinsSmall + 1] = {0.f, 10.f, 20.f, 30.f, 40.f, 50.f, 70.f, 100.f};
//constexpr int kNCentBinsSmall = 100;
//constexpr float kCentBinsSmall[kNCentBinsSmall + 1] = {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f, 24.f, 25.f, 26.f, 27.f, 28.f, 29.f, 30.f, 31.f, 32.f, 33.f, 34.f, 35.f, 36.f, 37.f, 38.f, 39.f, 40.f, 41.f, 42.f, 43.f, 44.f, 45.f, 46.f, 47.f, 48.f, 49.f, 50.f, 51.f, 52.f, 53.f, 54.f, 55.f, 56.f, 57.f, 58.f, 59.f, 60.f, 61.f, 62.f, 63.f, 64.f, 65.f, 66.f, 67.f, 68.f, 69.f, 70.f, 71.f, 72.f, 73.f, 74.f, 75.f, 76.f, 77.f, 78.f, 79.f, 80.f, 81.f, 82.f, 83.f, 84.f, 85.f, 86.f, 87.f, 88.f, 89.f, 90.f, 91.f, 92.f, 93.f, 94.f, 95.f, 96.f, 97.f, 98.f, 99.f, 100.f};
//constexpr int kNCentBinsSmall = 50;
//constexpr float kCentBinsSmall[kNCentBinsSmall + 1] = {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f, 24.f, 25.f, 26.f, 27.f, 28.f, 29.f, 30.f, 31.f, 32.f, 33.f, 34.f, 35.f, 36.f, 37.f, 38.f, 39.f, 40.f, 41.f, 42.f, 43.f, 44.f, 45.f, 46.f, 47.f, 48.f, 49.f, 50.f};
constexpr int kNTrklBinsSmall = 50;
constexpr float kTrklBinsSmall[kNTrklBinsSmall + 1] = {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f, 24.f, 25.f, 26.f, 27.f, 28.f, 29.f, 30.f, 31.f, 32.f, 33.f, 34.f, 35.f, 36.f, 37.f, 38.f, 39.f, 40.f, 41.f, 42.f, 43.f, 44.f, 45.f, 46.f, 47.f, 48.f, 49.f, 50.f};

// systematics vars
constexpr int kNVar = 42;
constexpr int kVar[kNVar + 1] = {364, 81, 87, 93, 99, 105, 111, 117, 123, 129, 135, 141, 147, 153, 159, 327, 333, 339, 345, 351, 357, 363, 369, 375, 381, 387, 393, 399, 567, 573, 579, 585, 591, 597, 603, 609, 615, 621, 627, 633, 639, 645};

#endif
