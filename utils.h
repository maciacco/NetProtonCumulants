#ifndef __UTILS_H__
#define __UTILS_H__

#include <TObject.h>
#include <TColor.h>

struct miniTrack : public TObject{
  float fPt;
  uint8_t fEtaMask;
  int fSelMask;
  float fOuterPID;
  float fGenPt;
  uint8_t fGenEtaMask;
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

const char *kEffPrFile = "results/prEff";
constexpr float kDummyEffPr = 1.;
constexpr const char* kSubsampleFlag = "_";

constexpr int kLimitSample = false;
constexpr int kLimitedSample = 100;
const char *kDataDir = "data";
const char *kResDir = "results";
const char *kCalibDir = "calib";
constexpr bool kUseIndex = true;
constexpr bool isMC = true;

constexpr int N_SAMPLE = 7;
constexpr int kNCentBins = 10;
constexpr float kCentBins[kNCentBins + 1] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};

constexpr int kNEtaBins = 1;
constexpr int kNBinsPt = 20;
constexpr int kNBinsPID = 32;
constexpr float kMinEta = -0.8f;
constexpr float kDeltaEta = 1.6f;
constexpr float kMinPt = 0.2f;
constexpr float kDeltaPt = 0.05f;
constexpr float kMinPID = -4.f;
constexpr float kDeltaPID = 0.25f;

constexpr int kNTPCcls = 3;
constexpr int kNChi2TPC = 3;
constexpr int kNDCAxy = 3;
constexpr int kNDCAz = 3;
constexpr int kNITSPID = 3;
constexpr int kNTPCPID = 3;

constexpr int kMaxCent = 100;

constexpr uint8_t kEtaCut = 8;
constexpr float kPtLowLimitPr = 0.4f;
constexpr float kTOFptCut = 1.2f;

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

//constexpr int kNCentBinsSmall = 10;
//constexpr float kCentBinsSmall[kNCentBinsSmall + 1] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
constexpr int kNCentBinsSmall = 100;
constexpr float kCentBinsSmall[kNCentBinsSmall + 1] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80., 81., 82., 83., 84., 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 96., 97., 98., 99., 100.};

#endif
