#include "utils.h"

#include <TStopwatch.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNtupleD.h>
#include <Riostream.h>
#include <TProfile.h>

void ProcessTupleSingleParticleHighOrder(const int smpl = 0, const int iVarMin = 364, const int iVarMax = 365, const bool FILL_MC = false)
{
  TStopwatch w;
  w.Start();

  int sample = smpl;
  int skippedVar = 0;

  for (int iVar = iVarMin; iVar < iVarMax; ++iVar){
    bool inVars = false;
    for (int iV{0}; iV < kNVar; ++iV) {
      if (iVar == kVar[iV]) {
        inVars = true;
        break;
      }
    }
    if (!inVars) continue;
    const char* dir = kTriggerSel == 0x1 ? "MB" : "HM";
    TFile *fin = TFile::Open(Form("%s/%s/LHC18ppTrig_HM%d_var_%d.root", kResDir, dir, sample, iVar));
    if (!fin || fin->TestBit(TFile::kZombie)){
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }
    TFile fout(Form("%s/%s/output_sys_singleParticle_%d_%d.root", kResDir, dir, sample, iVar), "recreate");

    TNtupleD *tuple_qmoment = (TNtupleD*)fin->Get(Form("evtTuple_%d", iVar));
    TH1D *hCent = (TH1D*)fin->Get(Form("hCent_%d", smpl));

    if (!tuple_qmoment)
    {
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }
    tuple_qmoment->SetCacheSize(0);

    TNtupleD *tuple_moment_gen;
    if (FILL_MC) {
      tuple_moment_gen = (TNtupleD*)fin->Get(Form("evtTupleGen_%d", iVar));
      tuple_moment_gen->SetCacheSize(0);
    }

    int evt[10] = {0};
    int centrality;

    double *arg;
    double total_event = tuple_qmoment->GetEntriesFast();

    double *arg_gen;
    double total_event_gen = 0.;
    TProfile *N1p = nullptr;
    TProfile *N1 = nullptr;
    TProfile *N2 = nullptr;
    TProfile *N3 = nullptr;
    TProfile *N4 = nullptr;
    TProfile *N5 = nullptr;
    TProfile *N6 = nullptr;
    if (FILL_MC) {
      total_event_gen = tuple_moment_gen->GetEntriesFast();
      // non-central moments of generated distribution
      N1p = new TProfile("N1p", "N1p", kNCentBinsSmall, kCentBinsSmall);
      N1 = new TProfile("N1", "N1", kNCentBinsSmall, kCentBinsSmall);
      N2 = new TProfile("N2", "N2", kNCentBinsSmall, kCentBinsSmall);
      N3 = new TProfile("N3", "N3", kNCentBinsSmall, kCentBinsSmall);
      N4 = new TProfile("N4", "N4", kNCentBinsSmall, kCentBinsSmall);
      N5 = new TProfile("N5", "N5", kNCentBinsSmall, kCentBinsSmall);
      N6 = new TProfile("N6", "N6", kNCentBinsSmall, kCentBinsSmall);
    }

    // full formula
    // 1st order
    TProfile* q1_1_1 = new TProfile("q1_1_1", "q1_1_1", kNCentBinsSmall, kCentBinsSmall);

    // 2nd order
    TProfile* q2_1_1 = new TProfile("q2_1_1", "q2_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_1 = new TProfile("q1_2_1", "q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_2 = new TProfile("q1_2_2", "q1_2_2", kNCentBinsSmall, kCentBinsSmall);

    // 3rd order
    TProfile* q3_1_1 = new TProfile("q3_1_1", "q3_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_2_1 = new TProfile("q1_1_1_x_q1_2_1", "q1_1_1_x_q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_2_2 = new TProfile("q1_1_1_x_q1_2_2", "q1_1_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_3_1 = new TProfile("q1_3_1", "q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_3_2 = new TProfile("q1_3_2", "q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_3_3 = new TProfile("q1_3_3", "q1_3_3", kNCentBinsSmall, kCentBinsSmall);

    // // 4th order
    TProfile* q4_1_1 = new TProfile("q4_1_1", "q4_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_2_1 = new TProfile("q2_1_1_x_q1_2_1", "q2_1_1_x_q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_2_2 = new TProfile("q2_1_1_x_q1_2_2", "q2_1_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_3_1 = new TProfile("q1_1_1_x_q1_3_1", "q1_1_1_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_2_1 = new TProfile("q2_2_1", "q2_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_2_2 = new TProfile("q2_2_2", "q2_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_3_2 = new TProfile("q1_1_1_x_q1_3_2", "q1_1_1_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_3_3 = new TProfile("q1_1_1_x_q1_3_3", "q1_1_1_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_1_x_q1_2_2 = new TProfile("q1_2_1_x_q1_2_2", "q1_2_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_4_1 = new TProfile("q1_4_1", "q1_4_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_4_2 = new TProfile("q1_4_2", "q1_4_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_4_3 = new TProfile("q1_4_3", "q1_4_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_4_4 = new TProfile("q1_4_4", "q1_4_4", kNCentBinsSmall, kCentBinsSmall);

    // 5th order
    TProfile* q5_1_1 = new TProfile("q5_1_1", "q5_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q3_1_1_x_q1_2_1 = new TProfile("q3_1_1_x_q1_2_1", "q3_1_1_x_q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q3_1_1_x_q1_2_2 = new TProfile("q3_1_1_x_q1_2_2", "q3_1_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_3_1 = new TProfile("q2_1_1_x_q1_3_1", "q2_1_1_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_3_2 = new TProfile("q2_1_1_x_q1_3_2", "q2_1_1_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_3_3 = new TProfile("q2_1_1_x_q1_3_3", "q2_1_1_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_2_2_x_q1_1_1 = new TProfile("q2_2_2_x_q1_1_1", "q2_2_2_x_q1_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_2_1_x_q1_1_1 = new TProfile("q2_2_1_x_q1_1_1", "q2_2_1_x_q1_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_2_1_x_q1_2_2 = new TProfile("q1_1_1_x_q1_2_1_x_q1_2_2", "q1_1_1_x_q1_2_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_4_1 = new TProfile("q1_1_1_x_q1_4_1", "q1_1_1_x_q1_4_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_4_2 = new TProfile("q1_1_1_x_q1_4_2", "q1_1_1_x_q1_4_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_4_3 = new TProfile("q1_1_1_x_q1_4_3", "q1_1_1_x_q1_4_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_4_4 = new TProfile("q1_1_1_x_q1_4_4", "q1_1_1_x_q1_4_4", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_1_x_q1_3_1 = new TProfile("q1_2_1_x_q1_3_1", "q1_2_1_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_1_x_q1_3_2 = new TProfile("q1_2_1_x_q1_3_2", "q1_2_1_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_1_x_q1_3_3 = new TProfile("q1_2_1_x_q1_3_3", "q1_2_1_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_2_x_q1_3_1 = new TProfile("q1_2_2_x_q1_3_1", "q1_2_2_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_2_x_q1_3_2 = new TProfile("q1_2_2_x_q1_3_2", "q1_2_2_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_2_x_q1_3_3 = new TProfile("q1_2_2_x_q1_3_3", "q1_2_2_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_5_1 = new TProfile("q1_5_1", "q1_5_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_5_2 = new TProfile("q1_5_2", "q1_5_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_5_3 = new TProfile("q1_5_3", "q1_5_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_5_4 = new TProfile("q1_5_4", "q1_5_4", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_5_5 = new TProfile("q1_5_5", "q1_5_5", kNCentBinsSmall, kCentBinsSmall);

    // 6th order
    TProfile* q6_1_1 = new TProfile("q6_1_1", "q6_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q4_1_1_x_q1_2_1 = new TProfile("q4_1_1_x_q1_2_1", "q4_1_1_x_q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q4_1_1_x_q1_2_2 = new TProfile("q4_1_1_x_q1_2_2", "q4_1_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q3_1_1_x_q1_3_1 = new TProfile("q3_1_1_x_q1_3_1", "q3_1_1_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q3_1_1_x_q1_3_2 = new TProfile("q3_1_1_x_q1_3_2", "q3_1_1_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q3_1_1_x_q1_3_3 = new TProfile("q3_1_1_x_q1_3_3", "q3_1_1_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_2_2_x_q1_2_1 = new TProfile("q2_1_1_x_q1_2_2_x_q1_2_1", "q2_1_1_x_q1_2_2_x_q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q2_2_1 = new TProfile("q2_1_1_x_q2_2_1", "q2_1_1_x_q2_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q2_2_2 = new TProfile("q2_1_1_x_q2_2_2", "q2_1_1_x_q2_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q3_2_1 = new TProfile("q3_2_1", "q3_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q3_2_2 = new TProfile("q3_2_2", "q3_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_4_1 = new TProfile("q2_1_1_x_q1_4_1", "q2_1_1_x_q1_4_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_4_2 = new TProfile("q2_1_1_x_q1_4_2", "q2_1_1_x_q1_4_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_4_3 = new TProfile("q2_1_1_x_q1_4_3", "q2_1_1_x_q1_4_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_1_1_x_q1_4_4 = new TProfile("q2_1_1_x_q1_4_4", "q2_1_1_x_q1_4_4", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_2_1_x_q1_2_2 = new TProfile("q2_2_1_x_q1_2_2", "q2_2_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_2_2_x_q1_2_1 = new TProfile("q2_2_2_x_q1_2_1", "q2_2_2_x_q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_2_1_x_q1_3_1 = new TProfile("q1_1_1_x_q1_2_1_x_q1_3_1", "q1_1_1_x_q1_2_1_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_2_1_x_q1_3_2 = new TProfile("q1_1_1_x_q1_2_1_x_q1_3_2", "q1_1_1_x_q1_2_1_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_2_1_x_q1_3_3 = new TProfile("q1_1_1_x_q1_2_1_x_q1_3_3", "q1_1_1_x_q1_2_1_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_2_2_x_q1_3_1 = new TProfile("q1_1_1_x_q1_2_2_x_q1_3_1", "q1_1_1_x_q1_2_2_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_2_2_x_q1_3_2 = new TProfile("q1_1_1_x_q1_2_2_x_q1_3_2", "q1_1_1_x_q1_2_2_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_2_2_x_q1_3_3 = new TProfile("q1_1_1_x_q1_2_2_x_q1_3_3", "q1_1_1_x_q1_2_2_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_5_1 = new TProfile("q1_1_1_x_q1_5_1", "q1_1_1_x_q1_5_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_5_2 = new TProfile("q1_1_1_x_q1_5_2", "q1_1_1_x_q1_5_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_5_3 = new TProfile("q1_1_1_x_q1_5_3", "q1_1_1_x_q1_5_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_5_4 = new TProfile("q1_1_1_x_q1_5_4", "q1_1_1_x_q1_5_4", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_1_1_x_q1_5_5 = new TProfile("q1_1_1_x_q1_5_5", "q1_1_1_x_q1_5_5", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_1_x_q1_4_1 = new TProfile("q1_2_1_x_q1_4_1", "q1_2_1_x_q1_4_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_1_x_q1_4_2 = new TProfile("q1_2_1_x_q1_4_2", "q1_2_1_x_q1_4_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_1_x_q1_4_3 = new TProfile("q1_2_1_x_q1_4_3", "q1_2_1_x_q1_4_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_1_x_q1_4_4 = new TProfile("q1_2_1_x_q1_4_4", "q1_2_1_x_q1_4_4", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_2_x_q1_4_1 = new TProfile("q1_2_2_x_q1_4_1", "q1_2_2_x_q1_4_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_2_x_q1_4_2 = new TProfile("q1_2_2_x_q1_4_2", "q1_2_2_x_q1_4_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_2_x_q1_4_3 = new TProfile("q1_2_2_x_q1_4_3", "q1_2_2_x_q1_4_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_2_2_x_q1_4_4 = new TProfile("q1_2_2_x_q1_4_4", "q1_2_2_x_q1_4_4", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_3_1 = new TProfile("q2_3_1", "q2_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_3_1_x_q1_3_2 = new TProfile("q1_3_1_x_q1_3_2", "q1_3_1_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_3_1_x_q1_3_3 = new TProfile("q1_3_1_x_q1_3_3", "q1_3_1_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_3_2 = new TProfile("q2_3_2", "q2_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_3_2_x_q1_3_3 = new TProfile("q1_3_2_x_q1_3_3", "q1_3_2_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q2_3_3 = new TProfile("q2_3_3", "q2_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_6_1 = new TProfile("q1_6_1", "q1_6_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_6_2 = new TProfile("q1_6_2", "q1_6_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_6_3 = new TProfile("q1_6_3", "q1_6_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_6_4 = new TProfile("q1_6_4", "q1_6_4", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_6_5 = new TProfile("q1_6_5", "q1_6_5", kNCentBinsSmall, kCentBinsSmall);
    TProfile* q1_6_6 = new TProfile("q1_6_6", "q1_6_6", kNCentBinsSmall, kCentBinsSmall);

    int readError = 0;
    for (int j = 0; j < total_event; j++)
    {
      if (tuple_qmoment->GetEntry(j) < 0) {readError = -999; break;}
      arg = tuple_qmoment->GetArgs();
      centrality = kMultV0M ? arg[0] : arg[13];

      double qPr_p[]{arg[1], arg[3], arg[5], arg[7], arg[9], arg[11]};
      double qPr_n[]{arg[2], arg[4], arg[6], arg[8], arg[10], arg[12]};

      auto qa_b_c = [&](int const a, int const b, int const c) -> double
      {
        double sgn = (b % 2) == 0 ? 1. : -1.;
        return std::pow(/* qPr_p[c - 1] + sgn * */ qPr_n[c - 1], a);
      };

      // full formula
      // 1st order
      q1_1_1->Fill(centrality, qa_b_c(1, 1, 1));

      // 2nd order
      q2_1_1->Fill(centrality, qa_b_c(2, 1, 1));
      q1_2_1->Fill(centrality, qa_b_c(1, 2, 1));
      q1_2_2->Fill(centrality, qa_b_c(1, 2, 2));

      // 3rd order
      q3_1_1->Fill(centrality, qa_b_c(3, 1, 1));
      q1_1_1_x_q1_2_1->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 2, 1));
      q1_1_1_x_q1_2_2->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 2, 2));
      q1_3_1->Fill(centrality, qa_b_c(1, 3, 1));
      q1_3_2->Fill(centrality, qa_b_c(1, 3, 2));
      q1_3_3->Fill(centrality, qa_b_c(1, 3, 3));

      // 4th order
      q4_1_1->Fill(centrality, qa_b_c(4, 1, 1));
      q2_1_1_x_q1_2_1->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 2, 1));
      q2_1_1_x_q1_2_2->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 2, 2));
      q1_1_1_x_q1_3_1->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 3, 1));
      q2_2_1->Fill(centrality, qa_b_c(2, 2, 1));
      q2_2_2->Fill(centrality, qa_b_c(2, 2, 2));
      q1_1_1_x_q1_3_2->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 3, 2));
      q1_1_1_x_q1_3_3->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 3, 3));
      q1_2_1_x_q1_2_2->Fill(centrality, qa_b_c(1, 2, 1) * qa_b_c(1, 2, 2));
      q1_4_1->Fill(centrality, qa_b_c(1, 4, 1));
      q1_4_2->Fill(centrality, qa_b_c(1, 4, 2));
      q1_4_3->Fill(centrality, qa_b_c(1, 4, 3));
      q1_4_4->Fill(centrality, qa_b_c(1, 4, 4));

      // 5th order
      q5_1_1->Fill(centrality, qa_b_c(5, 1, 1));
      q3_1_1_x_q1_2_1->Fill(centrality, qa_b_c(3, 1, 1) * qa_b_c(1, 2, 1));
      q3_1_1_x_q1_2_2->Fill(centrality, qa_b_c(3, 1, 1) * qa_b_c(1, 2, 2));
      q2_1_1_x_q1_3_1->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 3, 1));
      q2_1_1_x_q1_3_2->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 3, 2));
      q2_1_1_x_q1_3_3->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 3, 3));
      q2_2_2_x_q1_1_1->Fill(centrality, qa_b_c(2, 2, 2) * qa_b_c(1, 1, 1));
      q2_2_1_x_q1_1_1->Fill(centrality, qa_b_c(2, 2, 1) * qa_b_c(1, 1, 1));
      q1_1_1_x_q1_2_1_x_q1_2_2->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 2, 1) * qa_b_c(1, 2, 2));
      q1_1_1_x_q1_4_1->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 4, 1));
      q1_1_1_x_q1_4_2->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 4, 2));
      q1_1_1_x_q1_4_3->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 4, 3));
      q1_1_1_x_q1_4_4->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 4, 4));
      q1_2_1_x_q1_3_1->Fill(centrality, qa_b_c(1, 2, 1) * qa_b_c(1, 3, 1));
      q1_2_1_x_q1_3_2->Fill(centrality, qa_b_c(1, 2, 1) * qa_b_c(1, 3, 2));
      q1_2_1_x_q1_3_3->Fill(centrality, qa_b_c(1, 2, 1) * qa_b_c(1, 3, 3));
      q1_2_2_x_q1_3_1->Fill(centrality, qa_b_c(1, 2, 2) * qa_b_c(1, 3, 1));
      q1_2_2_x_q1_3_2->Fill(centrality, qa_b_c(1, 2, 2) * qa_b_c(1, 3, 2));
      q1_2_2_x_q1_3_3->Fill(centrality, qa_b_c(1, 2, 2) * qa_b_c(1, 3, 3));
      q1_5_1->Fill(centrality, qa_b_c(1, 5, 1));
      q1_5_2->Fill(centrality, qa_b_c(1, 5, 2));
      q1_5_3->Fill(centrality, qa_b_c(1, 5, 3));
      q1_5_4->Fill(centrality, qa_b_c(1, 5, 4));
      q1_5_5->Fill(centrality, qa_b_c(1, 5, 5));

      // 6th order
      q6_1_1->Fill(centrality, qa_b_c(6, 1, 1));
      q4_1_1_x_q1_2_1->Fill(centrality, qa_b_c(4, 1, 1) * qa_b_c(1, 2, 1));
      q4_1_1_x_q1_2_2->Fill(centrality, qa_b_c(4, 1, 1) * qa_b_c(1, 2, 2));
      q3_1_1_x_q1_3_1->Fill(centrality, qa_b_c(3, 1, 1) * qa_b_c(1, 3, 1));
      q3_1_1_x_q1_3_2->Fill(centrality, qa_b_c(3, 1, 1) * qa_b_c(1, 3, 2));
      q3_1_1_x_q1_3_3->Fill(centrality, qa_b_c(3, 1, 1) * qa_b_c(1, 3, 3));
      q2_1_1_x_q1_2_2_x_q1_2_1->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 2, 2) * qa_b_c(1, 2, 1));
      q2_1_1_x_q2_2_1->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(2, 2, 1));
      q2_1_1_x_q2_2_2->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(2, 2, 2));
      q3_2_1->Fill(centrality, qa_b_c(3, 2, 1));
      q3_2_2->Fill(centrality, qa_b_c(3, 2, 2));
      q2_1_1_x_q1_4_1->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 4, 1));
      q2_1_1_x_q1_4_2->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 4, 2));
      q2_1_1_x_q1_4_3->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 4, 3));
      q2_1_1_x_q1_4_4->Fill(centrality, qa_b_c(2, 1, 1) * qa_b_c(1, 4, 4));
      q2_2_1_x_q1_2_2->Fill(centrality, qa_b_c(2, 2, 1) * qa_b_c(1, 2, 2));
      q2_2_2_x_q1_2_1->Fill(centrality, qa_b_c(2, 2, 2) * qa_b_c(1, 2, 1));
      q1_1_1_x_q1_2_1_x_q1_3_1->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 2, 1) * qa_b_c(1, 3, 1));
      q1_1_1_x_q1_2_1_x_q1_3_2->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 2, 1) * qa_b_c(1, 3, 2));
      q1_1_1_x_q1_2_1_x_q1_3_3->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 2, 1) * qa_b_c(1, 3, 3));
      q1_1_1_x_q1_2_2_x_q1_3_1->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 2, 2) * qa_b_c(1, 3, 1));
      q1_1_1_x_q1_2_2_x_q1_3_2->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 2, 2) * qa_b_c(1, 3, 2));
      q1_1_1_x_q1_2_2_x_q1_3_3->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 2, 2) * qa_b_c(1, 3, 3));
      q1_1_1_x_q1_5_1->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 5, 1));
      q1_1_1_x_q1_5_2->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 5, 2));
      q1_1_1_x_q1_5_3->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 5, 3));
      q1_1_1_x_q1_5_4->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 5, 4));
      q1_1_1_x_q1_5_5->Fill(centrality, qa_b_c(1, 1, 1) * qa_b_c(1, 5, 5));
      q1_2_1_x_q1_4_1->Fill(centrality, qa_b_c(1, 2, 1) * qa_b_c(1, 4, 1));
      q1_2_1_x_q1_4_2->Fill(centrality, qa_b_c(1, 2, 1) * qa_b_c(1, 4, 2));
      q1_2_1_x_q1_4_3->Fill(centrality, qa_b_c(1, 2, 1) * qa_b_c(1, 4, 3));
      q1_2_1_x_q1_4_4->Fill(centrality, qa_b_c(1, 2, 1) * qa_b_c(1, 4, 4));
      q1_2_2_x_q1_4_1->Fill(centrality, qa_b_c(1, 2, 2) * qa_b_c(1, 4, 1));
      q1_2_2_x_q1_4_2->Fill(centrality, qa_b_c(1, 2, 2) * qa_b_c(1, 4, 2));
      q1_2_2_x_q1_4_3->Fill(centrality, qa_b_c(1, 2, 2) * qa_b_c(1, 4, 3));
      q1_2_2_x_q1_4_4->Fill(centrality, qa_b_c(1, 2, 2) * qa_b_c(1, 4, 4));
      q2_3_1->Fill(centrality, qa_b_c(2, 3, 1));
      q1_3_1_x_q1_3_2->Fill(centrality, qa_b_c(1, 3, 1) * qa_b_c(1, 3, 2));
      q1_3_1_x_q1_3_3->Fill(centrality, qa_b_c(1, 3, 1) * qa_b_c(1, 3, 3));
      q2_3_2->Fill(centrality, qa_b_c(2, 3, 2));
      q1_3_2_x_q1_3_3->Fill(centrality, qa_b_c(1, 3, 2) * qa_b_c(1, 3, 3));
      q2_3_3->Fill(centrality, qa_b_c(2, 3, 3));
      q1_6_1->Fill(centrality, qa_b_c(1, 6, 1));
      q1_6_2->Fill(centrality, qa_b_c(1, 6, 2));
      q1_6_3->Fill(centrality, qa_b_c(1, 6, 3));
      q1_6_4->Fill(centrality, qa_b_c(1, 6, 4));
      q1_6_5->Fill(centrality, qa_b_c(1, 6, 5));
      q1_6_6->Fill(centrality, qa_b_c(1, 6, 6));
    }
    if (readError < 0) {
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }

    if (FILL_MC) {
      for (int j_gen = 0; j_gen < total_event_gen; j_gen++)
      {
        if (tuple_moment_gen->GetEntry(j_gen) < 0) {readError = -999; break;}
        arg_gen = tuple_moment_gen->GetArgs();
        centrality = kMultV0M ? arg_gen[0] : arg_gen[3];

        double qP1_p = arg_gen[1];
        double qP1_n = arg_gen[2];

        // full formula (gen)
        N1p->Fill(centrality, /* qP1_p + */ qP1_n);
        N1->Fill(centrality, /* qP1_p - */ qP1_n);
        N2->Fill(centrality, std::pow(/* qP1_p - */ qP1_n, 2));
        N3->Fill(centrality, std::pow(/* qP1_p - */ qP1_n, 3));
        N4->Fill(centrality, std::pow(/* qP1_p - */ qP1_n, 4));
        N5->Fill(centrality, std::pow(/* qP1_p - */ qP1_n, 5));
        N6->Fill(centrality, std::pow(/* qP1_p - */ qP1_n, 6));
      }
      if (readError < 0) {
        std::cout << "no input, skip" << std::endl;
        skippedVar += 1;

        delete fin;
        continue;
      }
    }

    fout.mkdir(Form("var_%d", iVar));
    fout.cd(Form("var_%d", iVar));

    // full formula
    // 1st order
    q1_1_1->Write();

    // 2nd order
    q2_1_1->Write();
    q1_2_1->Write();
    q1_2_2->Write();

    // 3rd order
    q3_1_1->Write();
    q1_1_1_x_q1_2_1->Write();
    q1_1_1_x_q1_2_2->Write();
    q1_3_1->Write();
    q1_3_2->Write();
    q1_3_3->Write();

    // 4th order
    q4_1_1->Write();
    q2_1_1_x_q1_2_1->Write();
    q2_1_1_x_q1_2_2->Write();
    q1_1_1_x_q1_3_1->Write();
    q2_2_1->Write();
    q2_2_2->Write();
    q1_1_1_x_q1_3_2->Write();
    q1_1_1_x_q1_3_3->Write();
    q1_2_1_x_q1_2_2->Write();
    q1_4_1->Write();
    q1_4_2->Write();
    q1_4_3->Write();
    q1_4_4->Write();

    // 5th order
    q5_1_1->Write();
    q3_1_1_x_q1_2_1->Write();
    q3_1_1_x_q1_2_2->Write();
    q2_1_1_x_q1_3_1->Write();
    q2_1_1_x_q1_3_2->Write();
    q2_1_1_x_q1_3_3->Write();
    q2_2_2_x_q1_1_1->Write();
    q2_2_1_x_q1_1_1->Write();
    q1_1_1_x_q1_2_1_x_q1_2_2->Write();
    q1_1_1_x_q1_4_1->Write();
    q1_1_1_x_q1_4_2->Write();
    q1_1_1_x_q1_4_3->Write();
    q1_1_1_x_q1_4_4->Write();
    q1_2_1_x_q1_3_1->Write();
    q1_2_1_x_q1_3_2->Write();
    q1_2_1_x_q1_3_3->Write();
    q1_2_2_x_q1_3_1->Write();
    q1_2_2_x_q1_3_2->Write();
    q1_2_2_x_q1_3_3->Write();
    q1_5_1->Write();
    q1_5_2->Write();
    q1_5_3->Write();
    q1_5_4->Write();
    q1_5_5->Write();

    // 6th order
    q6_1_1->Write();
    q4_1_1_x_q1_2_1->Write();
    q4_1_1_x_q1_2_2->Write();
    q3_1_1_x_q1_3_1->Write();
    q3_1_1_x_q1_3_2->Write();
    q3_1_1_x_q1_3_3->Write();
    q2_1_1_x_q1_2_2_x_q1_2_1->Write();
    q2_1_1_x_q2_2_1->Write();
    q2_1_1_x_q2_2_2->Write();
    q3_2_1->Write();
    q3_2_2->Write();
    q2_1_1_x_q1_4_1->Write();
    q2_1_1_x_q1_4_2->Write();
    q2_1_1_x_q1_4_3->Write();
    q2_1_1_x_q1_4_4->Write();
    q2_2_1_x_q1_2_2->Write();
    q2_2_2_x_q1_2_1->Write();
    q1_1_1_x_q1_2_1_x_q1_3_1->Write();
    q1_1_1_x_q1_2_1_x_q1_3_2->Write();
    q1_1_1_x_q1_2_1_x_q1_3_3->Write();
    q1_1_1_x_q1_2_2_x_q1_3_1->Write();
    q1_1_1_x_q1_2_2_x_q1_3_2->Write();
    q1_1_1_x_q1_2_2_x_q1_3_3->Write();
    q1_1_1_x_q1_5_1->Write();
    q1_1_1_x_q1_5_2->Write();
    q1_1_1_x_q1_5_3->Write();
    q1_1_1_x_q1_5_4->Write();
    q1_1_1_x_q1_5_5->Write();
    q1_2_1_x_q1_4_1->Write();
    q1_2_1_x_q1_4_2->Write();
    q1_2_1_x_q1_4_3->Write();
    q1_2_1_x_q1_4_4->Write();
    q1_2_2_x_q1_4_1->Write();
    q1_2_2_x_q1_4_2->Write();
    q1_2_2_x_q1_4_3->Write();
    q1_2_2_x_q1_4_4->Write();
    q2_3_1->Write();
    q1_3_1_x_q1_3_2->Write();
    q1_3_1_x_q1_3_3->Write();
    q2_3_2->Write();
    q1_3_2_x_q1_3_3->Write();
    q2_3_3->Write();
    q1_6_1->Write();
    q1_6_2->Write();
    q1_6_3->Write();
    q1_6_4->Write();
    q1_6_5->Write();
    q1_6_6->Write();

    if (FILL_MC) {
      N1p->Write();
      N1->Write();
      N2->Write();
      N3->Write();
      N4->Write();
      N5->Write();
      N6->Write();
    }

    fout.Close();
    fin->Close();
    delete fin;
  }

  std::cout << "skipped = " << skippedVar << std::endl;
  w.Stop();
  w.Print();

}
