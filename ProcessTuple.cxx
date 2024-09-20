#include "utils.h"

#include <TStopwatch.h>
#include <TFile.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <Riostream.h>

#define FILL_MC

double qa_b(double const p, double const n, int const a, int const b){
  double sgn = (b % 2) == 0 ? 1. : -1.;
  return powI(p + sgn * n, a);
}

void ProcessTuple(int smpl = 0, int iVarMin = 0, int iVarMax = 3)
{
  TStopwatch w;
  w.Start();

  int sample = smpl;
  int skippedVar = 0;

  for (int iVar = iVarMin; iVar < iVarMax; ++iVar){

    TFile *fin = TFile::Open(Form("%s/LHC18%d_var_%d.root", kResDir, sample, iVar));
    if (!fin || fin->TestBit(TFile::kZombie)){
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }
    TFile fout(Form("%s/output_sys_%d_%d.root", kResDir, sample, iVar), "recreate");

    TNtupleD *tuple_qmoment = (TNtupleD*)fin->Get(Form("evtTuple_%d", iVar));
    if (!tuple_qmoment)
    {
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }
    tuple_qmoment->SetCacheSize(0);

    #ifdef FILL_MC
      TNtupleD *tuple_moment_gen = (TNtupleD*)fin->Get(Form("evtTupleGen_%d", iVar));
      tuple_moment_gen->SetCacheSize(0);
    #endif // FILL_MC

    int evt[10] = {0};
    int centrality;

    double *arg;
    double total_event = tuple_qmoment->GetEntriesFast();

    #ifdef FILL_MC
      double *arg_gen;
      double total_event_gen = tuple_moment_gen->GetEntriesFast();
      // non-central moments of generated distribution
      TProfile *N1p = new TProfile("N1p", "N1p", kNCentBinsSmall, kCentBinsSmall);
      TProfile *N1 = new TProfile("N1", "N1", kNCentBinsSmall, kCentBinsSmall);
      TProfile *N2 = new TProfile("N2", "N2", kNCentBinsSmall, kCentBinsSmall);
      TProfile *N3 = new TProfile("N3", "N3", kNCentBinsSmall, kCentBinsSmall);
      TProfile *N4 = new TProfile("N4", "N4", kNCentBinsSmall, kCentBinsSmall);
      TProfile *N5 = new TProfile("N5", "N5", kNCentBinsSmall, kCentBinsSmall);
      TProfile *N6 = new TProfile("N6", "N6", kNCentBinsSmall, kCentBinsSmall);
    #endif // FILL_MC

    // full formula
    // 1st order
    TProfile *q1_1_1 = new TProfile("q1_1_1", "q1_1_1", kNCentBinsSmall, kCentBinsSmall);

    // 2nd order
    TProfile *q2_1_1 = new TProfile("q2_1_1", "q2_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_1 = new TProfile("q1_2_1", "q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_2 = new TProfile("q1_2_2", "q1_2_2", kNCentBinsSmall, kCentBinsSmall);

    // 3rd order
    TProfile *q3_1_1 = new TProfile("q3_1_1", "q3_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_2_1 = new TProfile("q1_1_1_x_q1_2_1", "q1_1_1_x_q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_2_2 = new TProfile("q1_1_1_x_q1_2_2", "q1_1_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_3_1 = new TProfile("q1_3_1", "q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_3_2 = new TProfile("q1_3_2", "q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_3_3 = new TProfile("q1_3_3", "q1_3_3", kNCentBinsSmall, kCentBinsSmall);

    // // 4th order
    TProfile *q4_1_1 = new TProfile("q4_1_1", "q4_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_1_1_x_q1_2_1 = new TProfile("q2_1_1_x_q1_2_1", "q2_1_1_x_q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_1_1_x_q1_2_2 = new TProfile("q2_1_1_x_q1_2_2", "q2_1_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_3_1 = new TProfile("q1_1_1_x_q1_3_1", "q1_1_1_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_2_1 = new TProfile("q2_2_1", "q2_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_2_2 = new TProfile("q2_2_2", "q2_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_3_2 = new TProfile("q1_1_1_x_q1_3_2", "q1_1_1_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_3_3 = new TProfile("q1_1_1_x_q1_3_3", "q1_1_1_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_1_x_q1_2_2 = new TProfile("q1_2_1_x_q1_2_2", "q1_2_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_4_1 = new TProfile("q1_4_1", "q1_4_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_4_2 = new TProfile("q1_4_2", "q1_4_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_4_3 = new TProfile("q1_4_3", "q1_4_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_4_4 = new TProfile("q1_4_4", "q1_4_4", kNCentBinsSmall, kCentBinsSmall);

    // 5th order
    TProfile *q5_1_1 = new TProfile("q5_1_1", "q5_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q3_1_1_x_q1_2_1 = new TProfile("q3_1_1_x_q1_2_1", "q3_1_1_x_q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q3_1_1_x_q1_2_2 = new TProfile("q3_1_1_x_q1_2_2", "q3_1_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_1_1_x_q1_3_1 = new TProfile("q2_1_1_x_q1_3_1", "q2_1_1_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_1_1_x_q1_3_2 = new TProfile("q2_1_1_x_q1_3_2", "q2_1_1_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_1_1_x_q1_3_3 = new TProfile("q2_1_1_x_q1_3_3", "q2_1_1_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_2_2_x_q1_1_1 = new TProfile("q2_2_2_x_q1_1_1", "q2_2_2_x_q1_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_2_1_x_q1_1_1 = new TProfile("q2_2_1_x_q1_1_1", "q2_2_1_x_q1_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_2_1_x_q1_2_2 = new TProfile("q1_1_1_x_q1_2_1_x_q1_2_2", "q1_1_1_x_q1_2_1_x_q1_2_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_4_1 = new TProfile("q1_1_1_x_q1_4_1", "q1_1_1_x_q1_4_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_4_2 = new TProfile("q1_1_1_x_q1_4_2", "q1_1_1_x_q1_4_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_4_3 = new TProfile("q1_1_1_x_q1_4_3", "q1_1_1_x_q1_4_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_1_1_x_q1_4_4 = new TProfile("q1_1_1_x_q1_4_4", "q1_1_1_x_q1_4_4", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_1_x_q1_3_1 = new TProfile("q1_2_1_x_q1_3_1", "q1_2_1_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_1_x_q1_3_2 = new TProfile("q1_2_1_x_q1_3_2", "q1_2_1_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_1_x_q1_3_3 = new TProfile("q1_2_1_x_q1_3_3", "q1_2_1_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_2_x_q1_3_1 = new TProfile("q1_2_2_x_q1_3_1", "q1_2_2_x_q1_3_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_2_x_q1_3_2 = new TProfile("q1_2_2_x_q1_3_2", "q1_2_2_x_q1_3_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_2_x_q1_3_3 = new TProfile("q1_2_2_x_q1_3_3", "q1_2_2_x_q1_3_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_5_1 = new TProfile("q1_5_1", "q1_5_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_5_2 = new TProfile("q1_5_2", "q1_5_2", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_5_3 = new TProfile("q1_5_3", "q1_5_3", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_5_4 = new TProfile("q1_5_4", "q1_5_4", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_5_5 = new TProfile("q1_5_5", "q1_5_5", kNCentBinsSmall, kCentBinsSmall);

    // 6th order
    TProfile *q6_1_1 = new TProfile("q6_1_1", "q6_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q4_1_1_x_q1_2_1 = new TProfile("q4_1_1_x_q1_2_", "q4_1_1_x_q1_2_", kNCentBinsSmall, kCentBinsSmall);


    int readError = 0;
    for (int j = 0; j < total_event; j++)
    {
      if (tuple_qmoment->GetEntry(j) < 0) {readError = -999; break;}
      arg = tuple_qmoment->GetArgs();
      centrality = arg[0];

      double qP1_p = arg[1];
      double qP1_n = arg[2];
      double qP2_p = arg[3];
      double qP2_n = arg[4];
      double qP3_p = arg[5];
      double qP3_n = arg[6];
      double qP4_p = arg[7];
      double qP4_n = arg[8];
      double qP5_p = arg[9];
      double qP5_n = arg[10];
      double qP6_p = arg[11];
      double qP6_n = arg[12];

      // full formula
      // 1st order
      q1_1_1->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1));

      // 2nd order
      q2_1_1->Fill(centrality, qa_b(qP1_p, qP1_n, 2, 1));
      q1_2_1->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 2));
      q1_2_2->Fill(centrality, qa_b(qP2_p, qP2_n, 1, 2));

      // 3rd order
      q3_1_1->Fill(centrality, qa_b(qP1_p, qP1_n, 3, 1));
      q1_1_1_x_q1_2_1->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP1_p, qP1_n, 1, 2));
      q1_1_1_x_q1_2_2->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP2_p, qP2_n, 1, 2));
      q1_3_1->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 3));
      q1_3_2->Fill(centrality, qa_b(qP2_p, qP2_n, 1, 3));
      q1_3_3->Fill(centrality, qa_b(qP3_p, qP3_n, 1, 3));

      // 4th order
      q4_1_1->Fill(centrality, qa_b(qP1_p, qP1_n, 4, 1));
      q2_1_1_x_q1_2_1->Fill(centrality, qa_b(qP1_p, qP1_n, 2, 1) * qa_b(qP1_p, qP1_n, 1, 2));
      q2_1_1_x_q1_2_2->Fill(centrality, qa_b(qP1_p, qP1_n, 2, 1) * qa_b(qP2_p, qP2_n, 1, 2));
      q1_1_1_x_q1_3_1->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP1_p, qP1_n, 1, 3));
      q2_2_1->Fill(centrality, qa_b(qP1_p, qP1_n, 2, 2));
      q2_2_2->Fill(centrality, qa_b(qP2_p, qP2_n, 2, 2));
      q1_1_1_x_q1_3_2->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP2_p, qP2_n, 1, 3));
      q1_1_1_x_q1_3_3->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP3_p, qP3_n, 1, 3));
      q1_2_1_x_q1_2_2->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 2) * qa_b(qP2_p, qP2_n, 1, 2));
      q1_4_1->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 4));
      q1_4_2->Fill(centrality, qa_b(qP2_p, qP2_n, 1, 4));
      q1_4_3->Fill(centrality, qa_b(qP3_p, qP3_n, 1, 4));
      q1_4_4->Fill(centrality, qa_b(qP4_p, qP4_n, 1, 4));

      // 5th order
      q5_1_1->Fill(centrality, qa_b(qP1_p, qP1_n, 5, 1));
      q3_1_1_x_q1_2_1->Fill(centrality, qa_b(qP1_p, qP1_n, 3, 1) * qa_b(qP1_p, qP1_n, 1, 2));
      q3_1_1_x_q1_2_2->Fill(centrality, qa_b(qP1_p, qP1_n, 3, 1) * qa_b(qP2_p, qP2_n, 1, 2));
      q2_1_1_x_q1_3_1->Fill(centrality, qa_b(qP1_p, qP1_n, 2, 1) * qa_b(qP1_p, qP1_n, 1, 3));
      q2_1_1_x_q1_3_2->Fill(centrality, qa_b(qP1_p, qP1_n, 2, 1) * qa_b(qP2_p, qP2_n, 1, 3));
      q2_1_1_x_q1_3_3->Fill(centrality, qa_b(qP1_p, qP1_n, 2, 1) * qa_b(qP3_p, qP3_n, 1, 3));
      q2_2_2_x_q1_1_1->Fill(centrality, qa_b(qP2_p, qP2_n, 2, 2) * qa_b(qP1_p, qP1_n, 1, 1));
      q2_2_1_x_q1_1_1->Fill(centrality, qa_b(qP1_p, qP1_n, 2, 2) * qa_b(qP1_p, qP1_n, 1, 1));
      q1_1_1_x_q1_2_1_x_q1_2_2->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP1_p, qP1_n, 1, 2) * qa_b(qP2_p, qP2_n, 1, 2));
      q1_1_1_x_q1_4_1->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP1_p, qP1_n, 1, 4));
      q1_1_1_x_q1_4_2->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP2_p, qP2_n, 1, 4));
      q1_1_1_x_q1_4_3->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP3_p, qP3_n, 1, 4));
      q1_1_1_x_q1_4_4->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 1) * qa_b(qP4_p, qP4_n, 1, 4));
      q1_2_1_x_q1_3_1->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 2) * qa_b(qP1_p, qP1_n, 1, 3));
      q1_2_1_x_q1_3_2->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 2) * qa_b(qP2_p, qP2_n, 1, 3));
      q1_2_1_x_q1_3_3->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 2) * qa_b(qP3_p, qP3_n, 1, 3));
      q1_2_2_x_q1_3_1->Fill(centrality, qa_b(qP2_p, qP2_n, 1, 2) * qa_b(qP1_p, qP1_n, 1, 3));
      q1_2_2_x_q1_3_2->Fill(centrality, qa_b(qP2_p, qP2_n, 1, 2) * qa_b(qP2_p, qP2_n, 1, 3));
      q1_2_2_x_q1_3_3->Fill(centrality, qa_b(qP2_p, qP2_n, 1, 2) * qa_b(qP3_p, qP3_n, 1, 3));
      q1_5_1->Fill(centrality, qa_b(qP1_p, qP1_n, 1, 5));
      q1_5_2->Fill(centrality, qa_b(qP2_p, qP2_n, 1, 5));
      q1_5_3->Fill(centrality, qa_b(qP3_p, qP3_n, 1, 5));
      q1_5_4->Fill(centrality, qa_b(qP4_p, qP4_n, 1, 5));
      q1_5_5->Fill(centrality, qa_b(qP5_p, qP5_n, 1, 5));
    }
    if (readError < 0) {
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }

    #ifdef FILL_MC
    for (int j_gen = 0; j_gen < total_event_gen; j_gen++)
    {
      if (tuple_moment_gen->GetEntry(j_gen) < 0) {readError = -999; break;}
      arg_gen = tuple_moment_gen->GetArgs();
      centrality = arg_gen[0];

      double qP1_p = arg_gen[1];
      double qP1_n = arg_gen[2];

      // full formula (gen)
      N1p->Fill(centrality, qP1_p + qP1_n);
      N1->Fill(centrality, qP1_p - qP1_n);
      N2->Fill(centrality, powI(qP1_p - qP1_n, 2));
      N3->Fill(centrality, powI(qP1_p - qP1_n, 3));
      N4->Fill(centrality, powI(qP1_p - qP1_n, 4));
      N5->Fill(centrality, powI(qP1_p - qP1_n, 5));
      N6->Fill(centrality, powI(qP1_p - qP1_n, 6));
    }
    if (readError < 0) {
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }
    #endif // FILL_MC

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

    #ifdef FILL_MC
      N1p->Write();
      N1->Write();
      N2->Write();
      N3->Write();
      N4->Write();
      N5->Write();
      N6->Write();
    #endif // FILL_MC

    fout.Close();
    fin->Close();
    delete fin;
  }

  std::cout << "skipped = " << skippedVar << std::endl;
  w.Stop();
  w.Print();

}
