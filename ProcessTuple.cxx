#include "utils.h"

#include <TStopwatch.h>
#include <TFile.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <Riostream.h>

#define FILL_MC

void ProcessTuple(int smpl = 9, int iVarMin = 0, int iVarMax = 3)
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

    TNtuple *tuple_qmoment = (TNtuple*)fin->Get(Form("evtTuple_%d", iVar));
    if (!tuple_qmoment)
    {
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }
    tuple_qmoment->SetCacheSize(0);

    #ifdef FILL_MC
      TNtuple *tuple_moment_gen = (TNtuple*)fin->Get(Form("evtTupleGen_%d", iVar));
      tuple_moment_gen->SetCacheSize(0);
    #endif // FILL_MC

    int evt[10] = {0};
    int centrality;

    float *arg;
    float total_event = tuple_qmoment->GetEntriesFast();

    #ifdef FILL_MC
      float *arg_gen;
      float total_event_gen = tuple_moment_gen->GetEntriesFast();
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
    TProfile *q1_1_1 = new TProfile("q1_1_1", "q1_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q2_1_1 = new TProfile("q2_1_1", "q2_1_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_1 = new TProfile("q1_2_1", "q1_2_1", kNCentBinsSmall, kCentBinsSmall);
    TProfile *q1_2_2 = new TProfile("q1_2_2", "q1_2_2", kNCentBinsSmall, kCentBinsSmall);

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

      // full formula
      q1_1_1->Fill(centrality, qP1_p - qP1_n);
      q2_1_1->Fill(centrality, powI(qP1_p - qP1_n, 2.));
      q1_2_1->Fill(centrality, qP1_p + qP1_n);
      q1_2_2->Fill(centrality, qP2_p + qP2_n);
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

      // full formula
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
    q1_1_1->Write();
    q2_1_1->Write();
    q1_2_1->Write();
    q1_2_2->Write();

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
