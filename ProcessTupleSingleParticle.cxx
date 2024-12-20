#include "utils.h"

#include <TStopwatch.h>
#include <TFile.h>
#include <TProfile.h>
#include <TNtupleD.h>
#include <Riostream.h>

void ProcessTupleSingleParticle(int smpl = 0, int iVarMin = 364, int iVarMax = 365, bool dummy = false)
{
  TStopwatch w;
  w.Start();

  int sample = smpl;
  int skippedVar = 0;

  for (int iVar = iVarMin; iVar < iVarMax; ++iVar){

    TFile *fin = TFile::Open(Form("%s/LHC18pp%d_var_%d.root", kResDir, sample, iVar));
    if (!fin || fin->TestBit(TFile::kZombie)){
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }
    TFile fout(Form("%s/output_sys_singleParticle_%d_%d.root", kResDir, sample, iVar), "recreate");

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

    int evt[10] = {0};
    int centrality;

    double *arg;
    double total_event = tuple_qmoment->GetEntriesFast();

    //std::cout << total_event << "\t" << total_event2 << std::endl;

    TProfile *q1_pr_p = new TProfile("q1_pr_p", "q1_pr_p", kNCentBins, kCentBins);
    TProfile *q2_pr_p = new TProfile("q2_pr_p", "q2_pr_p", kNCentBins, kCentBins);
    TProfile *q1_pr_n = new TProfile("q1_pr_n", "q1_pr_n", kNCentBins, kCentBins);
    TProfile *q2_pr_n = new TProfile("q2_pr_n", "q2_pr_n", kNCentBins, kCentBins);
    TProfile *q1square_pr_p = new TProfile("q1square_pr_p", "q1square_pr_p", kNCentBins, kCentBins);
    TProfile *q1square_pr_n = new TProfile("q1square_pr_n", "q1square_pr_n", kNCentBins, kCentBins);
    TProfile *q1_pr_pr_pn = new TProfile("q1_pr_pr_pn", "q1_pr_pr_pn", kNCentBins, kCentBins);

    int readError = 0;
    for (int j = 0; j < total_event; j++)
    {
      if (tuple_qmoment->GetEntry(j) < 0) {readError = -999; break;}
      arg = tuple_qmoment->GetArgs();
      centrality = arg[0];
      double q1_p = arg[1];
      double q2_p = arg[3];
      double q1_n = arg[2];
      double q2_n = arg[4];

      q1_pr_p->Fill(centrality, q1_p);
      q2_pr_p->Fill(centrality, q2_p);
      q1_pr_n->Fill(centrality, q1_n);
      q2_pr_n->Fill(centrality, q2_n);
      q1square_pr_p->Fill(centrality, q1_p * q1_p);
      q1square_pr_n->Fill(centrality, q1_n * q1_n);
      q1_pr_pr_pn->Fill(centrality, q1_p * q1_n);
    }
    if (readError < 0) {
      std::cout << "no input, skip" << std::endl;
      skippedVar += 1;

      delete fin;
      continue;
    }

    fout.mkdir(Form("var_%d", iVar));
    fout.cd(Form("var_%d", iVar));

    q1_pr_p->Write();
    q2_pr_p->Write();
    q1_pr_n->Write();
    q2_pr_n->Write();
    q1square_pr_p->Write();
    q1square_pr_n->Write();
    q1_pr_pr_pn->Write();

    fout.Close();
    fin->Close();
    delete fin;
  }

  std::cout << "skipped = " << skippedVar << std::endl;
  w.Stop();
  w.Print();

}