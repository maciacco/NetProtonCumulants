void InterpMult(){
  double pct[]{0.5, 3.5, 7.5, 12.5, 17.5, 25., 35., 45., 60., 85.};
  double mult[]{26.01, 19.99, 16.18, 13.78, 12.01, 10.03, 7.95, 6.32, 4.49, 2.54};
  double pct_err[]{0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double mult_err[]{0.24, 0.2, 0.2, 0.15, 0.15, 0.12, 0.13, 0.11, 0.08, 0.06, 0.04};
  TGraphErrors gg(10, pct, mult, pct_err, mult_err);

  // first interp
  TF1 f1("f1", "[0]*TMath::Power([1] + [2] * x, [3])", 0., 20.);
  f1.SetParLimits(0, -10., 20.);
  f1.SetParLimits(1, -10., 10.);
  f1.SetParLimits(2, -10., 10.);
  f1.SetParLimits(3, -10., 10.);
  f1.SetParameter(0, 9.97);
  f1.SetParameter(1, 0.1);
  f1.SetParameter(2, 0.03);
  f1.SetParameter(3, -0.5);
  gg.Fit("f1", "RM+", "", 0., 20.);

  for (int i{1}; i < 4; ++i) {
    std::cout << i * 5 + 2.5 << "\t" << f1.Eval(i * 5 + 2.5) << std::endl;
  }

  // second interp
  TF1 f2("f2", "[0]*TMath::Power([1] + [2] * x, [3])", 0., 20.);
  f2.SetParLimits(0, -10., 20.);
  f2.SetParLimits(1, -10., 10.);
  f2.SetParLimits(2, -10., 10.);
  f2.SetParLimits(3, -10., 10.);
  f2.SetParameter(0, 9.97);
  f2.SetParameter(1, 0.1);
  f2.SetParameter(2, 0.03);
  f2.SetParameter(3, -0.5);
  gg.Fit("f2", "RM+", "", 14., 100.);

  for (int i{4}; i < 20; ++i) {
    std::cout << i * 5 + 2.5 << "\t" << f2.Eval(i * 5 + 2.5) << std::endl;
  }


  TFile *fo = TFile::Open("f_int.root", "recreate");
  fo->cd();
  gg.Write();
  fo->Close();
}
