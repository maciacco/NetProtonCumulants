import ROOT

dict = {
'v2' : ['out_sys_18_binning_fine_3_mix_finalBinning_singleParticleHighOrder_k2k1.root', '#it{#kappa}_{2}(#it{V})/#LT#it{V}#GT^{2}'],
'v3' : ['out_sys_18_binning_fine_3_mix_finalBinning_singleParticleHighOrder_k3k1.root', '#it{#kappa}_{3}(#it{V})/#LT#it{V}#GT^{3}']
}

file = ROOT.TFile(dict['v3'][0])
graph = file.Get('g_81')
#print(f'x = {graph.GetPointY(0)}')

graph.SetMarkerStyle(20)
graph.SetMarkerSize(1.2)
graph.GetXaxis().SetTitle('#LTd#it{N}_{ch}/d#it{#eta}#GT')
graph.GetYaxis().SetTitle(dict['v3'][1])
graph.GetYaxis().SetTitleOffset(1.2)

tx = ROOT.TLatex()
tx.SetTextFont(45)
tx.SetTextSize(27)
tx.SetNDC()

canv = ROOT.TCanvas('c', 'c', 600, 600)
canv.SetLeftMargin(0.16)
canv.SetTopMargin(0.03)
canv.SetRightMargin(0.03)
canv.SetBottomMargin(0.16)
canv.cd()
graph.Draw('ape')
tx.DrawLatex(0.35, 0.87, 'ALICE Work in progress')
tx.DrawLatex(0.35, 0.8, 'pp #sqrt{#it{s}} = 13 TeV')
tx.DrawLatex(0.35, 0.73, '|#it{#eta}| < 0.8, 0.5 < #it{p}_{T} < 1.5 GeV/#it{c}')
canv.Print('evMix_v2.pdf')
