import ROOT

file = []
graph = []
color = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2]
names = ['#it{N}_{ch} #rightarrow #it{N}_{trkl}^{|#it{#eta}| < 0.8}', '#it{N}_{ch} #rightarrow #it{N}_{trkl}^{|#it{#eta}| < 1.2}', '#it{N}_{ch} #rightarrow #it{N}_{trk}^{|#it{#eta}| < 0.8}']

ROOT.gStyle.SetOptStat(0)

file.append(ROOT.TFile('out_sys_18_binning_fine_3_mix_finalBinning_singleParticleHighOrder_k2k1.root'))
file.append(ROOT.TFile('out_sys_18_binning_fine_all_mix_finalBinning_singleParticleHighOrder_k2k1.root'))
file.append(ROOT.TFile('out_sys_18_binning_fine_tracks_mix_finalBinning_singleParticleHighOrder_k2k1.root'))

frame = ROOT.TH1D('hframe', ';#LTd#it{N}_{ch}/d#it{#eta}#GT;#it{#kappa}_{2}(#it{V})/#LT#it{V}#GT^{2}', 1, 0, 21)
frame.GetYaxis().SetTitleOffset(1.2)
frame.GetYaxis().SetRangeUser(0, 0.6)

for i in range(0, 3):
    graph.append(file[i].Get('g_81'))
    graph[i].SetMarkerStyle(20)
    graph[i].SetMarkerSize(1.2)
    graph[i].SetLineColor(color[i])
    graph[i].SetMarkerColor(color[i])

tx = ROOT.TLatex()
tx.SetTextFont(45)
tx.SetTextSize(27)
tx.SetNDC()

leg = ROOT.TLegend(0.5, 0.45, 0.7, 0.7)
leg.SetTextSize(27)
leg.SetTextFont(45)

canv = ROOT.TCanvas('c', 'c', 600, 600)
canv.SetLeftMargin(0.16)
canv.SetTopMargin(0.03)
canv.SetRightMargin(0.03)
canv.SetBottomMargin(0.16)
canv.cd()
frame.Draw()
for g in zip(graph, names):
    g[0].Draw('pesame')
    leg.AddEntry(g[0], g[1], 'pe')
leg.Draw('same')
tx.DrawLatex(0.35, 0.87, 'ALICE Work in progress')
tx.DrawLatex(0.35, 0.8, 'pp #sqrt{#it{s}} = 13 TeV')
tx.DrawLatex(0.35, 0.73, '|#it{#eta}| < 0.8, 0.5 < #it{p}_{T} < 1.5 GeV/#it{c}')
canv.Print('evMix_v2_compare.pdf')
