import copy, math, os
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, TCanvas, TLegend, TGraph
from officialStyle import officialStyle
import numpy as np

l1_ptrange = np.arange(6, 11.5, 1).tolist() # 6, 11.5, 1
hlt_ptrange = np.arange(4, 11.5, 1).tolist() # 6, 11.5, 1

print('l1', l1_ptrange)
print('hlt', hlt_ptrange)

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)
parser.add_option('-w', '--weight', action="store_true", default=False, dest='weight')
parser.add_option("-l", "--lumi", default=1.0, type="float", help="target lumi. with [E34]", dest="lumi")
(options, args) = parser.parse_args()


def returnGraph(name, rates, effs):
    graph = TGraph()
    graph.SetName('pt' + str(name).replace('.0', ''))
    graph.SetTitle('L1 di-e p_{T} > ' + str(name).replace('.0','') + 'GeV')


    idx = 0
    for rate, eff in zip(rates, effs):

        graph.SetPoint(idx, eff, rate)

        idx += 1

    return copy.deepcopy(graph)


colours = [1, 2, 4, 6, 8, 13, 15]
styles = [1, 2, 4, 3, 5, 1, 1,1]


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def applyLegendSettings(leg):
    leg.SetBorderSize(0)
    leg.SetFillColor(10)
    leg.SetLineColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.05)
#    leg.SetTextFont(42)


file_rate = TFile('ratemap.root')
file_eff = None
if not options.weight: file_eff = TFile('effmap.root')
else: file_eff = TFile('effmap_weighted.root')

ratemap = file_rate.Get('rate')
effmap = file_eff.Get('gall')

print(type(ratemap), type(effmap))

graphs = []

for l1pt in l1_ptrange:

    rates = []
    effs = []
    
    for hltpt in hlt_ptrange:


        print(l1pt, hltpt, type(l1pt), type(hltpt))
        xbin = ratemap.GetXaxis().FindBin(l1pt)
        ybin = ratemap.GetYaxis().FindBin(hltpt)

        rate = ratemap.GetBinContent(xbin, ybin)
        eff = effmap.GetBinContent(xbin, ybin)

        rates.append(rate)
        effs.append(eff)

    graph = returnGraph(l1pt, rates, effs)
    graphs.append(copy.deepcopy(graph))



canvas = TCanvas()
canvas.SetLogy()
canvas.SetLogx()
canvas.SetGridx()
canvas.SetGridy()

frame_roc = TH2F('frame', 'frame', 100, 0.000003, 0.01, 1000, 1, 100000)

if not options.weight: frame_roc.GetXaxis().SetTitle('L1 x HLT Trigger eff.')
else: frame_roc.GetXaxis().SetTitle('L1 x HLT x Analysis eff.')
frame_roc.GetYaxis().SetTitle('Rate @ PU ~ 40')
frame_roc.Draw()

leg = TLegend(0.2, 0.64,0.5,0.86)

applyLegendSettings(leg)
leg.SetTextSize(0.03)


for idx, graph in enumerate(graphs):
    print(idx)
    graph.SetLineColor(colours[idx])
    graph.SetMarkerColor(colours[idx])
    graph.SetMarkerSize(1)
    graph.Draw('plsame')
    leg.AddEntry(graph, graph.GetTitle(), 'l')

    
leg.Draw()
if not options.weight:
    canvas.SaveAs('roc_hlt.gif')
    canvas.SaveAs('roc_hlt.pdf')
else:
    canvas.SaveAs('roc_hlt_weighted.gif')
    canvas.SaveAs('roc_hlt_weighted.pdf')
