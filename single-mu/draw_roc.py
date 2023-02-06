import copy, math, os
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, TCanvas, TLegend, TGraph, TF1
from officialStyle import officialStyle
from DisplayManager import add_CMS, add_Preliminary
import numpy as np
from common import path

l1_ptrange = np.arange(6, 11.5, 1).tolist() 
hlt_ptrange = np.arange(4, 11.5, 1).tolist() 

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
parser.add_option('-p', '--pr', action="store_true", default=False, dest='pr')
(options, args) = parser.parse_args()

# Analysis efficiency map
eff_histos = {}
if options.pr and options.weight:
    eff_file = TFile(path+'eff_maps/2022Jan19/eff.root')
    #eff_histos['Mu7_IP4_match']  = eff_file.Get('eff_pt1_vs_pt2_qsq7_weighted')
    #eff_histos['Mu8_IP5_match']  = eff_file.Get('eff_pt1_vs_pt2_qsq8_weighted') #@@ this eff is for Mu8_IP3!!
    eff_histos['Mu9_IP6_match']  = eff_file.Get('tot16_pt1_vs_pt2_qsq9_weighted')
    #eff_histos['Mu12_IP6_match'] = eff_file.Get('eff_pt1_vs_pt2_qsq12_weighted')
    print('Found analysis efficiencies!',eff_histos)

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



def calc(name, num, den):
    eff = float(num)/float(den)
    efferr = 0.
    if num!=0:
        efferr = math.sqrt(eff*(1-eff)/float(den))
    print(name.rjust(15), " numer: ",num, " denom: ", den, " eff: {0:.8f}".format(eff), " err: {0:.8f}".format(efferr))
    #print(name.replace('_', '\_').rjust(15), " & {0:.5f}".format(eff), " $\\pm$ {0:.5f}".format(efferr))
    return eff, efferr



hltmenus = {
#    'Mu12_IP6':{'xmin':20, 'xmax':40},
#    'Mu9_IP5':{'xmin':16, 'xmax':30},
#    'Mu7_IP4':{'xmin':16, 'xmax':24},
#    'Mu8_IP5':{'xmin':16, 'xmax':24},
#    'Mu8_IP6':{'xmin':16, 'xmax':24},
#    'Mu9_IP6':{'xmin':20, 'xmax':30},
#
#    'Mu12_IP6_match':{'xmin':20, 'xmax':40},
#    'Mu9_IP5_match':{'xmin':16, 'xmax':30},
#    'Mu7_IP4_match':{'xmin':16, 'xmax':24},
#    'Mu8_IP5_match':{'xmin':16, 'xmax':24},
#    'Mu8_IP6_match':{'xmin':16, 'xmax':24},
    'Mu9_IP6_match':{'xmin':20, 'xmax':30},
}

ensureDir('plots')


file = TFile(path+'single-mu/gen_for_singlemu.root')
tree = file.Get('tree')

ratefile = TFile(path+'single-mu/obs_rate_summary.root')

save2file = []

# eff. evaluation 
gen_pt  = 'gen_e1_pt > 0.5 && gen_e2_pt > 0.5' 
gen_eta = 'abs(gen_e1_eta) < 2.4 && abs(gen_e2_eta) < 2.4' 
gen_q2  = '(gen_mass*gen_mass) > 1.1 && (gen_mass*gen_mass) < 6.25'



for hlt, hltdict in hltmenus.items():

    if options.pr:
        # Denominator
        den = tree.GetEntries()
        # Numerator
        num = 0.
        numstr = hlt+'==1'
        if not options.weight :
            num = tree.GetEntries(numstr)
        else :
            gen_pt  = 'gen_e1_pt > 2.0 && gen_e2_pt > 2.0'
            gen_eta = 'abs(gen_e1_eta) < 2.5 && abs(gen_e2_eta) < 2.5'
            eff_histo = eff_histos.get(hlt,None)
            if eff_histo is not None :
                print(eff_histo)
                ybins = eff_histo.GetYaxis().GetNbins()
                xbins = eff_histo.GetXaxis().GetNbins()
                cntr=0
                for ybin in range(1,ybins+1):
                    for xbin in range(1,xbins+1):
                        if xbin > ybin : continue
                        print("".join(["  ","#",str(cntr)," out of ",str(ybins*xbins/2.)]))
                        cntr+=1
                        eff = eff_histo.GetBinContent(xbin,ybin)
                        y_down = eff_histo.GetYaxis().GetBinLowEdge(ybin)
                        y_up = eff_histo.GetYaxis().GetBinLowEdge(ybin) + eff_histo.GetYaxis().GetBinWidth(ybin)
                        x_down = eff_histo.GetXaxis().GetBinLowEdge(xbin)
                        x_up = eff_histo.GetXaxis().GetBinLowEdge(xbin) + eff_histo.GetXaxis().GetBinWidth(xbin)
                        gen_e1_pt = 'gen_e1_pt > ' + str(y_down)
                        if ybin < ybins : gen_e1_pt += ' && ' + 'gen_e1_pt < ' + str(y_up)
                        gen_e2_pt = 'gen_e2_pt > ' + str(x_down)
                        if xbin < xbins : gen_e2_pt += ' && ' + 'gen_e2_pt < ' + str(x_up)
                        newgencut = ' && '.join([gen_e1_pt, gen_e2_pt, gen_pt, gen_eta])
                        entry = float(tree.GetEntries(numstr + ' && ' + newgencut))
                        num += entry*eff
            
        calc(hlt, num, den)
#        continue

    
    if hlt.find('match')!=-1: continue 

    # read rate 

    graph = ratefile.Get('HLT_' + hlt + '_part0')

    fit = TF1('fit', '[0]*x + [1]')

    graph.Fit(fit, '', '', hltdict['xmin'], hltdict['xmax'])

    fitted = graph.GetFunction('fit')
    
    p1 = fitted.GetParameter(0)
    p2 = fitted.GetParameter(1)

    fitted = TF1('fitted', str(p1) + '*x + ' + str(p2), hltdict['xmin'], hltdict['xmax'])
    fitted.SetLineWidth(2)
    fitted.SetLineColor(2)

    fitted_all = TF1('fitted', str(p1) + '*x + ' + str(p2), 0, 60)
    fitted_all.SetLineStyle(2)
    fitted_all.SetLineColor(1)
    fitted_all.SetName(hlt)
    fitted_all.SetTitle(hlt)

    canvas = TCanvas()
    
    frame = TH2F('frame_' + hlt, 'frame_' + hlt, 60,0,60,100,0,3)
    frame.GetXaxis().SetTitle('PU')
    frame.GetYaxis().SetTitle('unprescaled rate / # of bunches (Hz)')
    frame.GetXaxis().SetNdivisions(506)
    frame.GetYaxis().SetNdivisions(506)
    frame.Draw()

    graph.Draw('psame')

    fitted.Draw('same')
    fitted_all.Draw('same')

    l2=add_Preliminary()
    l2.Draw("same")
    l3=add_CMS()
    l3.Draw("same")

    leg = TLegend(0.15, 0.79, 0.5, 0.89)
    applyLegendSettings(leg)
    leg.AddEntry(fitted_all, hlt, '')
    leg.Draw()

    canvas.SaveAs('plots/' + hlt + '.gif')
    canvas.SaveAs('plots/' + hlt + '.pdf')

    save2file.append(fitted_all)
    


#
#file_rate = TFile('ratemap.root')
#file_eff = TFile('effmap.root')
#
#ratemap = file_rate.Get('rate')
#effmap = file_eff.Get('gall')
#
#print(type(ratemap), type(effmap))
#
#graphs = []
#
#for l1pt in l1_ptrange:
#
#    rates = []
#    effs = []
#    
#    for hltpt in hlt_ptrange:
#
#
#        print(l1pt, hltpt, type(l1pt), type(hltpt))
#        xbin = ratemap.GetXaxis().FindBin(l1pt)
#        ybin = ratemap.GetYaxis().FindBin(hltpt)
#
#        rate = ratemap.GetBinContent(xbin, ybin)
#        eff = effmap.GetBinContent(xbin, ybin)
#
#        rates.append(rate)
#        effs.append(eff)
#
#    graph = returnGraph(l1pt, rates, effs)
#    graphs.append(copy.deepcopy(graph))
#
#
#
#canvas = TCanvas()
#canvas.SetLogy()
#canvas.SetLogx()
#canvas.SetGridx()
#canvas.SetGridy()
#
#frame_roc = TH2F('frame', 'frame', 100, 0.000003, 0.01, 1000, 1, 100000)
#
#frame_roc.GetXaxis().SetTitle('L1 x HLT Trigger eff.')
#frame_roc.GetYaxis().SetTitle('Rate @ PU ~ 40')
#frame_roc.Draw()
#
#leg = TLegend(0.2, 0.64,0.5,0.86)
#
#applyLegendSettings(leg)
#leg.SetTextSize(0.03)
#
#
#for idx, graph in enumerate(graphs):
#    print(idx)
#    graph.SetLineColor(colours[idx])
#    graph.SetMarkerColor(colours[idx])
#    graph.SetMarkerSize(1)
#    graph.Draw('plsame')
#    leg.AddEntry(graph, graph.GetTitle(), 'l')
#
#    
#leg.Draw()
#canvas.SaveAs('roc_hlt.gif')
#canvas.SaveAs('roc_hlt.pdf')

if len(save2file)!=0:
    ensureDir('root')
    out = TFile('../root/single-mu/obs_rate_summary_fit.root', 'recreate')
    for fit in save2file:
        fit.Write()

    out.Write()
    out.Close()
