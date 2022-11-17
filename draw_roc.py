import copy, math, os, sys
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, TCanvas, TLegend, TGraph
from officialStyle import officialStyle
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label, applyLegendSettings, applyLegendSettings2
import numpy as np
from common import common_path, hlt_threshold_dict
import json

#l1_ptrange = np.arange(5, 10.9, 1.0).tolist() 
#hlt_ptrange = np.arange(4, 10.9, 1.0).tolist() 

l1_ptrange = np.arange(4, 11.5, 0.5).tolist() 
hlt_ptrange = np.arange(4, 11.0, 0.5).tolist() 

print('l1', l1_ptrange)
print('hlt', hlt_ptrange)

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)
parser.add_option('-w', '--weight', action="store_true", default=True, dest='weight', help="weight by analysis eff")
parser.add_option('-r', '--rates', action="store_true", default=False, dest='rates', help="write trigger rates to JSON")
parser.add_option('-c', '--corrected', action="store_true", default=False, dest='corrected', help="apply trigger rates correction factor")
(options, args) = parser.parse_args()


drdict = {
    4.0:0.9,
    4.5:0.9,
    5.0:0.9,
    5.5:0.8,
    6.0:0.8,
    6.5:0.8,
    7.0:0.8,
    7.5:0.7,
    8.0:0.7,
    8.5:0.7,
    9.0:0.7,
    9.5:0.6,
    10.0:0.6,
    10.5:0.6,
    11.0:0.6,
    11.5:0.5,
    12.0:0.5,
    12.5:0.5,
    13.0:0.5,
    13.5:0.4,
    14.0:0.4,
}   


effrefs = {
    #'Mu12_IP6':0.00011,
    'Mu9_IP6':0.00000371 if options.weight else 0.00025312
    #'Mu9_IP5':0.00028,
    #'Mu8_IP5':0.00040,
    #'Mu7_IP4':0.00066
    }

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

# 12, 10, 9, 8, 7
colours_l1 = [1, 12,  9, 7, 5]
styles = [1, 2, 4, 3, 5, 1, 1,1]


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)





def createPdf(rootfile, pname, xt, yt, prec):
    canvas = TCanvas(pname)

    canvas.cd()
    histo = rootfile.Get(pname)

    histo.Draw("TEXT COL")
    histo.GetXaxis().SetTitle(xt) #"sub-leading gen. electron p_{T} (GeV)")
    histo.GetYaxis().SetTitle(yt) #"leading gen. electron p_{T} (GeV)")
    gStyle.SetPaintTextFormat('.' + str(prec) + 'f');


    histo.SetStats(0)
    histo.SetMarkerSize(0.7)

    l2=add_Preliminary()
    l2.Draw("same")
    l3=add_CMS()
    l3.Draw("same")

    canvas.SaveAs("plots/hlt_" + pname + ".pdf")
    canvas.SaveAs("plots/hlt_" + pname + ".gif")

    canvas.SetLogz()
    canvas.SaveAs("plots/hlt_" + pname + "_log.pdf")
    canvas.SaveAs("plots/hlt_" + pname + "_log.gif")

    return canvas



def createROCPdf(effmap, l1_file_rate, file_rate, file_ref, npu, name):

    graphs = []
    graph_envelope = TGraph()
    graph_envelope.SetName('envelope')
    graph_envelope.SetTitle('envelope')
    graph_envelope_inv = TGraph()
    graph_envelope_inv.SetName('envelope_inv')
    graph_envelope_inv.SetTitle('envelope_inv')

    ii = 0
    
    rates_dict = {} # built separately for each npu (per method call)

    corrs_dict = {}
    if options.corrected:
        filename = common_path+'rates/corrections_' + str(npu) + '.json'
        infile = open(filename,'r')
        corrs_dict = json.load(infile)
        infile.close()
        #print(corrs_dict)
        #quit()

    for l1pt in l1_ptrange:

        rates = []
        effs = []
    
        for hltpt in hlt_ptrange:


            #print(l1pt, hltpt, type(l1pt), type(hltpt))
            xbin = effmap.GetXaxis().FindBin(l1pt)
            ybin = effmap.GetYaxis().FindBin(hltpt)
            
#            rate = ratemap.GetBinContent(xbin, ybin)
            #print('file_rate:',file_rate)
            #print('getting ... l1_' + str(l1pt).replace('.','p') + '_hlt' + str(hltpt).replace('.','p') + '_fit')
            rate = file_rate.Get('l1_' + str(l1pt).replace('.','p') + '_hlt_' + str(hltpt).replace('.','p') + '_fit').Eval(npu)
            eff = effmap.GetBinContent(xbin, ybin)
#            print('rate, eff=', rate, eff)

            # Record rates to write to JSON
            if l1pt not in rates_dict.keys() : rates_dict[l1pt] = {}
            rates_dict[l1pt][hltpt] = rate

            # Extract corrections to rates from JSON
            corr = 1.
            if options.corrected:
                l1ptstr = str(l1pt)
                hltptstr = str(hltpt)
                if l1ptstr in corrs_dict.keys() and hltptstr in corrs_dict[l1ptstr].keys():
                    corr = corrs_dict[l1ptstr][hltptstr]
                    if corr is not None: rate *= corr
                    else: print("null value",l1ptstr,hltptstr)
                else: print("unknown thresholds",l1ptstr,hltptstr)

            rates.append(rate)
            effs.append(eff)

            ##### 
            #if hltpt == l1pt - 1.:
            if hltpt == hlt_threshold_dict.get(l1pt,4.0):

#                print('(l1, hlt) = ', l1pt, hltpt)

                graph_rep = TGraph()
                graph_rep.SetName('rep_l1pt' + str(l1pt).replace('.','p') + '_hltpt' + str(hltpt).replace('.','p'))
                graph_rep.SetTitle('rep_l1pt' + str(l1pt).replace('.','p') + '_hltpt' + str(hltpt).replace('.','p'))
                
                graph_rep.SetPoint(0, eff, rate)

                graph_rep.SetMarkerColor(1)
                graph_rep.SetMarkerSize(4)
                graph_rep.SetMarkerStyle(42)

                graphs.append(copy.deepcopy(graph_rep))

                graph_envelope.SetPoint(ii, eff, rate)
                graph_envelope_inv.SetPoint(ii, rate, eff)
                ii += 1


        graph = returnGraph(l1pt, rates, effs)
        graph.SetMarkerSize(1)
        graph.SetName('pt' + str(l1pt).replace('.','p'))


#        import pdb; pdb.set_trace()
#        print('doubleE' + str(l1pt) + ', dR < ' + str(drdict[l1pt]) + '0')
#        l1_rate_file = l1_file_rate.Get('doubleE' + str(l1pt) + ', dR < ' + str(drdict[l1pt]) + '0')
#        l1_rate = l1_rate_file.Eval(npu*0.0357338 - 0.0011904)
#        l1_rate *= 0.001

        #print('l1_file_rate:',l1_file_rate)
        #print('getting:','L1_DoubleEG' + str(l1pt).replace('.','p').replace('p0','') + 'er1p22_dR_' + str(drdict[l1pt]).replace('.','p'))
        l1_rate_file = l1_file_rate.Get('L1_DoubleEG' + str(l1pt).replace('.','p').replace('p0','') + 'er1p22_dR_' + str(drdict[l1pt]).replace('.','p'))
        l1_rate = l1_rate_file.Eval(npu)

        graph.SetTitle('pt' + str(l1pt).replace('.','p') + ' (' + '{0:.1f}'.format(l1_rate) + 'kHz)')

        
        graph_inv = returnGraph(l1pt, effs, rates)
        graph_inv.SetName('inv_pt' + str(l1pt).replace('.','p'))
        graph_inv.SetTitle('inv_pt' + str(l1pt).replace('.','p'))

        graphs.append(copy.deepcopy(graph))
        graphs.append(copy.deepcopy(graph_inv))

    graphs.append(copy.deepcopy(graph_envelope))
    graphs.append(copy.deepcopy(graph_envelope_inv))

    graphs_ref = []

    if options.rates:
        filename = common_path+'rates/rates_' + str(npu) + '.json'
        outfile = open(filename,'w')
        json_string = json.dump(rates_dict,outfile)
        #outfile.write(json_string)
        outfile.close()

    for iref, refname in enumerate(['Mu9_IP6']):#['Mu12_IP6', 'Mu9_IP6', 'Mu9_IP5', 'Mu8_IP5', 'Mu7_IP4']):

        graph_ref = TGraph()
        graph_ref.SetName('ref_' + refname)
        graph_ref.SetTitle('ref_' + refname)

        rate_ref = file_ref.Get(refname)

#        print(refname, 0, effrefs[refname], rate_ref.Eval(options.pu)*2544.)
        graph_ref.SetPoint(0, effrefs[refname], rate_ref.Eval(npu)*2544.)
        if iref==0:
            graph_ref.SetMarkerStyle(30)
        else:
            graph_ref.SetMarkerStyle(29)

        graph_ref.SetMarkerColor(colours_l1[iref])
        graph_ref.SetMarkerSize(3)
        graphs_ref.append(graph_ref)


    # no envelope
    makeCanvas(name, options.weight, False, graphs, graphs_ref, npu)

    # with envelope
    makeCanvas(name, options.weight, True, graphs, graphs_ref, npu)



def makeCanvas(name, weight, envelope, graphs, graphs_ref, npu):

    canvas = TCanvas(name)
    canvas.SetLogy()
    canvas.SetLogx()
    canvas.SetGridx()
    canvas.SetGridy()

    if not weight: 
        frame_roc = TH2F(name, name, 100, 3.e-6, 3.e-3, 1000, 0.1, 100000)
        frame_roc.GetXaxis().SetTitle('Trigger efficiency (L1 and HLT)')
    else : 
        frame_roc = TH2F(name, name, 100, 3.e-7, 3.e-4, 1000, 0.1, 100000)
        frame_roc.GetXaxis().SetTitle('Signal A#times#varepsilon (L1, HLT, Reco, Analysis)')


    frame_roc.GetYaxis().SetTitle('HLT Trigger Rate (Hz)')
    frame_roc.Draw()

    if envelope == True: leg = TLegend(0.18, 0.15,0.4,0.25)
    else:                leg = TLegend(0.18, 0.15,0.4,0.65)

    applyLegendSettings(leg)
    leg.SetTextSize(0.03)

    
    if not envelope:

        col_ = 1

        for idx, graph in enumerate(graphs):
            if graph.GetName().find('rep')==-1 and graph.GetName().find('envelope')==-1 and graph.GetName().find('inv')==-1:

                col = col_
                if col >=10:
                    col += 1

                    
                graph.SetLineColor(col)
                graph.SetMarkerColor(col)
                leg.AddEntry(graph, 'L1 p_{T} > ' + graph.GetTitle().replace('pt','').replace('p','.'), 'l')

                col_ += 1

                graph.Draw('plsame')

#        leg.Draw()

    else:

        for idx, graph in enumerate(graphs):
            if graph.GetName()!='envelope': continue
            
#            print(idx, graph.GetName())

            graph.SetLineStyle(0)
            graph.SetLineWidth(0)
            graph.SetMarkerColor(1)
            graph.SetMarkerSize(1)
            graph.SetMarkerStyle(20)
            graph.Draw('plsame')
            leg.AddEntry(graph, 'Pairwise (L1,HLT) thresholds', 'l')

    leg.Draw()

    leg2 = TLegend(0.19, 0.81,0.3,0.85)
    applyLegendSettings2(leg2)
    leg2.SetTextSize(0.032)
    leg2.SetTextFont(62)
    leg2.AddEntry(frame_roc, 'L = {0:.1f}'.format(float(npu*0.0357338 - 0.0011904)) + 'E34 (PU ~ ' + str(npu) + ')', '')
    leg2.Draw()


    for graph_ in graphs_ref:
        graph_.Draw('psame')
        leg.AddEntry(graph_, graph_.GetName().replace('ref_',''), 'p')


    l2=add_Preliminary()
    l2.Draw("same")
    l3=add_CMS()
    l3.Draw("same")


    suffix1 = '_envelope' if envelope else ''
    suffix2 = '_weighted' if weight else ''

    canvas.SaveAs('plots/' + name + '_' + str(npu).replace('.','p') +suffix1+suffix2+ '.gif')
    canvas.SaveAs('plots/' + name + '_' + str(npu).replace('.','p') +suffix1+suffix2+ '.pdf')
        

    if not envelope:
        file = TFile(common_path + 'ee/' + name + '_pu' + str(int(npu)) + '.root', 'recreate')

        for idx, graph in enumerate(graphs):
            graph.Write()

        file.Write()
        file.Close()

#    for idx, graph in enumerate(graphs_inv):
#        graph.Write()

        
#    graph_envelope.Write()
#    graph_envelope_inv.Write()





ensureDir('plots')
ensureDir('root')

file_rate = TFile(common_path+'ee/ratemap4roc.root')

#createPdf(file_rate, 'rate', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 0)
#createPdf(file_rate, 'rate_mass', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 0)

file_eff = None

if not options.weight: file_eff = TFile(common_path+'ee/effmap4roc.root')
else: file_eff = TFile(common_path+'ee/effmap4roc_weighted.root')

effmap = file_eff.Get('gall')

file_ref = TFile(common_path+'single-mu/obs_rate_summary_fit.root')
#ref = file_ref.Get('Mu9_IP6')

l1_file_rate = TFile(common_path+'ee/l1_bandwidth_official.root')

#createPdf(file_eff, 'e1', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
#createPdf(file_eff, 'e2', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
#createPdf(file_eff, 'e1e2', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
#createPdf(file_eff, 'all', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)', 5)
#createPdf(file_eff, 'gall', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)',5 )
#createPdf(file_eff, 'gall_mass', 'Level-1 di-e X (GeV)', 'HLT di-e Y (GeV)',5 )

#sys.exit(1)

for npu in [56, 48, 42, 36, 30, 25, 17]:

    createROCPdf(effmap, l1_file_rate, file_rate, file_ref, npu, 'roc_hlt')

#createROCPdf(effmap, file_rate, file_ref, 'roc_mass_hlt')


