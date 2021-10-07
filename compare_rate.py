import copy, math, os
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, TCanvas, TColor, kLightTemperature, TGraphErrors
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label
from officialStyle import officialStyle
from array import array
import numpy as np
import itertools

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetTitleOffset(1.0,"X")
gStyle.SetTitleOffset(1.0,"Y")

l1_ptrange = np.arange(5, 12.5, 0.5).tolist() # 6, 11.5, 1.0
hlt_ptrange = np.arange(4, 12.5, 0.5).tolist() # 6, 11.5, 1.0

colours = [1, 2, 4, 6, 8, 13, 15]
styles = [1, 2, 4, 3, 5, 1, 1]


def overflow(hist):
#    import pdb; pdb.set_trace()
    lastp1 = hist.GetBinContent(hist.GetXaxis().GetNbins()+1)
    last = hist.GetBinContent(hist.GetXaxis().GetNbins())
    lastp1e = hist.GetBinError(hist.GetXaxis().GetNbins()+1)
    laste = hist.GetBinError(hist.GetXaxis().GetNbins())
    hist.SetBinContent(hist.GetXaxis().GetNbins(), last+lastp1)
    hist.SetBinError(hist.GetXaxis().GetNbins(), math.sqrt(math.pow(laste,2)+math.pow(lastp1e,2)))
    hist.SetBinContent(hist.GetXaxis().GetNbins()+1, 0)
    hist.SetBinError(hist.GetXaxis().GetNbins()+1, 0)

    firstp1 = hist.GetBinContent(1)
    first = hist.GetBinContent(0)
    firstp1e = hist.GetBinError(1)
    firste = hist.GetBinError(0)
    hist.SetBinContent(1, first+firstp1)
    hist.SetBinError(1, math.sqrt(math.pow(firste,2)+math.pow(firstp1e,2)))



def calc(num, den):
    eff = num/den

    efferr = 0.


    if num!=0:
        efferr = math.sqrt(eff*(1-eff)/num)

#    print num, den, eff, efferr

    return eff, efferr



def set_palette(name="palette", ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.20, 0.2, 0.2, 0.2, 0.2]
#        red   = [0.1, 0.2, 0.3, 0.4, 0.5]
        green = [0.1, 0.1, 0.1, 0.1, 0.1]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]
#        blue  = [0.1, 0.2, 0.3, 0.4, 0.5]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)



def applyHistStyle(h, i):
#    print h, i
    h.SetLineColor(colours[i])
    h.SetMarkerColor(colours[i])
    h.SetMarkerSize(0)
    h.SetLineStyle(styles[i])
    h.SetLineWidth(3)
    h.SetStats(False)

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def comparisonPlots(hists, titles, isLog=False, pname='sync.pdf', isEff = False, isRatio=False, isLegend=False, x=None, y=None, label="test"):

    display = DisplayManager(pname, isLog, isRatio, 0.6, 0.7, x, y, label)
    display.draw_legend = isLegend
    display.isEff = isEff
    
    display.Draw(hists, titles)


def sproducer(key, rootfile, name, ivar, addsel = '1'):

    hist = TH1F('h_' + key + '_' + name, 
                'h_' + key + '_' + name, 
                ivar['nbin'], ivar['xmin'], ivar['xmax'])

    hist.Sumw2()
    exp = '(' + ivar['sel'] + '&&' + addsel + ')'
        
    tree = rootfile.Get(ivar['tree'])

#    print ivar['var'] + ' >> ' + hist.GetName(), exp
    
    tree.Draw(ivar['var'] + ' >> ' + hist.GetName(), exp)
    hist.GetXaxis().SetTitle(ivar['title'])
    hist.GetYaxis().SetTitle('a.u.')
    hist.GetYaxis().SetNdivisions(506)
        
    return copy.deepcopy(hist)


xtit = "Generator-level electron p_{T} [GeV]"
xtit_b = "Generator-level B p_{T} [GeV]"

trackreq = 'hlt_pms2 < 10000'
#trackreq = '1'

vardict = {
    'hlt_pt':{'tree':'tree', 'var':'hlt_pt', 'nbin':30, 'xmin':0, 'xmax':30, 'title':'E/gamma et (GeV)', 'sel':'1'}, 
    'hlt_eta':{'tree':'tree', 'var':'hlt_eta', 'nbin':30, 'xmin':-1.5, 'xmax':1.5, 'title':'E/gamma eta', 'sel':'1'}, 
    'hlt_phi':{'tree':'tree', 'var':'hlt_phi', 'nbin':30, 'xmin':-math.pi, 'xmax':math.pi, 'title':'E/gamma phi', 'sel':'1'}, 
    'hlt_energy':{'tree':'tree', 'var':'hlt_energy', 'nbin':30, 'xmin':0, 'xmax':20, 'title':'E/gamma energy', 'sel':'1'}, 
    'hlt_rawEnergy':{'tree':'tree', 'var':'hlt_rawEnergy', 'nbin':30, 'xmin':0, 'xmax':20, 'title':'E/gamma raw Energy', 'sel':'1'}, 
    'hlt_phiWidth':{'tree':'tree', 'var':'hlt_phiWidth', 'nbin':30, 'xmin':0, 'xmax':0.25, 'title':'E/gamma phiWidth', 'sel':'1'}, 
    'hlt_nrClus':{'tree':'tree', 'var':'hlt_nrClus', 'nbin':10, 'xmin':0, 'xmax':10, 'title':'E/gamma nrClus', 'sel':'1'}, 
    'hlt_sigmalIEtaIEta':{'tree':'tree', 'var':'hlt_sigmaIEtaIEta', 'nbin':30, 'xmin':0, 'xmax':0.025, 'title':'sigmaIEtaIEta', 'sel':'1'}, 
    'hlt_sigmalIEtaIEtaNoise':{'tree':'tree', 'var':'hlt_sigmaIEtaIEtaNoise', 'nbin':30, 'xmin':0, 'xmax':0.025, 'title':'sigmaIEtaIEtaNoise', 'sel':'1'}, 
    'hlt_seednrCrystals':{'tree':'tree', 'var':'hlt_seednrCrystals', 'nbin':15, 'xmin':0, 'xmax':15, 'title':'seednrCrystals', 'sel':'1'}, 
    'hlt_ecalPFIsol':{'tree':'tree', 'var':'hlt_ecalPFIsol', 'nbin':30, 'xmin':0, 'xmax':20, 'title':'ecalPFIsol', 'sel':'1'}, 
    'hlt_hcalPFIsol':{'tree':'tree', 'var':'hlt_hcalPFIsol', 'nbin':30, 'xmin':0, 'xmax':20, 'title':'hcalPFIsol', 'sel':'1'}, 
    'hlt_hcalHForHoverE':{'tree':'tree', 'var':'hlt_hcalHForHoverE', 'nbin':30, 'xmin':0, 'xmax':15, 'title':'hcalHForHoverE', 'sel':'1'}, 
    'hlt_HE':{'tree':'tree', 'var':'hlt_hcalHForHoverE/hlt_energy', 'nbin':30, 'xmin':0, 'xmax':1, 'title':'H/E', 'sel':'1'}, 
    'hlt_pms2':{'tree':'tree', 'var':'hlt_pms2', 'nbin':30, 'xmin':0, 'xmax':10000, 'title':'pms2', 'sel':'1'}, 

    'hlt_dxy':{'tree':'tree', 'var':'hlt_dxy', 'nbin':50, 'xmin':-0.25, 'xmax':0.25, 'title':'dxy', 'sel':trackreq}, 
#    'hlt_trkpt':{'tree':'tree', 'var':'hlt_trkpt', 'nbin':50, 'xmin':-0.25, 'xmax':0.25, 'title':'dxy', 'sel':trackreq}, 
#    'hlt_trketa':{'tree':'tree', 'var':'hlt_', 'nbin':50, 'xmin':-0.25, 'xmax':0.25, 'title':'dxy', 'sel':trackreq}, 
#    'hlt_trkphi':{'tree':'tree', 'var':'hlt_dxy', 'nbin':50, 'xmin':-0.25, 'xmax':0.25, 'title':'dxy', 'sel':trackreq}, 
    'hlt_trkValidHits':{'tree':'tree', 'var':'hlt_trkValidHits', 'nbin':25, 'xmin':0, 'xmax':25, 'title':'trkvalidHits', 'sel':trackreq}, 
    'hlt_trkIsol':{'tree':'tree', 'var':'hlt_trkIsol', 'nbin':30, 'xmin':0, 'xmax':15, 'title':'trk Isol', 'sel':trackreq}, 
    'hlt_trkChi2':{'tree':'tree', 'var':'hlt_trkChi2', 'nbin':30, 'xmin':0, 'xmax':70, 'title':'trk Chi2', 'sel':trackreq}, 
    'hlt_trkMissHits':{'tree':'tree', 'var':'hlt_trkMissHits', 'nbin':3, 'xmin':0, 'xmax':3, 'title':'trackMissHits', 'sel':trackreq}, 
    'hlt_trkNrLayerIT':{'tree':'tree', 'var':'hlt_trkNrLayerIT', 'nbin':5, 'xmin':0, 'xmax':5, 'title':'trkNrLayerIT', 'sel':trackreq}, 
    'hlt_trkDEta':{'tree':'tree', 'var':'hlt_trkDEta', 'nbin':30, 'xmin':0, 'xmax':0.05, 'title':'trkDEta', 'sel':trackreq}, 
    'hlt_trkDEtaSeed':{'tree':'tree', 'var':'hlt_trkDEtaSeed', 'nbin':30, 'xmin':0, 'xmax':0.05, 'title':'trkDEtaSeed', 'sel':trackreq}, 
    'hlt_trkDPhi':{'tree':'tree', 'var':'hlt_trkDPhi', 'nbin':30, 'xmin':0, 'xmax':0.6, 'title':'trkDPhi', 'sel':trackreq}, 
    'hlt_invESeedInvP':{'tree':'tree', 'var':'hlt_invESeedInvP', 'nbin':30, 'xmin':0, 'xmax':0.8, 'title':'invESeedInvP', 'sel':trackreq},
    'hlt_invEInvP':{'tree':'tree', 'var':'hlt_invEInvP', 'nbin':30, 'xmin':0, 'xmax':0.8, 'title':'invEInvP', 'sel':trackreq},

    }

ensureDir('plots_compare/')




file = TFile('ratetest.root')
tree = file.Get('tree')


h_rate = TH2F('rate' , 'rate', len(l1_ptrange)-1, min(l1_ptrange), max(l1_ptrange), len(hlt_ptrange)-1, min(hlt_ptrange), max(hlt_ptrange))

h_rate.GetXaxis().SetTitle('L1 di-electron X GeV')
h_rate.GetYaxis().SetTitle('HLT leading-electron Y GeV')

def calcRate(tree, denstr, numstr):
    den = tree.GetEntries(denstr)
    num = tree.GetEntries(numstr)

    rate = 11200.*2254.*float(num)/float(den)

    print(numstr, '(', num, ')', denstr, '(', den, ')', '---->', rate)
    
    return rate


if True:
    for il1, l1_pt in enumerate(l1_ptrange):
        for ihlt, hlt_pt in enumerate(hlt_ptrange):
        
            sel_den = 'isgjson==1'
            
            sel_l1 = 'l1_doubleE' + str(l1_pt) + '==1'
            sel_hlt = 'doubleE' + str(hlt_pt) + '==1'

            sel_num = '&&'.join([sel_den, sel_l1, sel_hlt])

            h_rate.SetBinContent(il1+1, ihlt+1, calcRate(tree, sel_den, sel_num))


    ofile = TFile('ratemap.root', 'recreate')
    h_rate.Write()

    ofile.Write()
    ofile.Close()

