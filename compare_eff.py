import copy, math, os
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, TCanvas, TColor, kLightTemperature, TGraphErrors
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label
from officialStyle import officialStyle
from array import array
import numpy as np
import itertools
from common import common_path
from common import l1_ptmin, l1_ptmax, l1_ptstep, l1_ptlist
from common import hlt_ptmin, hlt_ptmax, hlt_ptstep, hlt_ptlist
from common import working_points
from common import timing
import time

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetTitleOffset(1.0,"X")
gStyle.SetTitleOffset(1.0,"Y")

from optparse import OptionParser, OptionValueError
usage = "usage: python compare_eff.py"
parser = OptionParser(usage)
parser.add_option('-w', '--weight', action="store_true", default=False, dest='weight')
parser.add_option('-p', '--plot', action="store_true", default=False, dest='plot')
parser.add_option("-o", "--option", action="store", type="int", default=0, dest="option", help="option: 0 = original, 1 = same as 0, but reprocessed, 2 = trg with q2 req, 3 = use gen q2 eff map")
(options, args) = parser.parse_args()

# Analysis efficiency map
eff_histo = None
eff_histo2 = None
if options.weight:
    if options.option == 0: # original
        eff_file = TFile(common_path+'eff_maps/0_original/eff.root')
        eff_histo = eff_file.Get('eff_pt1_vs_pt2_tot_weighted')
    elif options.option == 1: # reprocessed
        eff_file = TFile(common_path+'eff_maps/1_reproduced/eff.root')
        eff_histo = eff_file.Get('eff_pt1_vs_pt2_tot_weighted')
    else: # 2 or 3 ("trg with q2 req" or "use gen q2 eff map")
        eff_file = TFile(common_path+'eff_maps/2_add_qsq_req/eff.root')
        eff_histo = eff_file.Get('eff_pt1_vs_pt2_tot_weighted') # eff with reco q2 req. in numer
        #eff_histo = eff_file.Get('eff_pt1_vs_pt2_tot_norecoqsq') # eff with no reco q2 req. in numer
        if options.option == 3: # gen_q2_eff_map
            eff_file2 = TFile(common_path+'eff_maps/1_reproduced/eff2.root')
            eff_histo2 = eff_file2.Get('eff_pt1_vs_pt2_gen_qsq')
    print('Found analysis efficiencies!',eff_histo)

#l1_ptlist = np.arange(4, 11, 0.5).tolist()
#hlt_ptlist = np.arange(4, 11, 0.5).tolist()

colours = [1, 2, 4, 6, 8, 13, 15]
styles = [1, 2, 4, 3, 5, 1, 1]

drdict = {
    3.0:1.0,
    3.5:1.0,
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



set_palette()

file = TFile(common_path+'ee/gen_for_efficiency_evaluation.root')
tree = file.Get('tree')


#h_e1 = TH2F('e1' , 'e1', len(l1_ptlist)-1, min(l1_ptlist), max(l1_ptlist), len(hlt_ptlist)-1, min(hlt_ptlist), max(hlt_ptlist))
#
#h_e2 = TH2F('e2' , 'e2', len(l1_ptlist)-1, min(l1_ptlist), max(l1_ptlist), len(hlt_ptlist)-1, min(hlt_ptlist), max(hlt_ptlist))
#
#h_e1e2 = TH2F('e1e2' , 'e1e2', len(l1_ptlist)-1, min(l1_ptlist), max(l1_ptlist), len(hlt_ptlist)-1, min(hlt_ptlist), max(hlt_ptlist))
#
#h_all = TH2F('all' , 'all', len(l1_ptlist)-1, min(l1_ptlist), max(l1_ptlist), len(hlt_ptlist)-1, min(hlt_ptlist), max(hlt_ptlist))

l1_nbin = int((l1_ptmax - l1_ptmin) / l1_ptstep)
hlt_nbin = int((hlt_ptmax - hlt_ptmin) / hlt_ptstep)
h_gall = TH2F('gall' , 'gall', l1_nbin,l1_ptmin,l1_ptmax,hlt_nbin,hlt_ptmin,hlt_ptmax)
#h_gall_mass = TH2F('gall_mass' , 'gall_mass', len(l1_ptlist)-1, min(l1_ptlist), max(l1_ptlist), len(hlt_ptlist)-1, min(hlt_ptlist), max(hlt_ptlist))


#h_e1.GetXaxis().SetTitle('L1 di-electron X GeV')
#h_e1.GetYaxis().SetTitle('HLT leading-electron Y GeV')
#
#h_e2.GetXaxis().SetTitle('L1 di-electron X GeV')
#h_e2.GetYaxis().SetTitle('HLT subleading-electron Y GeV')
#
#h_e1e2.GetXaxis().SetTitle('L1 di-electron X GeV')
#h_e1e2.GetYaxis().SetTitle('HLT subleading-electron Y GeV')
#
#h_all.GetXaxis().SetTitle('L1 di-electron X GeV')
#h_all.GetYaxis().SetTitle('HLT di-electron Y GeV')

h_gall.GetYaxis().SetTitle('HLT di-electron X GeV')
h_gall.GetYaxis().SetTitle('HLT di-electron Y GeV')

#h_gall_mass.GetYaxis().SetTitle('HLT di-electron X GeV')
#h_gall_mass.GetYaxis().SetTitle('HLT di-electron Y GeV')

def calcEff(tree, denstr, numstr):
    # Denominator
    den = tree.GetEntries(denstr)
    # Numerator
    #gen_q2 = 'gen_mass > 1.05 && gen_mass < 2.5'
    gen_q2 = 'gen_mass*gen_mass > 1.1 && gen_mass*gen_mass < 6.25' if options.option == 2 else "1"
    if not options.weight :
        num = tree.GetEntries(numstr + " && " + gen_q2)
    else :
        gen_pt  = 'gen_e1_pt > 2.0 && gen_e2_pt > 2.0'
        gen_eta = 'abs(gen_e1_eta) < 1.2 && abs(gen_e2_eta) < 1.2'
        num = 0.
        ybins = eff_histo.GetYaxis().GetNbins() #ybins = eff_histo2.GetYaxis().GetNbins()
        xbins = eff_histo.GetXaxis().GetNbins() #xbins = eff_histo2.GetXaxis().GetNbins()
        cntr=0
        unweighted=0
        for ybin in range(1,ybins+1):
            for xbin in range(1,xbins+1):
                if xbin > ybin : continue
                cntr+=1
                print("".join(["  ","Processing #",str(cntr)," out of ",str(sum([x+1 for x in range(xbins)]))," eff bins..."]),)
                eff = eff_histo.GetBinContent(xbin,ybin) #eff = eff_histo.GetBinContent(min(xbin,6),min(ybin,6))
                #print("eff",eff)
                eff2 = eff_histo2.GetBinContent(xbin,ybin) if options.option == 3 else 1.
                y_down = eff_histo.GetYaxis().GetBinLowEdge(ybin)
                y_up = eff_histo.GetYaxis().GetBinLowEdge(ybin) + eff_histo.GetYaxis().GetBinWidth(ybin)
                x_down = eff_histo.GetXaxis().GetBinLowEdge(xbin)
                x_up = eff_histo.GetXaxis().GetBinLowEdge(xbin) + eff_histo.GetXaxis().GetBinWidth(xbin)
                gen_e1_pt = 'gen_e1_pt > ' + str(y_down)
                if ybin < ybins : gen_e1_pt += ' && ' + 'gen_e1_pt < ' + str(y_up)
                gen_e2_pt = 'gen_e2_pt > ' + str(x_down)
                if xbin < xbins : gen_e2_pt += ' && ' + 'gen_e2_pt < ' + str(x_up)
                newcuts = [gen_e1_pt, gen_e2_pt, gen_pt, gen_eta]
                if options.option == 2: newcuts.append(gen_q2)
                newgencut = ' && '.join(newcuts)
                entry = float(tree.GetEntries(numstr + ' && ' + newgencut))
                eff=1.
                num += entry*eff*eff2
                unweighted += entry
                #print("   ybin",ybin,"xbin",xbin)
                #print("   e1_pt_bin",y_down,"e2_pt_bin",x_down)
                #print("   eff",eff,)
                #print("   counts",entry)
                #print("   unweighted",unweighted)
                #print("   num",num)

    # Efficiency
    eff = float(num)/float(den)
    #eff_unw = float(unweighted)/float(den)
    #print("denom:",den,"unweighted:",unweighted,"weighted:",num,"eff(unw):",eff_unw,"eff:",eff)
    print(numstr, '(', num, ')', denstr, '(', den, ')', '---->', eff)
    return eff

qcut = 'e1_hlt_pms2 < 10000 && e1_hlt_invEInvP < 0.2 && e1_hlt_trkDEtaSeed < 0.01 && e1_hlt_trkDPhi < 0.2 && e1_hlt_trkChi2 < 40 && e1_hlt_trkValidHits >= 5 && e1_hlt_trkNrLayerIT >= 2'

iloop = 0
nloop = len(l1_ptlist)*len(hlt_ptlist) if working_points==False else len(l1_ptlist)
start = time.time()
for il1, l1_pt in enumerate(l1_ptlist):
    for ihlt, hlt_pt in enumerate(hlt_ptlist):
        if working_points==True and il1 != ihlt : continue
        print("".join(["Processing #",str(iloop)," out of ",str(nloop)," triggers..."]))

        sel_inclusive = '1'
        sel_den = 'gen_e1_l1_dr < 0.2 && gen_e2_l1_dr < 0.2 && e1_l1_pt >= ' + str(l1_pt) + ' && e2_l1_pt >= ' + str(l1_pt) 
        sel_dr = 'l1_eedr < ' + str(drdict[l1_pt]) + ' && hlt_eedr < ' + str(drdict[hlt_pt])
        sel_dr_mass = 'l1_eedr < ' + str(drdict[l1_pt]) + ' && hlt_mee < 6'
        
        match_e1 = 'gen_e1_hlt_dr < 0.2 && e1_hlt_pt >= ' + str(hlt_pt)
        match_e2 = 'gen_e2_hlt_dr < 0.2 && e2_hlt_pt >= ' + str(hlt_pt)

        sel_e1 = '&&'.join([sel_den, match_e1])
        sel_e2 = '&&'.join([sel_den, match_e2])
        
        sel_e1e2 = '&&'.join([sel_den, match_e1, match_e2])
        sel_all = '&&'.join([sel_den, match_e1, match_e2, sel_dr, qcut, qcut.replace('e1', 'e2')])
        
        sel_all_mass = '&&'.join([sel_den, match_e1, match_e2, sel_dr_mass, qcut, qcut.replace('e1', 'e2')])

        xbin = h_gall.GetXaxis().FindBin(l1_pt+1.e-3)
        ybin = h_gall.GetYaxis().FindBin(hlt_pt+1.e-3)
        #h_e1.SetBinContent(il1+1, ihlt+1, calcEff(tree, sel_den, sel_e1))
        #h_e2.SetBinContent(il1+1, ihlt+1, calcEff(tree, sel_den, sel_e2))
        #h_e1e2.SetBinContent(il1+1, ihlt+1, calcEff(tree, sel_den, sel_e1e2))
        #h_all.SetBinContent(il1+1, ihlt+1, calcEff(tree, sel_den, sel_all_mass))
        h_gall.SetBinContent(xbin, ybin, calcEff(tree, sel_inclusive, sel_all_mass))
        #h_gall_mass.SetBinContent(il1+1, ihlt+1, calcEff(tree, sel_inclusive, sel_all_mass))

        timing(iloop,nloop,start,"loops")
        iloop+=1

ensureDir('root/')

if not options.weight: ofile = TFile(common_path+'ee/effmap4roc.root', 'recreate')
else:                  ofile = TFile(common_path+'ee/effmap4roc_weighted.root', 'recreate')
#    h_e1.Write()
#    h_e2.Write()
#    h_e1e2.Write()
#    h_all.Write()
h_gall.Write()
#    h_gall_mass.Write()

ofile.Write()
ofile.Close()

        
    

filedict = {'sig_e1': {'file':file, 'sel':'gen_e1_l1_dr < 0.2 && gen_e2_l1_dr < 0.2 && e1_l1_pt >= 6 && e2_l1_pt >= 6 && gen_e1_hlt_dr < 0.2'},
            'sig_e2': {'file':file, 'sel':'gen_e1_l1_dr < 0.2 && gen_e2_l1_dr < 0.2 && e1_l1_pt >= 6 && e2_l1_pt >= 6 && gen_e1_hlt_dr < 0.2 && gen_e2_hlt_dr < 0.2'},
            'data': {'file':TFile(common_path+'ee/hlt_data_perele_dist.root'), 'sel':'isgjson==1 && l1_doubleE6==1'}}


if options.plot:
    ensureDir(common_path+'plots/')

    for vkey, ivar in vardict.items():

        print(vkey)

        hists = []
        titles = []

        for stype in ['data', 'sig_e1', 'sig_e2']:
        

            ivar_new = copy.deepcopy(ivar)
            if stype.find('sig')!=-1:
                ivar_new['var'] = ivar_new['var'].replace('hlt_', stype.replace('sig_','') + '_hlt_')
                ivar_new['sel'] = ivar_new['sel'].replace('hlt_', stype.replace('sig_','') + '_hlt_')
            
            print(vkey, stype, ivar_new['var'], ivar_new['sel'])

            hist = sproducer(stype, filedict[stype]['file'], vkey, ivar_new, filedict[stype]['sel'])
    
            hists.append(copy.deepcopy(hist))
            titles.append(stype)


        for ii, ihist in enumerate(hists):
            applyHistStyle(ihist, ii)
            overflow(ihist)
            ihist.Scale(1./ihist.GetSumOfWeights())
            ihist.SetMaximum(ihist.GetBinContent(ihist.GetMaximumBin())*1.2)

        comparisonPlots(hists, titles, False, common_path+'plots/' + vkey + '.pdf', False, False, True)
        comparisonPlots(hists, titles, True, common_path+'plots/' + vkey + '_log.pdf', False, False, True)
