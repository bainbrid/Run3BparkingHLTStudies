import copy, math, os
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, TCanvas, TColor, kLightTemperature, TGraphErrors, kRed, TGaxis, gPad, kRed, TF1, TLegend
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
#gStyle.SetTitleOffset(1.0,"X")
#gStyle.SetTitleOffset(1.0,"Y")
gStyle.SetPadBottomMargin(0.22)

#l1_ptlist = np.arange(4, 12, 0.5).tolist() 
#hlt_ptlist = np.arange(4, 12, 0.5).tolist() 

const=float(2544.*11200)

colours = [1, 2, 4, 6, 8, 13, 15]
styles = [1, 2, 4, 3, 5, 1, 1]



def calc(num, den):
    eff = num/den

    efferr = 0.


    if num!=0:
        efferr = math.sqrt(eff*(1-eff)/num)

#    print num, den, eff, efferr

    return eff, efferr


def applyLegendSettings2(leg):     
    leg.SetBorderSize(0)           
    leg.SetFillColor(10)           
    leg.SetLineColor(0)            
    leg.SetFillStyle(0)            
    leg.SetTextSize(0.05)          
    leg.SetTextFont(42)     


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

ensureDir('plots/')
ensureDir('root/')

set_palette()


file = TFile(common_path+'ee/hlt_data_for_rate_evaluation.root')
tree = file.Get('tree')

l1_nbin = int((l1_ptmax - l1_ptmin) / l1_ptstep)
hlt_nbin = int((hlt_ptmax - hlt_ptmin) / hlt_ptstep)
#h_rate = TH2F('rate' , 'rate', len(l1_ptlist)-1, min(l1_ptlist), max(l1_ptlist), len(hlt_ptlist)-1, min(hlt_ptlist), max(hlt_ptlist))
h_rate = TH2F('rate' , 'rate', l1_nbin,l1_ptmin,l1_ptmax,hlt_nbin,hlt_ptmin,hlt_ptmax)

h_rate.GetXaxis().SetTitle('L1 di-electron X GeV')
h_rate.GetYaxis().SetTitle('HLT leading-electron Y GeV')


#h_rate_mass = TH2F('rate_mass' , 'rate_mass', len(l1_ptlist)-1, min(l1_ptlist), max(l1_ptlist), len(hlt_ptlist)-1, min(hlt_ptlist), max(hlt_ptlist))
h_rate_mass = TH2F('rate_mass' , 'rate_mass', l1_nbin,l1_ptmin,l1_ptmax,hlt_nbin,hlt_ptmin,hlt_ptmax)

h_rate_mass.GetXaxis().SetTitle('L1 di-electron X GeV')
h_rate_mass.GetYaxis().SetTitle('HLT leading-electron Y GeV')


def calcRate(tree, denstr, numstr):
    den = tree.GetEntries(denstr)
    num = tree.GetEntries(numstr)

    rate = 11200.*2254.*float(num)/float(den)

    print(numstr, '(', num, ')', denstr, '(', den, ')', '---->', rate)
    
    return rate


def effproducer(tree, hname, sel):

    print('producing ...', hname)

#    hist = TH2F(hname, hname, 60,0,60, 2,-0.5,1.5)
    hist = TH2F(hname, hname, 30,0,60, 2,-0.5,1.5)

    tree.Draw(sel + ':npu' +  '>> ' + hname)

    print(hname, '->', tree.GetEntries(sel), '/', tree.GetEntries())

#    import pdb; pdb.set_trace()

#    for ii in range(100000):
#        hist.Fill(0,0)

    hprof = hist.ProfileX()

    # this is a bit dirty but ok ... 
#    hprof.SetBinContent(1, 0.0000000039)
#    hprof.SetBinError(1, 0.00000001)



    hprof.Scale(const)
    hprof.GetXaxis().SetTitle('# of PU / inst. L (E34)')
    hprof.GetXaxis().SetTitleOffset(2.1)
    hprof.GetYaxis().SetTitleOffset(1.3)
    hprof.GetYaxis().SetMaxDigits(3)
    hprof.GetYaxis().SetTitle('HLT trigger rate (Hz)')
    hprof.GetYaxis().SetNdivisions(506)
    hprof.SetTitle(hname)
    hprof.SetName(hname)

#    print(hprof.GetBinContent(1))


    can = TCanvas('can_' + hname, 'can_' + hname)
    can.SetGridx()
    can.SetGridy()

    hprof.Draw()
    
#    fitfunc = TF1('parabolla', '[0]*x*x*x + [1] *x*x + [2]*x', 0,60)

#    fitfunc = TF1('expo_mod', 'exp([0] + [1]*x)', 0, 60)
#    fitfunc = TF1('pol3_mod', '[0]*x*x*x + [1]*x*x + [2]*x', 0, 60)
    fitfunc = TF1('pol3_mod', '[0]*x*x*x*x*x + [1]*x*x*x*x + [2]*x*x*x + [3]*x*x + [4]*x', 0, 60)
#    fitfunc.SetParLimits(1, 0, 15)
#    fitfunc.SetParLimits(0, -2, 2)
#    fitfunc.FixParameter(0, -2.3)

    fitfunc.SetParLimits(0, 0, 100)
    fitfunc.SetParLimits(1, 0, 100)
    fitfunc.SetParLimits(2, 0, 100)
    fitfunc.SetParLimits(3, 0, 100)
    fitfunc.SetParLimits(4, 0, 100)

#    hprof.Fit("expo_mod")
#    hprof.Fit("expo")
    hprof.Fit("pol3_mod")
#    hprof.Fit("parabolla")
    
#    func = hprof.GetFunction("pol3")
#    func = hprof.GetFunction("expo_mod")
    func = hprof.GetFunction("pol3_mod")
#    func = hprof.GetFunction("expo")
    print('fit param. = ',  func.GetParameter(0), func.GetParameter(1), func.GetParameter(2))

#    func = hprof.GetFunction("parabolla")
    func.SetLineColor(kRed)
    func.SetLineStyle(2)
    func.SetLineWidth(3)
    func.Draw('same')
    func.SetTitle(hname + '_fit')
    func.SetName(hname + '_fit')

#    print 'test1', gPad.GetUxmin(), gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUymax()
#    print 'test2', Double(60)*0.0357338

#    axis = TGaxis(gPad.GetUxmin(), gPad.GetUymin(),
#                  gPad.GetUxmax(), gPad.GetUymin(), 0, Double(60)*0.0357338, 510,"+L");

    delta = hprof.GetMaximum() - hprof.GetMinimum()
    y = hprof.GetMinimum() - 0.15*delta

    print( hprof.GetMaximum(), hprof.GetMinimum(), 'delta=', delta, 'y=', y)
    axis = TGaxis(0, y, 60, y, 0, 60.*0.0357338-0.0011904, 506, "+L")
    axis.SetLabelFont(42)
    axis.SetLabelSize(0.05)
#    gPad.GetUxmax(), gPad.GetUymin(), 0, Double(60)*0.0357338, 510,"+L")
#    axis.SetLineColor(kRed)
#    axis->SetLabelColor(kRed)
    axis.Draw();


    leg = TLegend(0.19, 0.81,0.3,0.86)
    applyLegendSettings2(leg)
    leg.SetTextSize(0.032)
    leg.SetTextFont(62)
    leg.AddEntry(hprof,   'L1: p_{T} > ' + hname.split('_')[1].replace('p','.') + 'GeV, HLT: p_{T} > ' + hname.split('_')[3].replace('p','.') + 'GeV', '')
    leg.Draw()

#    leg2 = TLegend(0.19, 0.76,0.3,0.8)
#    applyLegendSettings2(leg2)
#    leg2.SetTextSize(0.032)
#    leg2.SetTextFont(62)
#    leg2.AddEntry(hprof,   'exp(' + '{0:.2f}'.format(func.GetParameter(0)) + '+' + '{0:.2f}'.format(func.GetParameter(1)) + '*x)', '')
#    leg2.AddEntry(hprof,   '{0:.2f}'.format(func.GetParameter(0)) + '*x^{3} + ' + '{0:.2f}'.format(func.GetParameter(1)) + '*x^{2}' + '{0:.2f}'.format(func.GetParameter(1)) + '*x', '')
#    leg2.Draw()

    can.SaveAs('plots/' + hname + '.pdf')
    can.SaveAs('plots/' + hname + '.gif')

    return copy.deepcopy(hprof), copy.deepcopy(hist), copy.deepcopy(func)




hists = []

graphs = []

iloop = 0
nloop = len(l1_ptlist)*len(hlt_ptlist) if working_points==False else len(l1_ptlist)
start = time.time()
for il1, l1_pt in enumerate(l1_ptlist):
    for ihlt, hlt_pt in enumerate(hlt_ptlist):
        if working_points==True and il1 != ihlt : continue
        print("".join(["Processing #",str(iloop)," out of ",str(nloop)," triggers..."]))
        
        sel_den = 'isgjson==1'
            
        sel_l1 = 'l1_doubleE' + str(l1_pt).replace('.','p') + '==1'
        sel_hlt = 'doubleE' + str(hlt_pt).replace('.','p') + '==1'
        sel_hlt_mass = 'mass_doubleE' + str(hlt_pt).replace('.','p') + '==1'

        sel_num = '&&'.join([sel_den, sel_l1, sel_hlt])
        sel_num_mass = '&&'.join([sel_den, sel_l1, sel_hlt_mass])

        xbin = h_rate.GetXaxis().FindBin(l1_pt+1.e-3)
        ybin = h_rate.GetYaxis().FindBin(hlt_pt+1.e-3)
        h_rate.SetBinContent(xbin, ybin, calcRate(tree, sel_den, sel_num))
        xbin = h_rate_mass.GetXaxis().FindBin(l1_pt+1.e-3)
        ybin = h_rate_mass.GetYaxis().FindBin(hlt_pt+1.e-3)
        h_rate_mass.SetBinContent(xbin, ybin, calcRate(tree, sel_den, sel_num_mass))

        name = 'l1_' + str(l1_pt).replace('.','p') + '_hlt_'+ str(hlt_pt).replace('.', 'p')

        hprof, hist, func = effproducer(tree, name, sel_num_mass)

        graphs.append(hprof)
        graphs.append(func)

        timing(iloop,nloop,start,"loops")
        iloop+=1

            #####
            #####

###            sel_dep = '&&'.join([sel_l1, sel_hlt])
###
###            hname = 'l1_' + str(l1_pt) + '_hlt_' + str(hlt_pt)
###            hist = TH2F(hname, hname, 60,0,2, 2,-0.5,1.5)
###
###            tree.Draw(sel_dep + ':instL' +  '>> ' + hname, sel_den)
###
####            print hname, '->', tree.GetEntries(sel), '/', tree.GetEntries()
###
###            hprof = hist.ProfileX()
###            hprof.Scale(const)
###            #    hprof.GetXaxis().SetTitle('# of PU')
###            #    hprof.GetXaxis().SetTitleOffset(1.5)
###            hprof.GetYaxis().SetTitle('Trigger Rate')
###            hprof.GetYaxis().SetNdivisions(506)
###            hprof.SetTitle(hname)
###            hprof.SetName(hname)
###
###            can = TCanvas('can_' + hname, 'can_' + hname)
###            can.SetGridx()
###            can.SetGridy()
###            
###            hprof.Draw()
###            hprof.Fit("pol3")
###            
###            func = hprof.GetFunction("pol3")
###            func.SetLineColor(kRed)
###            func.SetLineStyle(2)
###            func.SetLineWidth(3)
###            func.Draw('same')
###            func.SetTitle(hname + '_fit')
###            
###            func.SetName(hname + '_fit')
###            can.SaveAs('plots/' + hname + '.pdf')
###            can.SaveAs('plots/' + hname + '.gif')
###
###            hists.append((hprof, func))
            
#            return copy.deepcopy(hprof), copy.deepcopy(hist), copy.deepcopy(func)





ofile = TFile(common_path+'ee/ratemap4roc.root', 'recreate')
h_rate.Write()
h_rate_mass.Write()

for graph in graphs:
    graph.Write()

ofile.Write()
ofile.Close()

