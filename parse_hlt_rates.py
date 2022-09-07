# Local environment:
# conda activate root.6.26.0
# python metrics.py

################################################################################
# Imports
from __future__ import print_function
from ROOT import TDatime, TGraph, TFile, TH1F, TCanvas, TLegend, gROOT, gStyle, TH2F, kGray,TF1
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label, applyLegendSettings
from officialStyle import officialStyle
from ctypes import c_double
import copy, os, sys
import numpy as np
import json
import array

################################################################################
# Common definitions 
from common import ensureDir, common_path, dr_dict, hackRate, hlt_threshold_dict, fB, Sigma_B, Br_kee

################################################################################
# Plotting
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

################################################################################
# Command line configuration

from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)
parser.add_option('-c', '--corrected', # obsolete now we have "official" rates?
                  action="store_true",
                  default=False,
                  dest='corrected',
                  help="apply trigger rates correction factor")
(options, args) = parser.parse_args()

################################################################################
# 

l1_pts  = np.arange(4.,11.5,0.5).tolist()
hlt_pts = np.arange(4.,11.0,0.5).tolist()

dct = {}
switch_npu = [17,25,30,36,42,48,56]
for npu in switch_npu:
    hlt_file = TFile(common_path+'ee/roc_hlt_pu'+str(npu)+'.root')
    for l1_pt in l1_pts:
        hlt_graph = hlt_file.Get('pt' + str(l1_pt).replace('.','p'))
        hlt_pt = hlt_threshold_dict.get(l1_pt,4.0) # Get HLT pT threshold from "L1 pT <-> HLT pT" map 
        hltpt_list = np.arange(4, 11, 0.5).tolist()
        index = hltpt_list.index(hlt_pt) if hlt_pt in hltpt_list else -1
        name = \
            'L1_DoubleEG'+str(l1_pt).replace('.','p').replace('p0','')+\
            'er1p22_dR_'+str(dr_dict[l1_pt]).replace('.','p')
        title = name
        if name not in dct : dct[name] = []
        if index >= 0 and index < hlt_graph.GetN():
            if hlt_pt != hlt_pts[index] : 
                print("UNEXPECTED HLT PT!!!")
            else:
                hlt_eff = hlt_graph.GetPointX(index)
                hlt_rate = hlt_graph.GetPointY(index)
                dct[name].append((npu,hlt_pt,hlt_rate,hlt_eff))
        else:
            pass
    hlt_file.Close()

output = TFile(common_path+'ee/HLT_rate_DoubleEle.root','RECREATE')
output.cd()
for l1_pt in l1_pts:
    print()
    print("L1 pT:",l1_pt)

    # Name
    name = \
        'L1_DoubleEG'+str(l1_pt).replace('.','p').replace('p0','')+\
        'er1p22_dR_'+str(dr_dict[l1_pt]).replace('.','p')
    print("name: ",name)

    # Print vs NPU (hlt rate evolves, while pT and eff are const)
    for (npu,hlt_pt,hlt_rate,hlt_eff) in dct[name]:
        print("{:30s} {:4.1f}   {:3.1f}   {:2.0f}  {:7.1f}   {:8.6f}".format(name,l1_pt,hlt_pt,npu,hlt_rate,hlt_eff*1.e4))

    # Extract
    [npu,hlt_pt,hlt_rate,hlt_eff] = map(list, zip(*dct[name]))

    # Title
    hlt_pt = hlt_threshold_dict.get(l1_pt,4.0) # Get HLT pT threshold from "L1 pT <-> HLT pT" map 
    title = name+', HLT_DoubleEle'+str(hlt_pt).replace('.','p').replace('p0','')+'_eta1p22_mMax6'
    title = title+', eff='+"{:8.6f}".format(hlt_eff[0]*1.e4)
    print("title:",title)

    # Raw graph
    graph = TGraph(len(switch_npu),
                   array.array('d',npu),
                   array.array('d',hlt_rate))
    graph.SetName(name+'_raw')
    graph.SetTitle(title)
    graph.Write()

    # Fit graph
    fitf = TF1("myfunc","pol4",17,60.)
    graph.Fit(fitf,"Q")

    # Graph from fit
    npus = range(17,61,1)
    rates = [ fitf.Eval(npu) for npu in npus ]
    gr = TGraph(len(npus),
                array.array('d',npus),
                array.array('d',rates))
    gr.SetName(name)
    gr.SetTitle(title)
    gr.Write()

output.Close()
