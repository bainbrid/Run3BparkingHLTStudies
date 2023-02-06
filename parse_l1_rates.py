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

# Prescale columns identified by "peak Linst" values
switch_cols = [ ( float("{:.1f}".format(x/10.)),
                  float("{:.1f}".format((x-1)/10.)) ) for x in range(25,5,-1) ]

# Corresponding NPU used to extract L1 total and di-ele rates and HLT rates
# Fit: nPU = 32.0[e-34] * Linst + 0.9, established from data with 2160b
# Assuming gradient of 33, scaled to 2748b, 33*(2160/2748) = 26
# Hence, at nominal operating parameters, 2E34 @ 2748b, NPU = 52 ...
# ... which is in line with expectations
switch_npu = [ float("{:.1f}".format(x[0]*26.)) for x in switch_cols ]

file_old = TFile(common_path+'ee/l1_bandwidth_official.root')
file_new = TFile(common_path+'ee/L1_rate_DoubleEG.root')
# Iterate through L1 pT thresholds
for ii,l1_pt in enumerate(l1_pts) :
    name = \
        'L1_DoubleEG'+str(l1_pt).replace('.','p').replace('p0','')+\
        'er1p22_dR_' + str(dr_dict[l1_pt]).replace('.','p')
    fit_old = file_old.Get(name)
    fit_new = file_new.Get(name)
    ratios = []
    for npu in switch_npu:
        value_old = fit_old.Eval(npu) * 1000. # rates given in kHz, determined with 2544b
        value_new = fit_new.Eval(npu) * 2544. # rates given in Hz/nbunches, scaled up to 2544b
        ratio = value_new / value_old
        ratios.append(ratio)
    dct[name] = ratios
file_old.Close()
file_new.Close()

output = TFile(common_path+'ee/L1_rate_DoubleEG_corr.root','RECREATE')
for name,ratios in dct.items():
    graph = TGraph(len(switch_npu),
                   array.array('d',switch_npu),
                   array.array('d',ratios))
    graph.SetName(name)
    graph.SetTitle(name)
    graph.Write()
output.Close()
