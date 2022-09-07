# Local environment:
# conda activate root.6.26.0
# python metrics.py

################################################################################
# Imports
from __future__ import print_function
from ROOT import TDatime, TGraph, TFile, TH1F, TCanvas, TLegend, gROOT, gStyle, TH2F, kGray
from DisplayManager import DisplayManager, add_Preliminary, add_CMS, add_label, applyLegendSettings
from officialStyle import officialStyle
from ctypes import c_double
import copy, os, sys
import numpy as np
import json

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

parser.add_option('-l', '--limit', # obsolete?
                  action="store_true",
                  default=False,
                  dest='limit',
                  help="limit HLT rate to X Hz") # either 2018 bandwidth (from Sara) or 300 Hz
(options, args) = parser.parse_args()

################################################################################
# Configuration

# Output dir for plots
ensureDir('plots/')

# Level-1 total rate estimate for full CMS menu
l1_file = TFile(common_path+'ee/l1_bandwidth.root')
l1_his = l1_file.Get('otherrate') # CMS L1 rate vs Linst

# Di-electron trigger rate estimates from data 
# Contains parameterised L1 di-electron rate vs nPU
ee_file = TFile(common_path+'ee/l1_bandwidth_official.root')

# Number of fills to consider 
nfills = 100

# Normalise to unit integrated luminosity
integrated_lumi = 1.

# L1 and HLT pT thresholds
l1_pts = np.arange(4.,11.0,0.5).tolist()
hlt_pts = np.arange(4.,11.0,0.5).tolist()

# List of L_inst values to consider
# Linst = nPU*0.0357338 - 0.0011904
switch_lumi = [#(2.2, 2.0),
               (2.0, 1.7), (1.7, 1.5), (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.6), (0.6, 0.)]

# List of nPU values to consider
# nPU = (Linst+0.0011904)/0.0357388
switch_npu = [#62,
              56, 48, 42, 36, 30, 25, 17]

# Maximum L1 trigger bandwidth
l1_max = 95000.

# Rate dedicated to new di-muon seeds in Run 3
dimuon = 0. # Redundant?

# Determine correction to L1 rate for total CMS menu
# i.e. ensure rate does not exceed l1_max @ 2E34 (i.e. remove any overspend)
idx_2e34 = 0
l1_rate_2e34 = l1_his.Eval(switch_lumi[idx_2e34][0])
l1_rate_corr = max(0.,l1_rate_2e34-l1_max) # correct overspend only

# Corrections to HLT rate, now redundant?
corrs_dict = {}
if options.corrected:
    for npu in [56, 48, 42, 36, 30, 25, 17]:
        filename = common_path+'rates/corrections_' + str(npu) + '.json'
        infile = open(filename,'r')
        dct = json.load(infile)
        corrs_dict[npu] = dct
        infile.close()

################################################################################
# Print various tables

def print_table(debug,value,l1_max=95000.,dimuon=0.,allocation=0.):
    print("Table:",value)

    # Open HLT files for different nPU values
    hlt_files = [ TFile(common_path+'ee/roc_hlt_pu'+str(npu)+'.root') for npu in switch_npu ] 

    # Store variables to build "menu"
    menu_dict = {}

    # Iterate through L1 pT thresholds
    for ii,l1_pt in enumerate(l1_pts) :

        # Get correct histogram for di-electron L1 trigger rate, parameterised vs nPU
        histo_ee_rate = ee_file.Get('L1_DoubleEG'+
                                    str(l1_pt).replace('.','p').replace('p0','')+
                                    'er1p22_dR_' + str(dr_dict[l1_pt]).replace('.','p'))

        if debug : print("L1 LOOP:",ii,l1_pt,histo_ee_rate) # DEBUG
        if not debug and not "menu" in value : print("{:4.1f}".format(l1_pt),end="")

        # Iterate through L_inst values (and correponding nPU value and HLT file)
        for jj,(lumi,npu,hlt_file) in enumerate(reversed(list(zip(switch_lumi,switch_npu,hlt_files)))) :
            
            peak_lumi = lumi[0] # Upper bound from switch_lumi (list of 2-tuple)

            # Evaluate and correct L1 rate for CMS menu based on peak L_inst
            l1_rate = l1_his.Eval(peak_lumi)
            l1_rate -= l1_rate_corr # remove any overspend @ 2E34 (i.e. "translate down")
            l1_rate = hackRate(l1_rate,peak_lumi) # <<< HACK HACK HACK

            # Determine L1 spare capacity 
            spare = l1_max - l1_rate # APPLY CMS-WIDE PRESCALE HERE (IF REQUIRED)!

            # Scale (i.e. "decay") di-electron rate allocation according to L_inst
            alloc = allocation*(l1_rate/(l1_rate_2e34-l1_rate_corr))
            
            # Determine rate of L1 di-electron trigger for given nPU
            ee_rate = histo_ee_rate.Eval(npu)*1000

            ##### MISSING???
            #ee_rate = hackRate(ee_rate,which_lumi) # <<< HACK HACK HACK 
            ##### MISSING???

            # Check if the L1 di-ele rate satisfies the dedicate rate allocation or spare capacity
            l1_ok = ee_rate < (alloc+1.e-6) or ee_rate < (spare+alloc+1.e-6) # spare + epsilon
                
            # Get appropriate "graph" for HLT rate (dependant on L1 pT and L_inst)
            hlt_roc = hlt_file.Get('inv_pt' + str(l1_pt).replace('.','p'))
            hlt_n = hlt_roc.GetN()

            if options.limit:

                # Extract HLT pT, rate, and efficiency if limited to 300 Hz (obsolete?)

                hlt_max = 300. #max_bw_hlt[peak_lumi] # <-- Limit to 300 Hz for prompt reco
                hlt_ok = True
                hlt_rate = -1
                hlt_eff = -1
                hlt_pt = -1
                for kk in range(hlt_n):
                    ip = kk#hlt_n - kk - 1
                    # print("test",kk,ip,hlt_pts[ip],hlt_roc.GetPointX(ip),hlt_max,hlt_roc.GetPointY(ip))
                    hlt_rate = hlt_roc.GetPointX(ip)
                    hlt_eff = hlt_roc.GetPointY(ip)
                    hlt_pt = hlt_pts[ip]
                    if options.corrected:
                        l1ptstr = str(l1_pt)
                        hltptstr = str(hlt_pt)
                        if npu in corrs_dict.keys() and \
                                l1ptstr in corrs_dict[npu].keys() and \
                                hltptstr in corrs_dict[npu][l1ptstr].keys():
                            corr = corrs_dict[npu][l1ptstr][hltptstr]
                            if corr is not None: hlt_rate *= corr
                            else: print("null value",npu,l1ptstr,hltptstr)
                        else: print("unknown thresholds",npu,l1ptstr,hltptstr)
                    if hlt_rate < hlt_max:
                        if debug : print("    HLT SCAN:",hlt_n,ip,hlt_pt,hlt_max,hlt_rate,hlt_eff)
                        break
                if hlt_eff == -1: hlt_ok = False #print('!!!! This cannot happen !!!!')

            else:

                # Extract HLT pT, rate, and efficiency

                hlt_ok = True
                hlt_rate = -1
                hlt_eff = -1
                hlt_pt = hlt_threshold_dict.get(l1_pt,4.0) # Get HLT pT threshold from "L1 pT <-> HLT pT" map 
                hltpt_list = np.arange(4, 11, 0.5).tolist()
                index = hltpt_list.index(hlt_pt)
                if index < hlt_roc.GetN():
                    hlt_rate = hlt_roc.GetPointX(index) # Get HLT rate for given L1 pT, HLT pT, and L_inst
                    hlt_eff = hlt_roc.GetPointY(index)  # Get HLT signal efficiency (L1+HLT+)
                    if options.corrected: # If corrections are needed for the HLT rates? Now redundant?
                        l1ptstr = str(l1_pt)
                        hltptstr = str(hlt_pt)
                        if npu in corrs_dict.keys() and \
                                l1ptstr in corrs_dict[npu].keys() and \
                                hltptstr in  corrs_dict[npu][l1ptstr].keys():
                            corr = corrs_dict[npu][l1ptstr][hltptstr]
                            if corr is not None: hlt_rate *= corr
                            else: print("null value",npu,l1ptstr,hltptstr)
                        else: print("unknown thresholds",npu,l1ptstr,hltptstr)
                else: print("cannot find hltpt")

            # Store various metrics to later build menu
            if peak_lumi not in menu_dict and l1_ok and hlt_ok : 
                menu_dict[peak_lumi] = (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,spare+alloc)

            # Print appropriate table given "value" argument to method 
            if not debug : 
                if value == "ee_rate" : 
                    string = " & {:5.1f}".format(ee_rate/1000.)
                    if ee_rate > l1_max : string = " &  -   "
                    if not l1_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "l1_rate" : 
                    string = " & {:5.1f}".format(l1_rate/1000.)
                    print(string,end="")
                elif value == "spare" : 
                    string = " & {:5.1f}".format(spare/1000.)
                    print(string,end="")
                elif value == "l1_pt" : 
                    string = " & {:5.1f}".format(l1_pt)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "hlt_pt" : 
                    string = " & {:5.1f}".format(hlt_pt)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "hlt_rate" : 
                    string = " & {:5.2f}".format(hlt_rate/1000.)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "eff" : 
                    string = " & {:5.2f}".format(hlt_eff*1.e4)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "eff_per_rate" :
                    string = " & {:5.3f}".format((hlt_eff*1.e4)/(ee_rate/1000.))
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "lint" : 
                    Lint = peak_lumi * 1.e-5 * 3600. * 12. # Assume 12-hour fill
                    string = " & {:5.2f}".format(Lint)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "counts_per_fill" : 
                    Lint = peak_lumi * 1.e-5 * 3600. * 12.
                    count = Lint * fB * Sigma_B * Br_kee * hlt_eff
                    string = " & {:5.1f}".format(count)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "counts_per_nfills" : 
                    Lint = peak_lumi * 1.e-5 * 3600. * 12.
                    count = Lint * fB * Sigma_B * Br_kee * hlt_eff * nfills
                    string = " & {:5.1f}".format(count)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif value == "counts_per_fb" : 
                    Lint = integrated_lumi # defined above (25.?)
                    count = Lint * fB * Sigma_B * Br_kee * hlt_eff
                    string = " & {:5.1f}".format(count)
                    if ee_rate > l1_max : string = " &  -   "
                    elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
                    print(string,end="")
                elif "menu" in value : 
                    pass
                else : 
                    print("???")

            # 
            if debug : print("  LUMI LOOP:",jj,peak_lumi,npu,\
               "L1:","{:.1f}".format(spare),"{:.1f}".format(ee_rate),\
               "HLT:",hlt_n,hlt_pt,"{:.1f}".format(hlt_rate),"{:.2f}".format(hlt_eff*1.e4),\
               hlt_ok)#,hlt_roc, # DEBUG

        if not debug and not "menu" in value : print("\\\\")

    # End of loop through L1 pT thresholds ...

    # Special case: print "menu"
    if value == "menu" :
        Lint = integrated_lumi # defined above
        #print("Lint:", Lint)
        #print("Nfills:", nfills)
        print("Linst   L1 pT   HLT pT    L1 rate HLT rate"
              "      AxE   L/fill   #/fill   #/Lint  AxE/kHz")  #/fills")
        for jj,lumi in enumerate(switch_lumi) :
            peak_lumi = lumi[0]
            if peak_lumi not in menu_dict : 
                print(peak_lumi," \\\\")
                continue
            (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,capacity) = menu_dict[peak_lumi]
            count = Lint * fB * Sigma_B * Br_kee * hlt_eff
            Lint_per_fill = peak_lumi * 1.e-5 * 3600. * 12.
            count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff
            hlt_eff_per_rate = (hlt_eff*1.e4) / (ee_rate/1000.)
            print(peak_lumi," & ",
                  "{:5.1f}".format(l1_pt),
                  " & {:5.1f}".format(hlt_pt),
                  " & & {:5.1f}".format(ee_rate/1000.),
                  " & {:5.2f}".format(hlt_rate/1000.),
                  " & {:5.2f}".format(hlt_eff*1.e4),
                  " & {:5.2f}".format(Lint_per_fill),
                  " & {:5.2f}".format(count_per_fill),
                  " & {:5.1f}".format(count),
                  " & {:5.3f}".format(hlt_eff_per_rate),
                  #" & {:5.1f}".format(count_per_fill*nfills),
                  " \\\\ ")

    # Special case: print "menu" that makes use of prescales
    elif value == "menu_prescaled" :
        peak_lumi_ = 0.6 # this is the default trigger
        (l1_pt_,hlt_pt_,ee_rate_,hlt_rate_,hlt_eff_,capacity_) = menu_dict[peak_lumi_]

        Lint = integrated_lumi # defined above
        #print("Lint:", Lint)
        #print("Nfills:", nfills)
        print("Linst   L1 pT   HLT pT    L1 rate HLT rate"
              "      AxE   L/fill   #/fill     #/fb  AxE/kHz prescale")
        for jj,(lumi,npu) in enumerate(zip(switch_lumi,switch_npu)) :
            peak_lumi = lumi[0]
            if peak_lumi not in menu_dict : 
                print(peak_lumi," \\\\")
                continue
            (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,capacity) = menu_dict[peak_lumi]
            #@@prescale = ee_rate_ / ee_rate if ee_rate > 0. else 1.e9

            histo_ee_rate = ee_file.Get('L1_DoubleEG'+
                                        str(l1_pt_).replace('.','p').replace('p0','')+
                                        'er1p22_dR_' + str(dr_dict[l1_pt_]).replace('.','p'))
            ee_rate_ = histo_ee_rate.Eval(npu)*1000

            prescale = ee_rate_/capacity if capacity > 0. else 1.
            prescale = max(prescale,1.) # Ensure prescale >= 1.

            count = Lint * fB * Sigma_B * Br_kee * (hlt_eff_/prescale)
            Lint_per_fill = peak_lumi * 1.e-5 * 3600. * 12.
            count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * (hlt_eff_/prescale)
            hlt_eff_per_rate = (hlt_eff_*1.e4/prescale) / (ee_rate_/1000./prescale)
            print(peak_lumi," & ",
                  "{:5.1f}".format(l1_pt_),
                  " & {:5.1f}".format(hlt_pt_),
                  " & & {:5.1f}".format(ee_rate_/1000./prescale),
                  " & {:5.2f}".format(hlt_rate_/1000./prescale),
                  " & {:5.2f}".format(hlt_eff_*1.e4/prescale),
                  " & {:5.2f}".format(Lint_per_fill),
                  " & {:5.2f}".format(count_per_fill),
                  " & {:5.1f}".format(count),
                  " & {:5.3f}".format(hlt_eff_per_rate),
                  " & {:5.2f}".format(prescale),
                  " \\\\ ")

    print()

################################################################################
# Print various tables here
if __name__ == "__main__":

    # Print individual metrics
    for value in ["ee_rate",
                  "l1_rate",
                  "spare",
                  "l1_pt",
                  "hlt_pt",
                  "hlt_rate",
                  "eff",
                  "eff_per_rate",
                  "lint",
                  "counts_per_fill",
                  "counts_per_fb",
                  #"counts_per_nfills",
                  ]: 
        print_table(debug=False,value=value,l1_max=l1_max,dimuon=dimuon,allocation=0.) # Assume 0 kHz allocation here

    # Print "menu" for different allociations (and with/out prescales)
    for allocation in [0.,5000.,10000.,20000.,90000]:
        print("allocation:",allocation)
        print_table(debug=False,value="menu",l1_max=l1_max,dimuon=dimuon,allocation=allocation)
        print("Prescales based on trigger used @ 0.6E34...!")
        print_table(debug=False,value="menu_prescaled",l1_max=l1_max,dimuon=dimuon,allocation=allocation)
