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
_l1_file = TFile(common_path+'ee/l1_bandwidth.root')
_l1_his = _l1_file.Get('otherrate') # CMS L1 rate vs Linst

# Di-electron L1 trigger rate estimates from data 
# Contains parameterised L1 di-electron rate vs nPU
_ee_file = TFile(common_path+'ee/l1_bandwidth_official.root')

# Normalise to integrated luminosity
_integrated_lumi = 86.4 # 12-hour fill @ 2E34 delivers 0.864/fb 

# L1 and HLT pT thresholds
_l1_pts = np.arange(4.,11.5,0.5).tolist()
_hlt_pts = np.arange(4.,11.0,0.5).tolist()

# Prescale columns labelled by peak Linst values
_switch_cols = [#(2.2, 2.0),
                (2.0, 1.7), 
                (1.7, 1.5),
                #(1.7, 1.6), (1.6, 1.5),
                (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.6), (0.6, 0.)]

# List of nPU values to consider
# nPU = (Linst+0.0011904)/0.0357388
_switch_npu = [#62,
               56, 48, 42, #@@ @ 1.6E34?,
               36, 30, 25, 17]
# OVERWRITE
#switch_npu = [ float("{:.1f}".format((lumi[0]+0.0011904)/0.0357388)) for lumi in switch_lumi ]
#switch_npu = [ float("{:.1f}".format(lumi[0]*(56.366/2.))) for lumi in switch_lumi ]
#switch_npu = [ float("{:.1f}".format(lumi[0]*(56./2.))) for lumi in switch_lumi ]

# List of L_inst values to consider
# Linst = nPU*0.0357338 - 0.0011904
_switch_lumi = [#(2.2, 2.0),
                (2.0, 1.7), 
                (1.7, 1.5),
                #(1.7, 1.6), (1.6, 1.5),
                (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.6), (0.6, 0.)]

print(_switch_cols)
print(_switch_npu)
print(_switch_lumi)
#quit()

# Maximum L1 trigger bandwidth
_l1_max_2018 = 95000.
_l1_max = 100000.

# Nominal number of bunches in LHC orbit
_nb_default = 2740.

# Determine correction to L1 rate for total CMS menu
# i.e. ensure rate does not exceed l1_max @ 2E34 (i.e. remove any overspend)
_idx_2e34 = 0
_l1_rate_2e34 = _l1_his.Eval(_switch_lumi[_idx_2e34][0])
_l1_rate_corr = max(0.,_l1_rate_2e34-_l1_max_2018) # correct overspend only

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
# Prints the appropriate "metric table" (according to the 'table' arg)
# Called within a nested loop over L1 pT and peak Linst (and/or NPU)
# i.e. prints only a row for the given L1 pT
def print_metric_table_row(table, # prints this metric
                           l1_pt, # prints row per L1 pT
                           ee_rate,l1_rate,spare,hlt_pt,hlt_rate,hlt_eff, # metrics
                           l1_ok,hlt_ok,l1_max, # conditions
                           peak_lumi # for Lint calc
                           ):
    if table == "ee_rate" : 
        string = " & {:5.1f}".format(ee_rate/1000.)
        if ee_rate > l1_max : string = " &  -   "
        if not l1_ok : string = string.replace(" & "," & \gr{")+"}"
        print(string,end="")
    elif table == "l1_rate" : 
        string = " & {:5.1f}".format(l1_rate/1000.)
        print(string,end="")
    elif table == "spare" : 
        string = " & {:5.1f}".format(spare/1000.)
        print(string,end="")
    elif table == "l1_pt" : 
        string = " & {:5.1f}".format(l1_pt)
        if ee_rate > l1_max : string = " &  -   "
        elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
        print(string,end="")
    elif table == "hlt_pt" : 
        string = " & {:5.1f}".format(hlt_pt)
        if ee_rate > l1_max : string = " &  -   "
        elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
        print(string,end="")
    elif table == "hlt_rate" : 
        string = " & {:5.2f}".format(hlt_rate/1000.)
        if ee_rate > l1_max : string = " &  -   "
        elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
        print(string,end="")
    elif table == "eff" : 
        string = " & {:5.2f}".format(hlt_eff*1.e4)
        if ee_rate > l1_max : string = " &  -   "
        elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
        print(string,end="")
    elif table == "eff_per_rate" :
        string = " & {:5.3f}".format((hlt_eff*1.e4)/(ee_rate/1000.))
        if ee_rate > l1_max : string = " &  -   "
        elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
        print(string,end="")
    elif table == "lint" : 
        Lint = peak_lumi * 1.e-5 * 3600. * 12. # Assume 12-hour fill
        string = " & {:5.2f}".format(Lint)
        if ee_rate > l1_max : string = " &  -   "
        elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
        print(string,end="")
    elif table == "counts_per_fill" : 
        Lint = peak_lumi * 1.e-5 * 3600. * 12.
        count = Lint * fB * Sigma_B * Br_kee * hlt_eff
        string = " & {:5.1f}".format(count)
        if ee_rate > l1_max : string = " &  -   "
        elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
        print(string,end="")
    # elif table == "counts_per_nfills" : 
        # Lint = peak_lumi * 1.e-5 * 3600. * 12.
        # count = Lint * fB * Sigma_B * Br_kee * hlt_eff * nfills
        # string = " & {:5.1f}".format(count)
        # if ee_rate > l1_max : string = " &  -   "
        # elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
        # print(string,end="")
    elif table == "counts_per_fb" : 
        Lint = _integrated_lumi # defined above (25.?)
        count = Lint * fB * Sigma_B * Br_kee * hlt_eff
        string = " & {:5.1f}".format(count)
        if ee_rate > l1_max : string = " &  -   "
        elif not l1_ok and hlt_ok : string = string.replace(" & "," & \gr{")+"}"
        print(string,end="")
    elif "menu" in table : 
        pass
    else :
        print("???")

################################################################################
# Print menu table

def print_menu(table,
               menu_dict,
               switch_lumi=_switch_lumi,
               switch_npu=_switch_npu,
               scaled_by_nb=None):

    ########## Print menu ##########
    if table == "menu" :
        default = [0.128, 0.405, 0.724, 0.968, 1.154, 1.383, 1.729] #@@???
        Lint = _integrated_lumi
        titles="Linst     NPU     L1 pT    HLT pT     "
        titles+="Spare     "
        titles+="L1 rate  HLT rate       AxE    L/fill    #/fill"
        if scaled_by_nb is not None: titles += "     Col"
        print(titles)
        for jj,(lumi,npu) in enumerate(zip(switch_lumi,switch_npu)) :
            peak_lumi = lumi[0]
            #if peak_lumi not in menu_dict :
            if npu not in menu_dict :
                print("{:4.2f}".format(peak_lumi)," & ",npu," \\\\")
                continue
            (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,capacity,count_per_fill_) = menu_dict[npu]
            count = Lint * fB * Sigma_B * Br_kee * hlt_eff
            Lint_per_fill = peak_lumi * 1.e-5 * 3600. * 12.
            count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff
            hlt_eff_per_rate = (hlt_eff*1.e4) / (ee_rate/1000.)
            orig_col = peak_lumi/scaled_by_nb if scaled_by_nb is not None else None
            print("".join(["{:4.2f}".format(peak_lumi),
                           "  &  {:4.1f}".format(npu),
                           "  &  {:5.1f}".format(l1_pt),
                           "  &  {:5.1f}".format(hlt_pt),
                           "  &  {:5.1f}".format(capacity/1000.),
                           "  & &  {:5.1f}".format(ee_rate/1000.),
                           "  &  {:5.2f}".format(hlt_rate/1000.),
                           "  &  {:5.2f}".format(hlt_eff*1.e4),
                           "  &  {:5.2f}".format(Lint_per_fill),
                           "  &  {:5.2f}".format(count_per_fill),
                           #"  &  {:6.1f}".format(count),
                           #"  &  {:5.3f}".format(hlt_eff_per_rate), # AxE/kHz
                           #"  &  {:5.2f}".format(hlt_eff*1.e4/default[jj]), # ???
                           "  &  {:3.1f}".format(orig_col) if orig_col is not None else "",
                           " \\\\ "]))

    ########## Print "hybrid" menu that makes use of prescales ##########
    elif table == "menu_prescaled" :
        #peak_lumi_ = 0.6 # this is the default trigger
        #(l1_pt_,hlt_pt_,ee_rate_,hlt_rate_,hlt_eff_,capacity_) = menu_dict[peak_lumi_]
        npu_ = 30 # this is the default trigger
        print(menu_dict[npu_])
        (l1_pt_,hlt_pt_,ee_rate_,hlt_rate_,hlt_eff_,capacity_,count_per_fill_) = menu_dict[npu_]

        Lint = _integrated_lumi # defined above
        print("Lint:", Lint)
        print("Linst   L1 pT   HLT pT    L1 rate HLT rate"
              "      AxE   L/fill   #/fill   #/Lint  AxE/kHz prescale")
        for jj,(lumi,npu) in enumerate(zip(switch_lumi,switch_npu)) :
            peak_lumi = lumi[0]
            #if peak_lumi not in menu_dict : 
            if npu not in menu_dict : 
                print(peak_lumi," \\\\") #@@ scaled_by_nb???
                continue
            #(l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,capacity) = menu_dict[peak_lumi]
            (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,capacity,count_per_fill) = menu_dict[npu]
            #@@prescale = ee_rate_ / ee_rate if ee_rate > 0. else 1.e9

            if npu > npu_ :
                (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,capacity,count_per_fill) = menu_dict[npu_]

#            histo_ee_rate = ee_file.Get('L1_DoubleEG'+
#                                        str(l1_pt).replace('.','p').replace('p0','')+
#                                        'er1p22_dR_' + str(dr_dict[l1_pt]).replace('.','p'))
#            ee_rate = histo_ee_rate.Eval(npu)*1000

            prescale = 1.
            if npu > npu_ :
                prescale = ee_rate_/capacity if capacity > 0. else 1.
                prescale = max(prescale,1.) # Ensure prescale >= 1.

            count = Lint * fB * Sigma_B * Br_kee * (hlt_eff/prescale)
            Lint_per_fill = peak_lumi * 1.e-5 * 3600. * 12. #@@ scaled_by_nb???
            count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * (hlt_eff/prescale)
            hlt_eff_per_rate = (hlt_eff*1.e4/prescale) / (ee_rate/1000./prescale)
            print(peak_lumi," & ", #@@ scaled_by_nb???
                  "{:5.1f}".format(l1_pt),
                  " & {:5.1f}".format(hlt_pt),
                  " & & {:5.1f}".format(ee_rate/1000./prescale),
                  " & {:5.2f}".format(hlt_rate/1000./prescale),
                  " & {:5.2f}".format(hlt_eff*1.e4/prescale),
                  " & {:5.2f}".format(Lint_per_fill),
                  " & {:5.2f}".format(count_per_fill),
                  " & {:5.1f}".format(count),
                  " & {:5.3f}".format(hlt_eff_per_rate),
                  " & {:5.2f}".format(prescale),
                  " \\\\ ")

#            histo_ee_rate = ee_file.Get('L1_DoubleEG'+
#                                        str(l1_pt_).replace('.','p').replace('p0','')+
#                                        'er1p22_dR_' + str(dr_dict[l1_pt_]).replace('.','p'))
#            ee_rate_ = histo_ee_rate.Eval(npu)*1000
#
#            prescale = ee_rate_/capacity if capacity > 0. else 1.
#            prescale = max(prescale,1.) # Ensure prescale >= 1.
#
#            count = Lint * fB * Sigma_B * Br_kee * (hlt_eff_/prescale)
#            Lint_per_fill = peak_lumi * 1.e-5 * 3600. * 12. #@@ scaled_by_nb???
#            count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * (hlt_eff_/prescale)
#            hlt_eff_per_rate = (hlt_eff_*1.e4/prescale) / (ee_rate_/1000./prescale)
#            print(peak_lumi," & ", #@@ scaled_by_nb???
#                  "{:5.1f}".format(l1_pt_),
#                  " & {:5.1f}".format(hlt_pt_),
#                  " & & {:5.1f}".format(ee_rate_/1000./prescale),
#                  " & {:5.2f}".format(hlt_rate_/1000./prescale),
#                  " & {:5.2f}".format(hlt_eff_*1.e4/prescale),
#                  " & {:5.2f}".format(Lint_per_fill),
#                  " & {:5.2f}".format(count_per_fill),
#                  " & {:5.1f}".format(count),
#                  " & {:5.3f}".format(hlt_eff_per_rate),
#                  " & {:5.2f}".format(prescale),
#                  " \\\\ ")
#    print()

################################################################################
# Loop through L1 pT and peak Linst scenarios and print tables (original method)

def print_tables_original(debug,
                          table,
                          l1_pts=_l1_pts,
                          l1_max=95000.,
                          allocation=0.,
                          hlt_max=-1.,
                          nbunches=2544,
                          nb_default=2544.):
    print("Table:",table)

    # Scale linearly according to nbunches (relative to 2544b)
    scaled_by_nb = (nbunches/nb_default)

    # Peak Linst values scaled according to # bunches (relative to 2544b from 2018)
    switch_lumi = [ ( float("{:.2f}".format(x1*scaled_by_nb)),
                      float("{:.2f}".format(x2*scaled_by_nb)) ) for (x1,x2) in _switch_cols ]

    # Store variables to build "menu"
    menu_dict = {}

    # Open HLT files, one per nPU value
    hlt_files = [ TFile(common_path+'ee/roc_hlt_pu'+str(npu)+'.root') for npu in _switch_npu ] 

    # Iterate through L1 pT thresholds
    for ii,l1_pt in enumerate(l1_pts) :

        # Get correct histogram for di-electron L1 trigger rate, parameterised vs nPU
        l1_name = \
            'L1_DoubleEG'+str(l1_pt).replace('.','p').replace('p0','')+\
            'er1p22_dR_' + str(dr_dict[l1_pt]).replace('.','p')
        ee_rate_his = _ee_file.Get(l1_name)
        
        if debug : print("L1 LOOP:",ii,l1_pt,ee_rate_his) # DEBUG
        if not debug and not "menu" in table : print("{:4.1f}".format(l1_pt),end="")

        # Iterate through L_inst values (and correponding nPU value and HLT file)
        for jj,(lumi,npu,hlt_file) in enumerate(reversed(list(zip(switch_lumi,_switch_npu,hlt_files)))) :

            # Peak Linst (first entry of 2-tupl from each entry in switch_lumi)
            peak_lumi = lumi[0] 

            # Evaluate and correct L1 rate for CMS menu based on peak L_inst
            l1_rate = _l1_his.Eval(peak_lumi)
            l1_rate -= _l1_rate_corr # remove any overspend @ 2E34 (i.e. "translate down")

            # Correct in case of a reduced number of bunches in machine
            l1_rate *= (nbunches/nb_default) #@@scaled_by_nb

            # Determine L1 spare capacity
            spare = l1_max - l1_rate # APPLY CMS-WIDE PRESCALE HERE (IF REQUIRED)!

            # Scale (i.e. "decay") di-electron rate allocation according to L_inst
            alloc = allocation*(l1_rate/(_l1_rate_2e34-_l1_rate_corr))
            
            # Determine rate of L1 di-electron trigger for given nPU
            ee_rate = ee_rate_his.Eval(npu)*1000

            # Correct in case of a reduced number of bunches in machine
            ee_rate *= scaled_by_nb

            # Check if the L1 di-ele rate satisfies the dedicate rate allocation or spare capacity
            l1_ok = ee_rate < (alloc+1.e-6) or ee_rate < (spare+alloc+1.e-6) # spare + epsilon

            #print("L1 pT:",l1_pt,"L! rate",ee_rate)#@@
            #if l1_ok : print(" L1 rate ok !!")#@@

            #if l1_ok == False : continue #@@

            # Get appropriate "graph" for HLT rate (dependant on L1 pT and L_inst)
            hlt_roc = hlt_file.Get('inv_pt' + str(l1_pt).replace('.','p'))
            hlt_n = hlt_roc.GetN()

            hlt_ok = False
            hlt_rate = -1
            hlt_eff = -1
            hlt_pt = -1

            # Get HLT pT threshold from "L1 pT <-> HLT pT" map
            hlt_pt = hlt_threshold_dict.get(l1_pt,4.0)
            index = _hlt_pts.index(hlt_pt) # Assumes 
            if index < hlt_n:
                indices = list(range(index,hlt_n))[0:1] # Currently: don't scan; only consider a single HLT pT threshold
                for idx in indices:
                    # print("test",idx,_hlt_pts[idx],hlt_roc.GetPointX(idx),hlt_max,hlt_roc.GetPointY(idx))
                    hlt_rate = hlt_roc.GetPointX(idx)
                    hlt_eff = hlt_roc.GetPointY(idx)
                    hlt_pt = _hlt_pts[idx]
                    # print("HLT pT:",hlt_pt,"HLT rate",hlt_rate,"eff:",hlt_eff)#@@
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
                    # Correct for reduced number of bunches in machine
                    hlt_rate *= scaled_by_nb
                    if hlt_max < 0. or hlt_rate < hlt_max:
                        hlt_ok = True
                        if debug : print("    HLT SCAN:",hlt_n,idx,hlt_pt,hlt_max,hlt_rate,hlt_eff)
                        break
            else: print("cannot find hltpt")
            if hlt_eff == -1: hlt_ok = False #print('!!!! This cannot happen !!!!')

#            print("TEST",
#                  "l1_pt",l1_pt,
#                  "peak_lumi",peak_lumi,
#                  "npu",npu,
#                  #"l1_total_raw {:.1f}".format(_l1_his_NEW.Eval(npu)*nbunches),
#                  #"file",subtract_ee[jj],
#                  #"ee_rate {:.1f}".format(subtract_ee_his_NEW.Eval(npu)*nbunches),
#                  #"l1_corrected {:.1f}".format(l1_rate_NEW),
#                  "l1_param {:.1f}".format(l1_rate),
#                  "spare {:.1f}".format(spare),
#                  "hlt_rate {:.1f}".format(hlt_rate),
#                  )

            # Store various metrics to later build menu
            if l1_ok and hlt_ok : 
                if npu in menu_dict:
                    Lint_per_fill = (peak_lumi*scaled_by_nb) * 1.e-5 * 3600. * 12.
                    count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff # NEW
                    count_per_fill_ = Lint_per_fill * fB * Sigma_B * Br_kee * menu_dict[npu][4] # STORED hlt_eff
                    if count_per_fill > count_per_fill_:
                     #if hlt_eff > menu_dict[npu][4]:
                        #print("OVERWRITING WITH HIGHER EFF!",peak_lumi,l1_pt,hlt_pt,hlt_eff,menu_dict[peak_lumi][4])
                        menu_dict[npu] = (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,spare+alloc,count_per_fill)
                else:
                    #print("ADDING PEAK LUMI AND EFF!",peak_lumi,l1_pt,hlt_pt,hlt_eff)
                    Lint_per_fill = (peak_lumi*scaled_by_nb) * 1.e-5 * 3600. * 12.
                    count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff # NEW
                    menu_dict[npu] = (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,spare+alloc,count_per_fill)

            # Print individual tables
            if not debug : print_metric_table_row(table,
                                                  l1_pt,
                                                  ee_rate,l1_rate,spare,hlt_pt,hlt_rate,hlt_eff,
                                                  l1_ok,hlt_ok,l1_max,
                                                  peak_lumi)
                    
            # Some debug statement
            if debug : print("  LUMI LOOP:",jj,peak_lumi,npu,\
               "L1:","{:.1f}".format(spare),"{:.1f}".format(ee_rate),\
               "HLT:",hlt_n,hlt_pt,"{:.1f}".format(hlt_rate),"{:.2f}".format(hlt_eff*1.e4),\
               hlt_ok)#,hlt_roc, # DEBUG

        if not debug and not "menu" in table : print("\\\\")

    # End of loop through L1 pT thresholds ...

    # Print menu table 
    print_menu(table,
               menu_dict,
               switch_lumi=switch_lumi, # local var (scaled_by_nb)
               switch_npu=_switch_npu)

################################################################################
# Loop through L1 pT and peak Linst scenarios and print tables, using Run 3 rates

def print_tables_run3rates(debug,
                           table,
                           l1_pts=_l1_pts,
                           l1_max=95000.,
                           allocation=0.,
                           hlt_max=-1.,
                           nbunches=2160,
                           nb_default=2160.):
    print("Table:",table)

    scaled_by_nb = (nbunches/nb_default)

    switch_cols = [(2.0, 1.7),(1.7, 1.5),(1.5, 1.3), (1.3, 1.1), (1.1, 0.9)]
    
    switch_npu = [ float("{:.1f}".format(x[0]*26.)) for x in switch_cols ]

    switch_lumi = [ ( float("{:.2f}".format(x1*scaled_by_nb)),
                      float("{:.2f}".format(x2*scaled_by_nb)) ) for (x1,x2) in switch_cols ]

    # Di-ele trigger rate to subtract from total rate vs peak Linst (high-to-low: 1.1 --> 2.0)
    keys = [ x1 for (x1,x2) in switch_cols ]
    entries = ["L1_DoubleEG10p5er1p22_dR_0p6",
               "L1_DoubleEG8p5er1p22_dR_0p7",
               #"L1_DoubleEG8er1p22_dR_0p7",# @ 1.6E34
               "L1_DoubleEG7er1p22_dR_0p8",
               "L1_DoubleEG6p5er1p22_dR_0p8",
               "L1_DoubleEG6er1p22_dR_0p8",]
    subtract_ee = dict(zip(keys,entries))

    # Store variables to build "menu"
    menu_dict = {}

    # Contains L1 di-electron rate (from data) vs nPU
    l1_file = TFile(common_path+'ee/L1_rate_total.root')

    # Contains L1 di-electron rate (from data) vs nPU
    ee_file = TFile(common_path+'ee/L1_rate_DoubleEG.root')

    # Open HLT files for different nPU values
    hlt_files = [ TFile(common_path+'ee/roc_hlt_pu'+str(npu)+'.root') for npu in _switch_npu ] 

    # Iterate through L1 pT thresholds
    for ii,l1_pt in enumerate(l1_pts) :

        # Get correct histogram for di-electron L1 trigger rate, parameterised vs nPU
        l1_name = \
            'L1_DoubleEG'+str(l1_pt).replace('.','p').replace('p0','')+\
            'er1p22_dR_' + str(dr_dict[l1_pt]).replace('.','p')
        ee_rate_his = ee_file.Get(l1_name)

        if debug : print("L1 LOOP:",ii,l1_pt,ee_rate_his) # DEBUG
        if not debug and not "menu" in table : print("{:4.1f}".format(l1_pt),end="")

        # Iterate through L_inst values (and correponding nPU value and HLT file)
        for jj,(column,lumi,npu,hlt_file) in enumerate(reversed(list(zip(switch_cols,switch_lumi,switch_npu,hlt_files)))) :

            # Peak Linst (first entry of 2-tupl from each entry in switch_lumi)
            peak_lumi = lumi[0]
            
            # Determine rate of L1 di-electron trigger for given nPU
            ee_rate = ee_rate_his.Eval(npu) * nbunches # units are already Hz

            # L1 total rate
            name = 'L1_rate_{:3.1f}'.format(column[0]).replace('.','p')
            l1_fit = l1_file.Get(name)
            l1_rate = l1_fit.Eval(npu)

            # Di-ele contribution to L1 total rate
            subtract_ee_his = ee_file.Get(subtract_ee.get(column[0]))
            subtract_ee_rate = subtract_ee_his.Eval(npu)
            
            # L1 total rate corrected for di-ele contribution
            l1_rate -= subtract_ee_rate
            l1_rate *= nbunches

            # Dedicated allocation
            alloc = allocation

            # Determine L1 spare capacity
            spare = l1_max - l1_rate

            print(", ".join(["L1 pT {:4.1f}".format(l1_pt),
                             "Col {:3.1f}".format(column[0]),
                             "Peak {:4.2f}".format(peak_lumi),
                             "NPU {:4.1f}".format(npu),
                             "L1tot (raw) {:8.1f}".format(l1_fit.Eval(npu)*nbunches),
                             "Seed {:28s} {:8.1f}".format(subtract_ee.get(column[0]),
                                                          subtract_ee_his.Eval(npu)*nbunches),
                             "L1tot (corr) {:8.1f}".format(l1_rate),
                             "spare {:8.1f}".format(spare)]),
                             "ee rate {:8.1f}".format(ee_rate))

            # Check if the L1 di-ele rate satisfies the dedicate rate allocation or spare capacity
            l1_ok = ee_rate < (alloc+1.e-6) or ee_rate < (spare+alloc+1.e-6) # spare + epsilon

            #print("L1 pT:",l1_pt,"L! rate",ee_rate)#@@
            #if l1_ok : print(" L1 rate ok !!")#@@

            #if l1_ok == False : continue #@@

            # Get appropriate "graph" for HLT rate (dependant on L1 pT and L_inst)
            hlt_roc = hlt_file.Get('inv_pt' + str(l1_pt).replace('.','p'))
            hlt_n = hlt_roc.GetN()

            hlt_ok = False
            hlt_rate = -1
            hlt_eff = -1
            hlt_pt = -1

            if options.limit:

                # Extract HLT pT, rate, and efficiency if limited to 300 Hz (obsolete?)
                #hlt_max = max_bw_hlt[peak_lumi]
                #hlt_ok = False
                #hlt_rate = -1
                #hlt_eff = -1
                #hlt_pt = -1
                hlt_pt = hlt_threshold_dict.get(l1_pt,4.0) # Get HLT pT threshold from "L1 pT <-> HLT pT" map 
                hltpt_list = np.arange(4, 11, 0.5).tolist()
                index = hltpt_list.index(hlt_pt)
                if index < hlt_n:
                    for kk in range(index,index+1): #hlt_n):# range(hlt_n): #@@ CONSIDER ONLY ONE HLT PT THRESHOLD PER L1 THRESHOLD
                        ip = kk#hlt_n - kk - 1
                        #print("test",kk,ip,_hlt_pts[ip],hlt_roc.GetPointX(ip),hlt_max,hlt_roc.GetPointY(ip))
                        hlt_rate = hlt_roc.GetPointX(ip)
                        hlt_eff = hlt_roc.GetPointY(ip)
                        hlt_pt = _hlt_pts[ip]
                        #print("HLT pT:",hlt_pt,"HLT rate",hlt_rate,"eff:",hlt_eff)#@@
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
                        # Correct for reduced number of bunches in machine
                        hlt_rate *= scaled_by_nb
                        if hlt_max < 0. or hlt_rate < hlt_max:
                            #print(" HLT rate ok!!",)#@@
                            hlt_ok = True
                            if debug : print("    HLT SCAN:",hlt_n,ip,hlt_pt,hlt_max,hlt_rate,hlt_eff)
                            break
                else: print("cannot find hltpt")
                if hlt_eff == -1: hlt_ok = False #print('!!!! This cannot happen !!!!')
                    
            else:

                # Extract HLT pT, rate, and efficiency
                #hlt_ok = True
                #hlt_rate = -1
                #hlt_eff = -1
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
                                hltptstr in corrs_dict[npu][l1ptstr].keys():
                            corr = corrs_dict[npu][l1ptstr][hltptstr]
                            if corr is not None: hlt_rate *= corr
                            else: print("null value",npu,l1ptstr,hltptstr)
                        else: print("unknown thresholds",npu,l1ptstr,hltptstr)
                    # Correct for reduced number of bunches in machine
                    hlt_rate *= scaled_by_nb
                    hlt_ok = True
                else: print("cannot find hltpt")

            # Store various metrics to later build menu
            if l1_ok and hlt_ok : 
                if npu in menu_dict:
                    Lint_per_fill = (peak_lumi*scaled_by_nb) * 1.e-5 * 3600. * 12.
                    count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff # NEW
                    count_per_fill_ = Lint_per_fill * fB * Sigma_B * Br_kee * menu_dict[npu][4] # STORED hlt_eff
                    if count_per_fill > count_per_fill_:
                     #if hlt_eff > menu_dict[npu][4]:
                        #print("OVERWRITING WITH HIGHER EFF!",peak_lumi,l1_pt,hlt_pt,hlt_eff,menu_dict[peak_lumi][4])
                        menu_dict[npu] = (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,spare+alloc,count_per_fill)
                else:
                    #print("ADDING PEAK LUMI AND EFF!",peak_lumi,l1_pt,hlt_pt,hlt_eff)
                    Lint_per_fill = (peak_lumi*scaled_by_nb) * 1.e-5 * 3600. * 12.
                    count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff # NEW
                    menu_dict[npu] = (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,spare+alloc,count_per_fill)

            # Print individual tables
            if not debug : print_metric_table_row(table,
                                                  l1_pt,
                                                  ee_rate,l1_rate,spare,hlt_pt,hlt_rate,hlt_eff,
                                                  l1_ok,hlt_ok,l1_max,
                                                  peak_lumi)
                    
            # Some debug statement
            if debug : print("  LUMI LOOP:",jj,peak_lumi,npu,\
               "L1:","{:.1f}".format(spare),"{:.1f}".format(ee_rate),\
               "HLT:",hlt_n,hlt_pt,"{:.1f}".format(hlt_rate),"{:.2f}".format(hlt_eff*1.e4),\
               hlt_ok)#,hlt_roc, # DEBUG

        if not debug and not "menu" in table : print("\\\\")

    # End of loop through L1 pT thresholds ...

    # Print menu table 
    print_menu(table,
               menu_dict,
               switch_lumi=switch_lumi,
               switch_npu=switch_npu)

################################################################################
# Loop through L1 pT and peak Linst scenarios and print tables, using Run 3 rates
#
def print_tables_run3rates_param(debug,
                                 table,
                                 l1_pts=_l1_pts,
                                 l1_max=95000.,
                                 allocation=5000.,
                                 hlt_max=-1.,
                                 nbunches=_nb_default,
                                 nb_default=_nb_default):
    print("Table:",table)
    
    # Scale by # bunches
    scaled_by_nb = (nbunches/nb_default)

#    # Prescale columns identified by "peak Linst" values
#    switch_cols = [ ( float("{:.1f}".format(x/10.)),
#                      float("{:.1f}".format((x-1)/10.)) ) for x in range(25,5,-1) ]
#    
#    # Uncomment to override with columns used for original menu
#    #switch_cols = _switch_cols
#
#    # Corresponding NPU used to extract L1 total and di-ele rates and HLT rates
#    # Fit: nPU = 32.0[e-34] * Linst + 0.9, established from data with 2160b
#    # Assuming gradient of 33, scaled to 2748b, 33*(2160/2748) = 26
#    # Hence, at nominal operating parameters, 2E34 @ 2748b, NPU = 52 ...
#    # ... which is in line with expectations
#    switch_npu = [ float("{:.1f}".format(x[0]*26.)) for x in switch_cols ]
#
#    # Corresponding peak Linst after scaling for # bunches in machine
#    # Fit: Linst = nPU/26 (from above)
#    # Then peak Linst scales linearly with number of bunches relative to nominal (2748b)
#    switch_lumi = [ ( float("{:.2f}".format(x1*scaled_by_nb)),
#                      float("{:.2f}".format(x2*scaled_by_nb)) ) for (x1,x2) in switch_cols ]
#
#    print("HERE")
#    print([x1 for (x1,x2) in switch_cols])
#    print(switch_npu)
#    print([x1 for (x1,x2) in switch_lumi])

    # Required prescale columns identified by "peak Linst" values
    switch_cols = [ ( float("{:.1f}".format(x/10.)),
                      float("{:.1f}".format((x-1)/10.)) ) for x in range(22,5,-1) ]

#    switch_cols = [ ( float("{:.1f}".format(x/10.)),
#                      float("{:.1f}".format((x-1)/10.)) ) for x in [20,17,15,13,11] ]

    # Corresponding peak Linst after *inverse* scaling for # bunches in machine
    # Fit: Linst = nPU/26 (from above)
    # Then peak Linst scales linearly with number of bunches relative to nominal (2748b)
    switch_lumi = [ ( float("{:.2f}".format(x1/scaled_by_nb)),
                      float("{:.2f}".format(x2/scaled_by_nb)) ) for (x1,x2) in switch_cols ]

    # Corresponding NPU used to extract L1 total and di-ele rates and HLT rates
    # Fit: nPU = 32.0[e-34] * Linst + 0.9, established from data with 2160b
    # Assuming gradient of 33, scaled to 2748b, 33*(2160/2748) = 26
    # Hence, at nominal operating parameters, 2E34 @ 2748b, NPU = 52 ...
    # ... which is in line with expectations
    switch_npu = [ float("{:.1f}".format(x[0]*26.)) for x in switch_lumi ]

    #print(switch_cols)
    #print(switch_npu)
    #print(switch_lumi)

#    #print("HERE")
#    import bisect
#    lumis = list(reversed([x1 for (x1,x2) in switch_cols]))
#    upper = bisect.bisect_left(lumis, 2.2)
#    lower = bisect.bisect_left(lumis, 0.6)
#    first = len(lumis) - (upper+1)
#    last  = len(lumis) - lower
#    #print(lumis)
#    #print(len(lumis),lower,upper,first,last)
#    #print(lumis[lower:upper])
#    #print("HERE")
#    switch_cols = switch_cols[first:last]]
#    switch_npu = switch_npu[first:last]
#    switch_lumi = switch_lumi[first:last]]
#    print([x1 for (x1,x2) in switch_cols[first:last]])
#    print(switch_npu[first:last])
#    print([x1 for (x1,x2) in switch_lumi[first:last]])

    # Contains parameterised L1 di-electron rate (from data) vs nPU
    ee_file = TFile(common_path+'ee/L1_rate_DoubleEG.root')

    # Contains ratios of L1 di-electron rates (2022 vs 2018) vs nPU (determined @ 2544b)
    ee_file_corr = TFile(common_path+'ee/L1_rate_DoubleEG_corr.root')

    # Contains parameterised HLT di-electron rate (from data) vs nPU
    hlt_file = TFile(common_path+'ee/HLT_rate_DoubleEle.root')

    # Store variables to build "menu"
    menu_dict = {}

    # Iterate through L1 pT thresholds
    for ii,l1_pt in enumerate(l1_pts) :

        # Get correct histogram for di-electron L1 trigger rate, parameterised vs nPU
        l1_name = \
            'L1_DoubleEG'+str(l1_pt).replace('.','p').replace('p0','')+\
            'er1p22_dR_' + str(dr_dict[l1_pt]).replace('.','p')
        ee_rate_his = ee_file.Get(l1_name)
        ee_rate_his_corr = ee_file_corr.Get(l1_name)

        if debug : print("L1 LOOP:",ii,l1_pt,ee_rate_his) # DEBUG
        if not debug and not "menu" in table : print("{:4.1f}".format(l1_pt),end="")

        # Iterate through L_inst values (and correponding nPU value and HLT file)
        for jj,(lumi,npu) in enumerate(reversed(list(zip(switch_cols,switch_npu)))) :

            peak_lumi = lumi[0]
            
            # Determine rate of L1 di-electron trigger for given nPU
            ee_rate = ee_rate_his.Eval(npu) * nbunches

            # Determine ratios of L1 di-electron rates (2022 vs 2018) vs nPU (determined @ 2544b)
            ee_rate_corr = ee_rate_his_corr.Eval(npu)

            # PARAMETERISED L1 total rate vs NPU (with di-electron subtracted)
            l1_rate_per_pu = 100000./52.

            # Estimated L1 rate
            l1_rate = npu * l1_rate_per_pu * scaled_by_nb
            spare = l1_max - l1_rate
            alloc = allocation 

#            print("TEST",
#                  "l1_pt",l1_pt,
#                  "peak_lumi",peak_lumi,
#                  "npu",npu,
#                  "l1_total {:.1f}".format(l1_rate),
#                  "ee_rate {:.1f}".format(ee_rate),
#                  "alloc {:.1f}".format(alloc),
#                  "spare {:.1f}".format(spare))

            # Check if the L1 di-ele rate satisfies the dedicate rate allocation or spare capacity
            l1_ok = ee_rate < (alloc+1.e-6) or ee_rate < (spare+alloc+1.e-6) # spare + epsilon

            # Don't allow use of L1 thresholds below 11 GeV for peak Linst > 2E34
            if peak_lumi >  2.0 and l1_pt < 11.0 : l1_ok = False
            #if _l1_max > 99999. and peak_lumi >  2.0 and l1_pt < 11.0 : l1_ok = False
            #if _l1_max > 94999. and peak_lumi >= 2.0 and l1_pt < 11.0 : l1_ok = False

            # Get appropriate "graph" for HLT rate (dependent on L1 pT)
            hlt_roc = hlt_file.Get(l1_name)
            hlt_n = hlt_roc.GetN()
            hlt_title = hlt_roc.GetTitle()

            hlt_ok = False
            hlt_rate = -1
            hlt_eff = -1
            hlt_pt = -1

            hlt_eff = float(hlt_title.split("eff=")[1])*1.e-4
            # Rate, scale by: # bunches; ratio of L1 rates (2022 c.f. 2018); (arbitrary) 20%
            hlt_rate = hlt_roc.Eval(npu) * scaled_by_nb * 1.5 # * ee_rate_corr * 1.2
            hlt_pt = hlt_threshold_dict.get(l1_pt,4.0)
            if hlt_max < 0. or hlt_rate < hlt_max: hlt_ok = True

#            print(", ".join(["L1 pT {:4.1f}".format(l1_pt),
#                             "Peak {:4.2f}".format(peak_lumi),
#                             "NPU {:4.1f}".format(npu),
#                             "L1 tot {:8.1f}".format(l1_rate),
#                             "spare {:8.1f}".format(spare),
#                             "ee rate {:8.1f}".format(ee_rate),
#                             "hlt_rate {:6.1f}".format(hlt_rate),
#                             "hlt_eff {:4.2f}".format(hlt_eff*1.e4),
#                             "l1 ok {:1.0f}".format(l1_ok),
#                             "hlt ok {:1.0f}".format(hlt_ok)]))

            # Store various metrics to later build menu
            if l1_ok and hlt_ok : 
                if npu in menu_dict:
                    Lint_per_fill = (peak_lumi*scaled_by_nb) * 1.e-5 * 3600. * 12.
                    count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff # NEW
                    count_per_fill_ = Lint_per_fill * fB * Sigma_B * Br_kee * menu_dict[npu][4] # STORED hlt_eff
                    if count_per_fill > count_per_fill_:
                     #if hlt_eff > menu_dict[npu][4]:
                        #print("OVERWRITING WITH HIGHER EFF!",peak_lumi,l1_pt,hlt_pt,hlt_eff,menu_dict[peak_lumi][4])
                        menu_dict[npu] = (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,spare+alloc,count_per_fill)
                else:
                    #print("ADDING PEAK LUMI AND EFF!",peak_lumi,l1_pt,hlt_pt,hlt_eff)
                    Lint_per_fill = (peak_lumi*scaled_by_nb) * 1.e-5 * 3600. * 12.
                    count_per_fill = Lint_per_fill * fB * Sigma_B * Br_kee * hlt_eff # NEW
                    menu_dict[npu] = (l1_pt,hlt_pt,ee_rate,hlt_rate,hlt_eff,spare+alloc,count_per_fill)
                    
            # Some debug statement
            if debug : print("  LUMI LOOP:",jj,peak_lumi,npu,\
               "L1:","{:.1f}".format(spare),"{:.1f}".format(ee_rate),\
               "HLT:",hlt_n,hlt_pt,"{:.1f}".format(hlt_rate),"{:.2f}".format(hlt_eff*1.e4),\
               hlt_ok)#,hlt_roc, # DEBUG

        if not debug and not "menu" in table : print("\\\\")

    # End of loop through L1 pT thresholds ...

    # Print menu table 
    print_menu(table,
               menu_dict,
               switch_lumi=switch_cols,
               switch_npu=switch_npu)
    #scaled_by_nb=scaled_by_nb)

################################################################################
# Configurations ...
def configurations(action):

    # Configure ...
    if action == "print_original_tables":
        # Print individual metrics
        for table in ["ee_rate","l1_rate","spare","l1_pt","hlt_pt","hlt_rate","eff",
                      "eff_per_rate","lint","counts_per_fill","counts_per_fb",]:
            print_tables_original(debug=False,
                                  table=table,
                                  l1_max=_l1_max_2018,
                                  allocation=0.)
        # Print "menu" for different allocations
        for allocation in [0.,5000.,10000.,20000.]:
            print("allocation:",allocation)
            print_tables_original(debug=False,
                                  table="menu",
                                  l1_pts=[4.5,5.5,6.0,6.5,7.0,8.0,10.5], # constrain to the six thresholds of the default menu
                                  l1_max=_l1_max_2018,
                                  hlt_max=-1.,
                                  allocation=allocation)

    elif action == "print_prescaled_menu":
        # Print original menu and hybrid (i.e. prescaled) menu
        allocation = 5000.
        print("THIS NEEDS SOME WORK ...")
        print("allocation:",allocation)
        print_tables_original(debug=False,
                              table="menu",
                              l1_max=_l1_max_2018,
                              hlt_max=-1.,
                              allocation=allocation)
        print("THE PRESCALED TABLE BELOW NEEDS SOME WORK ...")
        print_tables_original(debug=False,
                              table="menu_prescaled",
                              l1_max=_l1_max_2018,
                              hlt_max=-1.,
                              allocation=allocation)
        
    elif action == "print_original_menu_vs_nbunches":
        # Print "menu" for different nbunches 
        # NOTA BENE: rerun with different values of hlt_max!!!
        for nbunches in [600,900,1200,1500,1800,2100,2544]:
            print("nbunches:",nbunches)
            print_tables_original(debug=False,
                                  table="menu",
                                  l1_pts=[4.5,5.5,6.0,6.5,7.0,8.0,10.5], # constrain to the six thresholds of the default menu
                                  l1_max=_l1_max,
                                  allocation=5000.,
                                  hlt_max=2000.,
                                  nbunches=nbunches,
                                  nb_default=2544.)

    elif action == "print_tables_run3rates":
        nbunches = 2740
        print("nbunches:",nbunches)
        print_tables_run3rates(debug=False,
                               table="menu",
                               l1_pts=[4.5,5.5,6.0,6.5,7.0,8.0,10.5],
                               l1_max=_l1_max,
                               allocation=5000.,
                               hlt_max=-1.,
                               nbunches=nbunches,
                               nb_default=_nb_default)

    elif action == "print_tables_run3rates_param":
        #for nbunches in [2740, 2640, 2544, 2400, 2160, 1922]:
        for nbunches in [2740]:
            print("nbunches:",nbunches)
            print_tables_run3rates_param(debug=False,
                                         table="menu",
                                         #l1_pts=[4.5,5.5,6.0,6.5,7.0,8.0,10.5], #@@ temp!!
                                         #l1_pts=[6.0,6.5,7.0,8.0,8.5,10.5], #@@ temp!!
                                         #l1_pts=[6.],
                                         l1_max=_l1_max,
                                         allocation=5000.,
                                         hlt_max=-1.,
                                         nbunches=nbunches,
                                         nb_default=_nb_default)

    else:
        pass

################################################################################
# Print various tables here
if __name__ == "__main__":

    # Choose here ...
    action = [
        "print_original_tables",           # Prints metric tables and original menus
        "print_prescaled_menu",            # Prints original menu and "hybrid" (i.e. prescaled) menu
        "print_original_menu_vs_nbunches", # Prints original menu vs number of bunches
        "print_tables_run3rates",          # Prints metric tables and menus using Run 3 observed L1 trigger rates
        "print_tables_run3rates_param",    # Prints metric tables and menus using parameterised Run 3 observed L1 trigger rates
        ][2]

    print("########## Configuration ##########") 
    print("Action:",action) 
    print() 
    configurations(action)
