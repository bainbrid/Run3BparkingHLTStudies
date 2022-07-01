# Local environment:
# conda activate root.6.26.0
# python metrics.py

################################################################################
# Imports
from __future__ import print_function
from ROOT import TDatime, TGraph, TFile, TH1F, TCanvas, TLegend, gROOT, gStyle, TH2F, TGaxis, gPad, TLine, TPaveText, TLatex, kBlue, kRed, kBlack
from DisplayManager import DisplayManager, add_Preliminary, add_Private, add_CMS, add_lumi, applyLegendSettings
from officialStyle import officialStyle
from ctypes import c_double
import copy, os, sys
import numpy as np
import json

################################################################################
# Common definitions 
from common import ensureDir, common_path, dr_dict, hackRate, \
    hlt_threshold_dict, fB, Sigma_B, Br_kee, createLumiProfiles, scaleGraph

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
                  help="limit HLT rate to 300 Hz")
parser.add_option('-p', '--prescaled', # prescale-adjusted approach ...
                  action="store_true",
                  default=False,
                  dest='prescaled',
                  help="Adjust prescales rather than pT thresholds at L1")
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

# Normalise to integrated luminosity
integrated_lumi = 86.4 # 12-hour fill @ 2E34 delivers 0.864/fb 

# L1 pT thresholds
l1_ptrange = np.arange(4,11,0.5).tolist() 

# List of L_inst values to consider
# Linst = nPU*0.0357338 - 0.0011904
switch_lumi = [(2.2, 2.0),
               (2.0, 1.7), (1.7, 1.5), (1.5, 1.3), (1.3, 1.1), (1.1, 0.9), (0.9, 0.7),
               (0.7, 0.47), (0.47, 0.24), (0.24, 0.06), (0.06, 0.)]

# List of nPU values to consider
# nPU = (Linst+0.0011904)/0.0357388
switch_npu = [62, # NEED TO UPDATE FOR HIGHEST LINST?
              56, 48, 42, 36, 30, 25,
              17, 17, 17, 17] # NEED TO UPDATE FOR LOWER LINST

# Maximum L1 trigger bandwidth
l1_max = 95000.

# Rate (in Hz??) dedicated to new di-muon seeds in Run 3
#dimuon = 0. # Redundant?

# Rate (in kHz !!) dedicated to single-muon seed in Run 3
single_muon = 5.

# Determine correction to L1 rate for total CMS menu
# i.e. ensure rate does not exceed l1_max @ 2E34 (i.e. remove any overspend)
idx_2e34 = 1
l1_rate_2e34 = l1_his.Eval(switch_lumi[idx_2e34][0])
l1_rate_corr = max(0.,l1_rate_2e34-l1_max) # correct overspend only

idx_1p7e34 = 2
l1_rate_1p7e34 = l1_his.Eval(switch_lumi[idx_1p7e34][0])
l1_rate_1p7e34 -= l1_rate_corr

# List of peak L_inst (upper bounds from switch_lumi, list of 2-tuples ...)
peak_lumis = [ h for (h,l) in switch_lumi ]

# L1 prescales for entire CMS menu (options 1 and 2)
l1_rate_prescale = [1.94,
                    1.94, 1.57, 1.32, 1.09, 1., 1.,
                    1., 1., 1., 1.]
#l1_rate_prescale = [5.]*len(peak_lumis) # Option 1
l1_rate_prescale = [1.94,                # Option 2
                    1.94, 1.57, 1.32, 1.09, 1., 1.,
                    1., 1., 1., 1.]
l1_rate_prescale_map = dict(zip(peak_lumis,l1_rate_prescale)) # Not applied by default - check below !!

# This is the desired L1 pT threshold above which triggers are prescaled
pt_max = 6.0

# List of tuples recorded for later analysis
summary = []

# Create luminosity profiles ...
output='root/lumiprof.root'
createLumiProfiles(output,synthetic=True)
file_profile = TFile(output)

# Get the relevant luminosity profiles
profiles = [
# FALLING
    #file_profile.Get('falling_from_1p8e34_original'),
    #file_profile.Get('falling_from_1p8e34'),
    #file_profile.Get('falling_from_0p9e34'),
# LUMI-LEVELLED
    file_profile.Get('levelled_at_2p0e34'),
    file_profile.Get('levelled_at_1p7e34'),
    file_profile.Get('levelled_at_1p5e34'),
    file_profile.Get('levelled_at_1p3e34'),
    file_profile.Get('levelled_at_1p1e34'),
    file_profile.Get('levelled_at_0p9e34'),
    file_profile.Get('levelled_at_0p7e34'),
    file_profile.Get('levelled_at_0p4e34'),
    file_profile.Get('levelled_at_0p2e34'),
]

# Corrections to HLT rates - obsolete?
corrs_dict = {}
if options.corrected:
    for npu in [56, 48, 42, 36, 30, 25, 17]:
        filename = common_path+'rates/corrections_' + str(npu) + '.json'
        infile = open(filename,'r')
        dct = json.load(infile)
        corrs_dict[npu] = dct
        infile.close()

################################################################################
# Process the files ...

# Iterate through the different luminosity profiles
for idx,profile in enumerate(profiles):

    name = profile.GetName()
    print()
    print("name:",name)

    # Graph for spare capacity
    gspare = TGraph()
    gspare.SetName('spare')
    gspare.SetTitle('L1 spare')
    gspare.SetMarkerSize(0)
    gspare.SetMarkerColor(kRed+1)
    gspare.SetLineStyle(1)
    gspare.SetLineColor(kRed+1)
    gspare.SetLineWidth(3)

    # List of graphs (one graph per dedicated rate allocation)
    graphs = []

    # Dedicated L1 rate allocation
    allocations = [0,5000,10000,20000][1:2]

    # Iterate through different rate allocations
    for style,allocation in enumerate(allocations):

        # Graph for given rate allocation
        graph = TGraph()
        graph.SetName('name'+'_allocation_' + str(allocation))
        graph.SetTitle(str(int(allocation/100.)/10.))
        graph.SetMarkerSize(0)
        graph.SetMarkerColor(kBlue+2)
        graph.SetLineStyle(1)
        graph.SetLineColor(kBlue+len(allocations)-style)
        graph.SetLineWidth(style+1)

        # Init for graph showing spare capacity
        if style==0: gspare.SetPoint(0, 0., 0.)

        # List to record various variables related to estimates
        switches = []

        # Counters
        ls_cntr = 0     # lumi section
        total_lumi = 0  # L_int
        total_count = 0 # Estimated candidates

        # Loop through "luminosity profile" graph
        old_lumi = -1
        time_elapsed = 0.
        for ii in range(profile.GetN()):
            ls_cntr += 1

            # Get L_inst from profile
            Linst = max(profile.GetPointY(ii), 1.e-9)

            # Time interval (since last PointX, i.e. 1 LS = 23 secs) and total elapsed
            time_interval = profile.GetPointX(ii) - profile.GetPointX(ii-1)
            time_elapsed += time_interval

            # Count L_int
            total_lumi += Linst * 1.e-5 * time_interval # L_inst units: 1E-5 Hz/pb
            
            # Check which luminosity it is ... 
            which_lumi = -1
            which_npu = -1
            for index, sl in enumerate(switch_lumi):
                if sl[1] < Linst and Linst <= sl[0]: # Is L_inst within upper and lower bounds
                    which_lumi = sl[0]            # If yes, store L_inst
                    which_npu = switch_npu[index] # If yes, store nPU
                    break
            if which_lumi==-1: 
                print('   WARNING!!: no corresponding lumis!!!',"Linst",Linst)
                continue

            # Check if the "L_inst column" has changed 
            switch = True if ( ( which_lumi != old_lumi ) or 
                               ( ii == profile.GetN()-1 ) ) else False
            old_lumi = which_lumi

            # Evaluate and correct L1 rate for CMS menu based on peak L_inst
            l1_rate = l1_his.Eval(which_lumi)
            l1_rate -= l1_rate_corr # remove any overspend @ 2E34 (i.e. "translate down")
            l1_rate = hackRate(l1_rate,which_lumi) # <<< HACK HACK HACK

            # APPLY PRESCALE FOR L1 DI-ELE MENU IF USING PRESCALE-ADJUSTED APPROACH!!
            #l1_rate /= l1_rate_prescale_map.get(which_lumi,1.) 

            # Determine L1 spare capacity 
            spare = l1_max - l1_rate # APPLY CMS-WIDE PRESCALE HERE (IF REQUIRED)!

            # HACK TO ACCOUNT FOR SINGLE MUON RATE (NORMALISED TO 5 KHZ)
            mu_rate = single_muon * 1.e3 * l1_rate/l1_rate_1p7e34 if which_lumi < idx_1p7e34 else 0.

            # Scale (i.e. "decay") di-electron rate allocation according to L_inst
            alloc = allocation*(l1_rate/(l1_rate_2e34-l1_rate_corr)) # allocation scaled down by lumi

            # Add spare capacity to graph (only for first rate allocation in list)
            if style==0: gspare.SetPoint(ls_cntr, profile.GetPointX(ii), spare)

            # Determine rate of L1 di-electron trigger for given nPU
            which_l1pt = l1_ptrange[-1] 
            flag_park = False
            l1_ee_rate = None
            prescale = 1.
            for pt in l1_ptrange:
                histo_ee_rate = ee_file.Get('L1_DoubleEG'+
                                            str(pt).replace('.','p').replace('p0','')+
                                            'er1p22_dR_' + str(dr_dict[pt]).replace('.','p'))
                ee_rate = histo_ee_rate.Eval(which_npu)*1000
                ee_rate = hackRate(ee_rate,which_lumi) # <<< HACK HACK HACK

                # Check if the L1 di-ele rate satisfies the dedicate rate allocation or spare capacity
                l1_ok = ee_rate < (alloc+1.e-6) or ee_rate < (alloc+spare-mu_rate+1.e-6) # spare + epsilon
                if l1_ok:
                    which_l1pt = pt
                    flag_park = True
                    l1_ee_rate = ee_rate
                    break
            if not flag_park:
                #print('\t L1 not available')
                continue

            # If using prescale-adjusted method, OVERRIDE above settings ...
            if options.prescaled and which_l1pt > pt_max : # ... only if above pt_max

                # Old settings
                which_l1pt_ = which_l1pt
                flag_park_ = flag_park
                l1_ee_rate_ = l1_ee_rate
                prescale_ = prescale

                # New settings (reset)
                which_l1pt = l1_ptrange[-1] 
                flag_park = False
                l1_ee_rate = None
                prescale = 1.

                # Determine rate of L1 di-electron trigger for given nPU
                histo_ee_rate = ee_file.Get('L1_DoubleEG'+
                                            str(pt_max).replace('.','p').replace('p0','')+
                                            'er1p22_dR_' + str(dr_dict[pt_max]).replace('.','p'))
                ee_rate = histo_ee_rate.Eval(which_npu)*1000
                ee_rate = hackRate(ee_rate,which_lumi) # <<< HACK HACK HACK

                # Determine prescale to apply
                l1_rate_max = max(alloc,alloc+spare-mu_rate)-1.e-6 # max allowed? (minus epsilon)
                prescale = ee_rate/l1_rate_max if l1_rate_max > 0. else 1.
                prescale = max(prescale,1.) # Ensure prescale >= 1.
                ee_rate = ee_rate / prescale # Apply prescale

                # Check if the L1 di-ele rate satisfies the dedicate rate allocation or spare capacity
                l1_ok = ee_rate < (alloc+1.e-6) or ee_rate < (alloc+spare-mu_rate+1.e-6) # spare + epsilon
                if l1_ok:
                    which_l1pt = pt_max
                    flag_park = True
                    l1_ee_rate = ee_rate
                if not flag_park: 
                    #print('\t L1 not available')
                    continue

            # Determine HLT threshold ...
            hlt_file = TFile(common_path+'ee/roc_hlt_pu' + str(which_npu) + '.root')
            hlt_roc = hlt_file.Get('inv_pt' + str(which_l1pt).replace('.','p'))
            hlt_n = hlt_roc.GetN()

            if options.limit:

                # Extract corrections to rates from JSON

                hlt_max = 300. #max_bw_hlt[peak_lumi] # <-- Limit to 300 Hz for prompt reco
                hlt_ok = True
                hlt_rate = -1
                hlt_eff = -1
                which_hltpt = -1
                for kk in range(hlt_n):
                    ip = kk#hlt_n - kk - 1
                    # print("test",kk,ip,hlt_pts[ip],hlt_roc.GetPointX(ip),hlt_max,hlt_roc.GetPointY(ip))
                    hlt_rate = hlt_roc.GetPointX(ip) / prescale # Apply prescale
                    hlt_eff = hlt_roc.GetPointY(ip) / prescale # Apply prescale
                    which_hltpt = np.arange(4, 11, 0.5).tolist()[ip]
                    if options.corrected:
                        l1ptstr = str(which_l1pt)
                        hltptstr = str(which_hltpt)
                        if which_npu in corrs_dict.keys() and \
                                l1ptstr in corrs_dict[which_npu].keys() and \
                                hltptstr in corrs_dict[which_npu][l1ptstr].keys():
                            corr = corrs_dict[which_npu][l1ptstr][hltptstr]
                            if corr is not None: hlt_rate *= corr
                            else: print("null value",which_npu,l1ptstr,hltptstr)
                        else: print("unknown thresholds",which_npu,l1ptstr,hltptstr)
                    if hlt_rate < hlt_max:
                        break
                if hlt_eff == -1: hlt_ok = False #print('!!!! This cannot happen !!!!')

            else:

                # Extract HLT pT, rate, and efficiency

                hlt_rate = -1
                hlt_eff = -1
                which_hltpt = hlt_threshold_dict.get(which_l1pt,4.0)
                hltpt_list = np.arange(4, 11, 0.5).tolist()
                index = hltpt_list.index(which_hltpt)
                if index < hlt_roc.GetN():
                    hlt_rate = hlt_roc.GetPointX(index) / prescale # Apply prescale
                    hlt_eff = hlt_roc.GetPointY(index) / prescale # Apply prescale
                    if options.corrected:
                        l1ptstr = str(which_l1pt)
                        hltptstr = str(which_hltpt)
                        if which_npu in corrs_dict.keys() and \
                                l1ptstr in corrs_dict[which_npu].keys() and \
                                hltptstr in  corrs_dict[which_npu][l1ptstr].keys():
                            corr = corrs_dict[which_npu][l1ptstr][hltptstr]
                            if corr is not None: hlt_rate *= corr
                            else: print("null value",which_npu,l1ptstr,hltptstr)
                        else: print("unknown thresholds",which_npu,l1ptstr,hltptstr)
                else: print("cannot find hltpt")

            # Determine candidates per time interval and increment counter
            count = Linst * 1.e-5 * fB * Sigma_B * Br_kee * hlt_eff * time_interval
            total_count += count
            graph.SetPoint(ls_cntr, profile.GetPointX(ii), total_count)
            
            # Store various metrics for later analysis
            if switch:
                switches.append((time_elapsed,
                                 which_l1pt,
                                 l1_ee_rate,
                                 spare,
                                 which_lumi,
                                 total_lumi,
                                 total_count))
                print(
                    " ".join(["ii:",str("{:5.0f}".format(ii)),
                              "time:",str("{:5.0f}".format(time_elapsed)),
                              "Peak:",str("{:4.2f}".format(which_lumi)),
                              "Linst:",str("{:.3f}".format(Linst)),
                              "Lint:",str("{:.3f}".format(total_lumi)),
                              "alloc:",str("{:5.0f}".format(alloc)),
                              "spare:",str("{:5.0f}".format(spare)),
                              "pT:",str("{:4.1f}".format(which_l1pt)),
                              "L1 rate:",str("{:5.0f}".format(l1_ee_rate)),
                              "mu rate:",str("{:5.0f}".format(mu_rate)),
                              "HLT rate:",str("{:5.0f}".format(hlt_rate)),
                              "HLT eff:",str("{:.6f}".format(hlt_eff)),
                              "N:",str("{:6.3f}".format(total_count)),
                              str("prescale: {:5.1f}".format(prescale) if options.prescaled else ""),
                          ]))
 
        graphs.append(copy.deepcopy(graph))

        # Record peak lumis and L1 pTs
        peak_lumis = [ peak_lumi for _,_,_,_,peak_lumi,_,_ in switches ][:-1]
        l1pts = [ l1pt for _,l1pt,_,_,_,_,_ in switches ][:-1]

        # Record time spent at each ...
        start_t = [ time_elapsed for time_elapsed,_,_,_,_,_,_ in switches ]
        end_t = start_t[1:]
        diff_t = [ e-s for e,s in zip(end_t,start_t) ] 

        # Record lumis spent at each ...
        start_L = [ total_lumi for _,_,_,_,_,total_lumi,_ in switches ]
        end_L = start_L[1:]
        diff_L = [ e-s for e,s in zip(end_L,start_L) ] 

        # Record counts accumulated at each ...
        start_C = [ total_count for _,_,_,_,_,_,total_count in switches ]
        end_C = start_C[1:]
        diff_C = [ e-s for e,s in zip(end_C,start_C) ] 

        # Dict of peak_lumis:[diff_t,diff_L,diff_C]
        peaks = {}
        for pL,dt,dL,dC in zip(peak_lumis,diff_t,diff_L,diff_C) :
            if pL not in peaks.keys() : peaks[pL] = [0.,0.,0.]
            peaks[pL][0] += dt
            peaks[pL][1] += dL
            peaks[pL][2] += dC

        # Dict of pt_thresholds:[diff_t,diff_L,diff_C]
        thresholds = {}
        for pt,dt,dL,dC in zip(l1pts,diff_t,diff_L,diff_C) :
            if pt not in thresholds.keys() : thresholds[pt] = [0.,0.,0.]
            thresholds[pt][0] += dt
            thresholds[pt][1] += dL
            thresholds[pt][2] += dC

        print("SUMMARY:") 
        print('Lumi profile name:       ', name) 
        print('L1 total rate:           ', l1_max)
        print('Di-ele allocation:       ', allocation)
        total_time = profile.GetPointX(profile.GetN()-1)
        print('Time duration            ', total_time, 's')
        print('Total lumi:              ', total_lumi,"fb^-1")
        print('Total # or Kee events:   ', total_count)
        factor = integrated_lumi/(total_lumi)
        print('Target lumi:             ', integrated_lumi)
        print('Luminosity factor:       ', factor)
        print('Expected # of Kee events:', total_count*factor)

        # Currently append for each allocation
        summary.append((name,allocation,total_lumi,total_count,peaks,thresholds))

    # First, lumi profile (only), then add estimates
    for only_profile in [True,False]:

        # Margins
        gStyle.SetPadTopMargin(0.08)
        gStyle.SetPadLeftMargin(0.14)
        if only_profile: gStyle.SetPadRightMargin(0.03)
        else: gStyle.SetPadRightMargin(0.26) # Extra space for multiple axes

        # Create canvas
        canvas = TCanvas('canvas_' +  name)
        canvas.SetGridx()
        canvas.SetGridy()
        gPad.SetTicks(1,0) # No tick marks on RHS

        # Create frame for lumi profile
        ymax = 2.5 # Maximum for y-axis (Linst)
        frame = TH2F('frame_'+name, 'frame_'+name, 
                     profile.GetN(), 0, total_time, # time axis
                     100, 0, ymax ) # counts axis
        frame.GetXaxis().SetTitle('Time [s]')
        frame.GetYaxis().SetTitle('L_{inst} [10^{34} Hz/cm^{2}]')
        frame.SetTitleOffset(1.2)
        frame.Draw()

        # Draw lumi profile
        profile.SetLineColor(kBlack)
        profile.SetLineWidth(3)
        profile.SetMarkerSize(0)
        profile.Draw('same')
        canvas.Update()

        # Create legend
        t = TLatex()
        t.SetTextAlign(11)
        t.SetTextSize(0.035)
        if not only_profile: t.DrawLatexNDC(0.47,0.87,"Allocation [kHz]")
        upper = 0.86
        lower = upper-len(graphs)*0.04
        leg = TLegend(0.45,lower,0.7,upper)
        applyLegendSettings(leg)
        leg.SetTextSize(0.035)
                
        # Draw graphs and right-hand axes
        if not only_profile:

            # Spare capacity
            right_max = l1_max * ymax / 2.0
            scale = gPad.GetUymax()/right_max
            scaleGraph(gspare,scale)
            gspare.Draw("same")
            xshift = (gPad.GetUxmax()-gPad.GetUxmin())*0.25
            axis2 = TGaxis(gPad.GetUxmax()+xshift,
                           gPad.GetUymin(),
                           gPad.GetUxmax()+xshift,
                           gPad.GetUymax(),
                           0,right_max,510,"+L")
            axis2.SetLineColor(kRed+1)
            axis2.SetLabelColor(kRed+1)
            axis2.SetTitleOffset(1.4)
            axis2.SetTitle('Spare L1 trigger capacity [Hz]')
            axis2.SetTitleColor(kRed+1)
            axis2.Draw()

            # Cumu candidates
            right_max = 12.5
            scale = gPad.GetUymax()/right_max
            for ii, graph in enumerate(graphs[::-1]): # in reverse order
                scaleGraph(graph,scale)
                graph.Draw('lsame')
                leg.AddEntry(graph, graph.GetTitle(), 'l')
            axis1 = TGaxis(gPad.GetUxmax(),
                           gPad.GetUymin(),
                           gPad.GetUxmax(),
                           gPad.GetUymax(),
                           0,right_max,510,"+L")
            axis1.SetLineColor(kBlue+2)
            axis1.SetLabelColor(kBlue+2)
            axis1.SetTitleOffset(1.2)
            axis1.SetTitle('Cumulative number of Kee candidates')
            axis1.SetTitleColor(kBlue+2)
            axis1.Draw()

        # Draw legend
        if not only_profile: leg.Draw()

        # Add CMS labels
        l2=add_Private(text="#it{Private Work}",x=0.27)
        l2.Draw("same")
        l3=add_CMS(x=0.14)
        l3.Draw("same")
        l4=add_lumi(total_lumi,x=0.83 if only_profile else 0.6)
        l4.Draw("same")

        # Save canvas
        canvas.Update()
        filename = 'plots/'+name+'.pdf'
        if not only_profile: filename = filename.replace("plots/","plots/estimates_for_")
        canvas.SaveAs(filename)
        del canvas

# Print summary
for name,allocation,lumi,count,peaks,thresholds in summary:
    print('Name:',"{:30s}".format(name),
          'L1 alloc [Hz]:',"{:6.0f}".format(allocation),
          'Lint [/fb] per fill:',"{:6.3f}".format(lumi),
          'Counts, per fill:',"{:5.1f}".format(count),
          ', per /fb:',"{:5.1f}".format(count/lumi),
          ', per 25/fb:',"{:5.1f}".format(integrated_lumi*count/lumi)
    )

# Tables for estimates

dct = {}
for name,allocation,lumi,count,peaks,thresholds in summary:
    if name not in dct.keys() : dct[name] = {}
    dct[name][allocation] = (lumi,count)

print("Per fill (Lint/fill, Counts for each allocation):")
for name,counts in reversed(dct.items()):
    print(name," ",end="")
    print("{:0.2f}".format(counts[allocations[0]][0])," ",end="")
    for alloc in allocations: print("{:5.1f}".format(counts[alloc][1])," ",end="")
    print()

nfills = {}
#nfills['falling_from_1p8e34'] = 1
nfills['levelled_at_0p2e34'] = 3
nfills['levelled_at_0p4e34'] = 3
nfills['levelled_at_0p7e34'] = 3
nfills['levelled_at_0p9e34'] = 12
nfills['levelled_at_1p1e34'] = 12
nfills['levelled_at_1p3e34'] = 12
nfills['levelled_at_1p5e34'] = 12
nfills['levelled_at_1p7e34'] = 12
nfills['levelled_at_2p0e34'] = 12

print("All fills (Lint/fill, Nfills, Lint all fills, Counts for each allocation):")
Lint_total = 0.
total = [0.]*len(allocations)
for name,counts in reversed(dct.items()):
    print(name," ",end="")
    Lint = counts[allocations[0]][0]
    nfill = nfills.get(name,1.)
    Lint_total += Lint*nfill
    print("{:4.2f}".format(Lint)," ",end="")
    print("{:2.0f}".format(nfill)," ",end="")
    print("{:5.2f}".format(Lint*nfill)," ",end="")
    for ialloc,alloc in enumerate(allocations): 
        count = counts[alloc][1]
        tot = count*nfill
        print("{:6.1f}".format(tot)," ",end="")
        total[ialloc] += tot
    print()
print(" "*29,"{:5.2f} ".format(Lint_total),"  ".join(["{:6.1f}".format(t) for t in total]))

dct2 = {}
total_dt = 0.
total_dL = 0.
total_dC = 0.
for name,_,_,_,peaks,_ in reversed(summary):
    for peak_lumi,[dt,dL,dC] in peaks.items():
        if peak_lumi not in dct2.keys() : dct2[peak_lumi] = [0.,0.,0.]
        tot_dt = dt*nfills.get(name,1.)
        dct2[peak_lumi][0] += tot_dt
        total_dt += tot_dt
        tot_dL = dL*nfills.get(name,1.)
        dct2[peak_lumi][1] += tot_dL
        total_dL += tot_dL
        tot_dC = dC*nfills.get(name,1.)
        dct2[peak_lumi][2] += tot_dC
        total_dC += tot_dC
print("Peak Linst, Duration [s], Lint [/fb], Counts, (%)") 
for peak,[dt,dL,dC] in dct2.items():
    print("{:4.1f} {:8.0f} {:4.2f} {:6.2f} ({:4.2f}) {:5.1f} ({:4.2f})".format(peak,
                                                                               dt,dt/total_dt,
                                                                               dL,dL/total_dL,
                                                                               dC,dC/total_dC))
print("Tot: {:8.0f}      {:6.2f}        {:5.1f}".format(total_dt,total_dL,total_dC))

dct3 = {}
total_dt = 0.
total_dL = 0.
total_dC = 0.
for name,_,_,_,_,thresholds in reversed(summary):
    for threshold,[dt,dL,dC] in thresholds.items():
        if threshold not in dct3.keys() : dct3[threshold] = [0.,0.,0.]
        tot_dt = dt*nfills.get(name,1.)
        dct3[threshold][0] += tot_dt
        total_dt += tot_dt
        tot_dL = dL*nfills.get(name,1.)
        dct3[threshold][1] += tot_dL
        total_dL += tot_dL
        tot_dC = dC*nfills.get(name,1.)
        dct3[threshold][2] += tot_dC
        total_dC += tot_dC
print("pT threshold, Duration [s], Lint [/fb], Counts, (%)") 
for threshold,[dt,dL,dC] in dct3.items():
    print("{:4.1f} {:8.0f} {:4.2f} {:6.2f} ({:4.2f}) {:5.1f} ({:4.2f})".format(threshold,
                                                                               dt,dt/total_dt,
                                                                               dL,dL/total_dL,
                                                                               dC,dC/total_dC))
print("Tot: {:8.0f}      {:6.2f}        {:5.1f}".format(total_dt,total_dL,total_dC))

print("DOES IT WORK FOR MULTIPLE ALLOCATIONS!!!")
