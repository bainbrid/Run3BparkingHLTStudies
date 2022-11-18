from bisect import insort, bisect_left, bisect_right, bisect
from collections import deque,Counter
from datetime import datetime
from itertools import islice
import numpy as np
from numpy import median
import os, copy
from ROOT import TDatime, TFile, TGraph

################################################################################
# Define comon path to input files (i.e. either EOS via lxplus or local)

lxplus=True
common_path='./root/current/'
if lxplus==True: common_path='/eos/cms/store/group/phys_bphys/DiElectronX/RootFiles4Run3Parking/current/'

################################################################################

# BBbar inclusive cross section 
# from http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
Sigma_B = 4.6940e+11 # fb (femtobarn!!) 

# Fragmentation fraction for B+-
fB = 0.4

# Branching fraction for "at least one rare B->Kee decay" per event (simplified to 2 * 4.4E-7)
Br_kee = 2*4.5e-7

################################################################################
# Dictionary that maps pT threshold to deltaR requirement at Level-1
dr_dict = {
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

################################################################################
# Maximum HLT bandwidth in 2018 vs Linst, from Sara's presentation:
# https://indico.cern.ch/event/1032638/contributions/4336416/
max_bw_hlt = {
    2.2:1515,
    2.0:1515,
    1.7:1740,
    1.5:1929,
    1.3:2163.5,
    1.1:2463,
    0.9:2929,
    0.6:3791, # Default TSG column?
    0.7:3791,0.47:3791,0.24:3791 # NEED TO UPDATE FOR LOWER LINST???
}

################################################################################
# Pairwise (L1,HLT) thresholds to avoid "vertical regions" in ROCs
l1_threshold_list = np.arange(4, 11.5, 0.5).tolist()
hlt_threshold_list = [4.0,4.0,4.0,4.0,4.0, # L1: 4.0->6.0
                      4.5,5.0,5.0,5.0,5.5, # L1: 6.5->8.5
                      6.0,6.5,6.5,6.5,6.5  # L1: 9.0->11.0
                      ]
hlt_threshold_dict = dict(zip(l1_threshold_list,hlt_threshold_list))

################################################################################
# List of PU values ...
# ... that map to Linst values: 2.0, 1.7, 1.5, 1.3, 1.1, 0.9, 0.6E34
npu_list = [56, 48, 42, 36, 30, 25, 17]

################################################################################
# Parse .csv file to extract example "luminosity profile" for 2018

def extractLumiProfiles(original=False,max_duration=12*3600) :

    # Golden JSON
    #https://cmsoms.cern.ch/cms/runs/report?cms_run=324980&cms_run_sequence=GLOBAL-RUN
    #"324980": [[39, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]],
    golden = [[53, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]]

    # Parse csv file 
    times = []
    lumis = []
    for line in open('LumiData/LumiData_2018_20200401.csv', 'r'):
        if line.find('324980:7321')==-1: continue

        line = line.rstrip().split(',')
        ls = line[1].split(':')[0]

        flag = False
        for lrange in golden:
            if int(ls) >= min(lrange) and int(ls) <= max(lrange): flag = True
        if not flag: continue

        time = line[2].split(' ')
        time = TDatime(int(time[0].split('/')[2])+2000, 
                       int(time[0].split('/')[0]), 
                       int(time[0].split('/')[1]), 
                       int(time[1].split(':')[0]), 
                       int(time[1].split(':')[1]), 
                       int(time[1].split(':')[2]))
        times.append(time.Convert())

        Linst = float(line[5])*0.0001
        lumis.append(Linst)
    
    # Start at zero
    min_time = min(times)
    times = [number - min_time for number in times]

    # Check if times are sorted
    if(times != sorted(times)):
        print("Times not sorted!")
        quit()

    # Return originals, before smoothing or truncating
    if original: return times,lumis

    # Smooth with running median
    window = 11 # has to be odd
    lumis = RunningMedian(lumis,window) # shortens by window-1
    for i in range((window-1)/2): # pad
        lumis.insert(0,lumis[0])
        lumis.insert(-1,lumis[-1])

    # Truncate to 12 hours    
    times,lumis = zip(*filter(lambda time: 
                              time[0]<max_duration, 
                              zip(times,lumis)))

    return times,lumis

################################################################################
# Creates various luminosity profiles
# Can be based on a real but modified (e.g. smoothed) lumi profile (from above) ...
# ... or a "synthetic" profile based on exponential parameterisation of real profile

def createLumiProfiles(output='root/lumiprof.root',backup=False,synthetic=True) :

    if os.path.exists(output): 
        print("Warning! File already exists!")
        if backup : 
            print("Renaming...")
            today = datetime.today().strftime('%Y%m%d_%H%M%S')
            os.rename(output,output.replace(".root","_{:s}.root".format(today)))

    # Configure duration (fill, lumi-levelling)
    max_duration = 12*3600
    level_duration = 6*3600

    # Produce output file
    file = TFile(output, 'recreate')

    # Falling from 1.8E34 (original, not truncated!)
    times_orig,lumis_orig = extractLumiProfiles(original=True)
    graph = TGraph()
    graph.SetName('falling_from_1p8e34_original')
    graph.SetTitle('falling_from_1p8e34_original')
    idx = 0
    for time, lumi in zip(times_orig, lumis_orig):
        graph.SetPoint(idx, time, lumi)
        idx += 1
    graph.Write()

    times,lumis = None,None
    if synthetic:
        times = np.arange(0.,max_duration,23.) # 12-hr fill
        lumis = [ np.exp(0.6-(x*2.3e-5)) for x in np.arange(0.,48.*3600.,23.) ] # "48-hr fill"
        peak_lumi = lumis[0]
    else:
        times,lumis = extractLumiProfiles(original=False)
        peak_lumi = median(lumis[:10]) # Just to protect against outliers
    
    # Falling from 1.8E34 
    graph = TGraph()
    graph.SetName('falling_from_1p8e34')
    graph.SetTitle('falling_from_1p8e34')
    idx = 0
    for time, lumi in zip(times, lumis):
        graph.SetPoint(idx, time, lumi) 
        idx += 1
    graph.Write()

    # Falling from 0.9E34 
    graph = TGraph()
    graph.SetName('falling_from_0p9e34')
    graph.SetTitle('falling_from_0p9e34')
    idx = 0
    for time, lumi in zip(times, lumis):
        graph.SetPoint(idx, time, lumi/2.) #@@ halved!
        idx += 1
    graph.Write()

    # Determine median value for time difference
    diff = [ y-x for x, y in zip(times[0::], times[1::]) ]
    step = median(diff)
    lumi_range = np.arange(0, level_duration, step).tolist() 

    # Levelled at 2.0E34 
# ORIG
#    levelling = zip([2.0,1.7,1.5,1.3,1.1,0.9,0.6,0.45,0.30,0.15],
#                    ["2p0e34","1p7e34","1p5e34","1p3e34","1p1e34",
#                     "0p9e34","0p6e34","4p5e33","3p0e33","1p5e33"])
# FILIP
#    levelling = zip([2.1,1.8,1.4,0.94,
#                     0.7,0.47,0.24],
#                    ["2p1e34","1p8e34","1p4e34","0p9e34","0p7e34","0p4e34","0p2e34"])
# HYBRID
    levelling = zip([2.0,1.7,1.5,1.3,1.1,0.9,
                     0.7,0.47,0.24],
                    ["2p0e34","1p7e34","1p5e34","1p3e34","1p1e34","0p9e34",
                     "0p7e34","0p4e34","0p2e34"])
    for level,name in levelling:
        graph = TGraph()
        graph.SetName('levelled_at_'+name)
        graph.SetTitle('levelled_at_'+name)
        idx = 0
        for time in lumi_range:
            graph.SetPoint(idx, time, level) # set levelled values
            idx+=1
        times_adj = [number + level_duration for number in times] # start at end of levelling
        if synthetic:
            lumi_adj = max(0., level - peak_lumi) # only adjust upwards
            lumis_adj = [lumi + lumi_adj for lumi in lumis]
            index = len(lumis_adj) - bisect_left(lumis_adj[::-1], level) # find lumi
            if index is not None: lumis_adj = lumis_adj[index:] # start from there
        else:
            lumi_adj = level - peak_lumi
            lumis_adj = [lumi + lumi_adj for lumi in lumis] # adjust lumi to levelled
        for time, lumi in zip(times_adj, lumis_adj):
            if time > max_duration : continue
            graph.SetPoint(idx, time, max(0.,lumi))
            idx += 1
        graph.Write()
    
    # Write to file and close
    file.Write()
    file.Close()

################################################################################
# Utility methods

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def hackRate(rate,which_lumi):
    #linst=[0.6,0.45,0.30,0.15]
    linst=[0.7,0.47,0.24,0.06]
    idx=None
    try: idx = linst.index(which_lumi)
    except ValueError: return rate
    if idx is not None and idx>0: return rate * linst[idx]/linst[0]
    else : return rate

def hackRateNPU(rate,npu):
    #linst=[0.6,0.45,0.30,0.15]
    npus=[17,25,30,36,42,48,56]
    if npu < npus[0] : return rate * (npu/npus[0])
    else : return rate
#    idx=None
#    try: idx = npu.index(npus)
#    except ValueError: return rate
#    if idx is not None and idx>0: return rate * npus[idx]/npus[0]
#    else : return rate

def scaleGraph(graph,scale) :
    for i in range(graph.GetN()) :
        graph.SetPointY(i,graph.GetPointY(i)*scale)
    return graph

def RunningMedian(seq, M):

    # Running median, used to smooth lumi profiles
    # Taken from https://code.activestate.com/recipes/578480-running-median-mean-and-mode/

    seq = iter(seq)
    s = []   
    m = M // 2

    # Set up list s (to be sorted) and load deque with first window of seq
    s = [item for item in islice(seq,M)]    
    d = deque(s)

    # Simple lambda function to handle even/odd window sizes    
    median = lambda : s[m] if bool(M&1) else (s[m-1]+s[m])*0.5

    # Sort it in increasing order and extract the median ("center" of the sorted window)
    s.sort()    
    medians = [median()]   

    # Now slide the window by one point to the right for each new position (each pass through 
    # the loop). Stop when the item in the right end of the deque contains the last item in seq
    for item in seq:
        old = d.popleft()          # pop oldest from left
        d.append(item)             # push newest in from right
        del s[bisect_left(s, old)] # locate insertion point and then remove old 
        insort(s, item)            # insert newest such that new sort is not required        
        medians.append(median())  
    return medians
