#import os, math, sys
#from ROOT import TFile, TH1F, gROOT, TTree, Double, TChain, TLorentzVector, TVector3
#import numpy as num

from TreeProducer_eff import *
from DeltaR import deltaR, deltaPhi
import copy
import random
import numpy as np
import itertools
from ROOT import TLorentzVector
import argparse

def hlt_criteria(chain, idx):
    
    flag = False
    
#    if chain.eg_et[idx] > float(pt) \
    if chain.eg_pms2[idx] < 10000 \
       and chain.eg_invEInvP[idx] < 0.2 \
       and chain.eg_trkDEtaSeed[idx] < 0.01 \
       and chain.eg_trkDPhi[idx] < 0.2 \
       and chain.eg_trkChi2[idx] < 50 \
       and chain.eg_trkValidHits[idx] >= 10 \
       and chain.eg_trkNrLayerIT[idx] >= 2:

        flag = True


    return flag




def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

me = 0.000511

from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)


parser = argparse.ArgumentParser(description='example e/gamma HLT analyser')
parser.add_argument('inputs',nargs="+",default="input.root",help='input filename')
parser.add_argument('--output','-o',default="output.root",help='output filename')
parser.add_argument('--type','-t',default="mc",help='output filename')
args = parser.parse_args()


#parser.add_option("-o", "--out", default='ratetest.root', type="string", help="output filename", dest="out")

#parser.add_option("-f", "--file", default='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Trigger/job/data_4gev_22/Myroot_simple.root', type="string", help="file", dest="file")
#parser.add_option("-f", "--file", default='test.root', type="string", help="file", dest="file")


print('filename=', args.inputs)


#(options, args) = parser.parse_args()

#print(options)


out = TreeProducer_eff(args.output)

chain = ROOT.TChain('egHLTRun3Tree', 'tree')
#chain.AddFile(options.file)

for _file in args.inputs:
    chain.AddFile(_file)

Nevt = chain.GetEntries()




#from optparse import OptionParser, OptionValueError
#usage = "usage: python runTauDisplay_BsTauTau.py"
#parser = OptionParser(usage)
#
#parser.add_option("-o", "--out", default='efftest.root', type="string", help="output filename", dest="out")
#parser.add_argument('inputs',nargs="+",help='input filename')
##parser.add_option("-f", "--file", default='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Trigger/job/mc_widedPhiMaxLowEtGrad/Myroot.root', type="string", help="file", dest="file")
##parser.add_option("-f", "--file", default='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Trigger/job/Winter21_official_mc/Myroot.root', type="string", help="file", dest="file")
#
#(options, args) = parser.parse_args()
#
#print(options)
#
#
#name = options.out
#
#
#out = TreeProducer_eff(options.out)
#
#chain = ROOT.TChain('egHLTRun3Tree', 'tree')
#chain.AddFile(options.file)
#
#Nevt = chain.GetEntries()

print('Total Number of events = ', Nevt)
evtid = 0


for evt in range(Nevt):
    chain.GetEntry(evt)

    if evt%10000==0: print('{0:.2f}'.format(float(evt)/float(Nevt)*100.), '% processed')

    
    if(len(chain.eg_gen_pt)!=2): 
        #print('# of gen =', len(chain.eg_gen_pt), 'detected!!! continue ...')
        continue

    # ordering ... 

    idx1 = 0
    idx2 = 1
    
    if chain.eg_gen_pt[1] > chain.eg_gen_pt[0]:
        idx1 = 1
        idx2 = 0

    epairs = []


    # check over L1 objects and see if it matches ... 

    match_index_l1 = -1
    match_index_hlt = -1
    
    for gindex, idx in enumerate([idx1,idx2]):

        getattr(out, 'gen_e' + str(gindex+1) + '_pt')[0] = chain.eg_gen_pt[idx]
        getattr(out, 'gen_e' + str(gindex+1) + '_eta')[0] = chain.eg_gen_eta[idx]
        getattr(out, 'gen_e' + str(gindex+1) + '_phi')[0] = chain.eg_gen_phi[idx]
        getattr(out, 'gen_e' + str(gindex+1) + '_vz')[0] = chain.eg_gen_vz[idx]


        tlv = TLorentzVector()
        tlv.SetPtEtaPhiM(chain.eg_gen_pt[idx], chain.eg_gen_eta[idx], chain.eg_gen_phi[idx], me)
        epairs.append(copy.deepcopy(tlv))
        

        l1index = -1
        maxl1_dr = 99
        maxl1_deta = 99
        maxl1_dphi = 99

        for jj in range(len(chain.eg_l1eg_et)):

            if abs(chain.eg_l1eg_eta[jj]) > 1.218: continue

            dr = deltaR(chain.eg_gen_eta[idx], chain.eg_gen_phi[idx],
                        chain.eg_l1eg_eta[jj], chain.eg_l1eg_phi[jj])

            if maxl1_dr > dr and jj!=match_index_l1:

                maxl1_dr = dr
                maxl1_dphi = deltaPhi(chain.eg_gen_phi[idx], chain.eg_l1eg_phi[jj])
                maxl1_deta = chain.eg_gen_eta[idx] - chain.eg_l1eg_eta[jj]

                l1index = jj

                if gindex==0: match_index_l1 = jj


        getattr(out, 'gen_e' + str(gindex+1) + '_l1_dr')[0] = maxl1_dr
        getattr(out, 'gen_e' + str(gindex+1) + '_l1_deta')[0] = maxl1_deta
        getattr(out, 'gen_e' + str(gindex+1) + '_l1_dphi')[0] = maxl1_dphi

        if maxl1_dr!=99:
            getattr(out, 'e' + str(gindex+1) + '_l1_pt')[0] = chain.eg_l1eg_et[l1index]
            getattr(out, 'e' + str(gindex+1) + '_l1_eta')[0] = chain.eg_l1eg_eta[l1index]
            getattr(out, 'e' + str(gindex+1) + '_l1_phi')[0] = chain.eg_l1eg_phi[l1index]
        else:
            getattr(out, 'e' + str(gindex+1) + '_l1_pt')[0] = -1
            getattr(out, 'e' + str(gindex+1) + '_l1_eta')[0] = -1
            getattr(out, 'e' + str(gindex+1) + '_l1_phi')[0] = -1


        hltindex = -1
        maxhlt_dr = 99
        maxhlt_dphi = 99
        maxhlt_deta = 99

        for jj in range(len(chain.eg_et)):

            if abs(chain.eg_eta[jj]) > 1.2: continue
#            if not hlt_criteria(chain, jj): continue

            dr = deltaR(chain.eg_gen_eta[idx], chain.eg_gen_phi[idx],
                        chain.eg_eta[jj], chain.eg_phi[jj]) 
            
           
            if maxhlt_dr > dr and jj!=match_index_hlt:
#            if maxhlt_dr > dr:
                
                maxhlt_dr = dr
                maxhlt_dphi = deltaPhi(chain.eg_gen_phi[idx], chain.eg_phi[jj])
                maxhlt_deta = chain.eg_gen_eta[idx] - chain.eg_eta[jj]

                hltindex = jj

                if gindex==0: match_index_hlt = jj


        getattr(out, 'gen_e' + str(gindex+1) + '_hlt_dr')[0] = maxhlt_dr
        getattr(out, 'gen_e' + str(gindex+1) + '_hlt_dphi')[0] = maxhlt_dphi
        getattr(out, 'gen_e' + str(gindex+1) + '_hlt_deta')[0] = maxhlt_deta

        if maxhlt_dr!=99:
            getattr(out, 'e' + str(gindex+1) + '_hlt_pt')[0] = chain.eg_et[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_eta')[0] = chain.eg_eta[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_phi')[0] = chain.eg_phi[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_energy')[0]  = chain.eg_energy[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_rawEnergy')[0]  = chain.eg_rawEnergy[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_phiWidth')[0]  = chain.eg_phiWidth[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_nrClus')[0]  = chain.eg_nrClus[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_seedId')[0]  = chain.eg_seedId[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_seedDet')[0]  = chain.eg_seedDet[hltindex]
            #getattr(out, 'e' + str(gindex+1) + '_hlt_seednrCrystals')[0]  = chain.eg_seednrCrystals[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_sigmaIEtaIEta')[0]  = chain.eg_sigmaIEtaIEta[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_sigmaIEtaIEtaNoise')[0]  = chain.eg_sigmaIEtaIEtaNoise[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_ecalPFIsol')[0]  = chain.eg_ecalPFIsol[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_hcalPFIsol')[0]  = chain.eg_hcalPFIsol[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkIsol')[0]  = chain.eg_trkIsol[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkChi2')[0]  = chain.eg_trkChi2[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkMissHits')[0]  = chain.eg_trkMissHits[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkNrLayerIT')[0]  = chain.eg_trkNrLayerIT[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkValidHits')[0]  = chain.eg_trkValidHits[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_hcalHForHoverE')[0]  = chain.eg_hcalHForHoverE[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_invESeedInvP')[0]  = chain.eg_invESeedInvP[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_invEInvP')[0]  = chain.eg_invEInvP[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkDEta')[0]  = chain.eg_trkDEta[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkDEtaSeed')[0]  = chain.eg_trkDEtaSeed[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkDPhi')[0]  = chain.eg_trkDPhi[hltindex]
            getattr(out, 'e' + str(gindex+1) + '_hlt_pms2')[0]  = chain.eg_pms2[hltindex]
            #getattr(out, 'e' + str(gindex+1) + '_hlt_dxy')[0]  = chain.eg_gsfdxy[hltindex]
            #getattr(out, 'e' + str(gindex+1) + '_hlt_vz')[0]  = chain.eg_track_vz[hltindex]


        else:
            getattr(out, 'e' + str(gindex+1) + '_hlt_pt')[0] = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_eta')[0] = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_phi')[0] = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_energy')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_rawEnergy')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_phiWidth')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_nrClus')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_seedId')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_seedDet')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_seednrCrystals')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_sigmaIEtaIEta')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_sigmaIEtaIEtaNoise')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_ecalPFIsol')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_hcalPFIsol')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkIsol')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkChi2')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkMissHits')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkNrLayerIT')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkValidHits')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_hcalHForHoverE')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_invESeedInvP')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_invEInvP')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkDEta')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkDEtaSeed')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_trkDPhi')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_pms2')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_dxy')[0]  = -1
            getattr(out, 'e' + str(gindex+1) + '_hlt_vz')[0]  = -1



    if getattr(out, 'gen_e1_l1_dr')[0]!=99 and getattr(out, 'gen_e2_l1_dr')[0]!=99:
        dr = deltaR(getattr(out, 'e1_l1_eta')[0], getattr(out, 'e1_l1_phi')[0], 
                    getattr(out, 'e2_l1_eta')[0], getattr(out, 'e2_l1_phi')[0])

        tlv1 = TLorentzVector()
        tlv2 = TLorentzVector()
        
        tlv1.SetPtEtaPhiM(getattr(out, 'e1_l1_pt'), getattr(out, 'e1_l1_eta'), getattr(out, 'e1_l1_phi'), me)
        tlv2.SetPtEtaPhiM(getattr(out, 'e2_l1_pt'), getattr(out, 'e2_l1_eta'), getattr(out, 'e2_l1_phi'), me)


        getattr(out, 'l1_eedr')[0] = dr 
        getattr(out, 'l1_mee')[0] = (tlv1 + tlv2).M()

    else:
        getattr(out, 'l1_eedr')[0] = -1
        getattr(out, 'l1_mee')[0] = -1


    if getattr(out, 'gen_e1_hlt_dr')[0]!=99 and getattr(out, 'gen_e2_hlt_dr')[0]!=99:
        dr = deltaR(getattr(out, 'e1_hlt_eta')[0], getattr(out, 'e1_hlt_phi')[0], 
                    getattr(out, 'e2_hlt_eta')[0], getattr(out, 'e2_hlt_phi')[0])

        tlv1 = TLorentzVector()
        tlv2 = TLorentzVector()
        
        tlv1.SetPtEtaPhiM(getattr(out, 'e1_hlt_pt'), getattr(out, 'e1_hlt_eta'), getattr(out, 'e1_hlt_phi'), me)
        tlv2.SetPtEtaPhiM(getattr(out, 'e2_hlt_pt'), getattr(out, 'e2_hlt_eta'), getattr(out, 'e2_hlt_phi'), me)
        
#        print('check', tlv1.E(), getattr(out, 'e1_hlt_energy'))

        getattr(out, 'hlt_eedr')[0] = dr 
        getattr(out, 'hlt_mee')[0] = (tlv1+tlv2).M()

    else:
        getattr(out, 'hlt_eedr')[0] = -1
        getattr(out, 'hlt_mee')[0] = -1
        
    out.gen_mass[0] = (epairs[0] + epairs[1]).M()
    out.gen_dr[0] = epairs[0].DeltaR(epairs[1])
        
        
        

#    # check exclusive dR matching ... 
#
#    print(match_index_l1,match_index_hlt)
#    for gindex in [1]:
#
#        l1index = -1
#        maxl1_dr = 99
#        
#        for jj in range(len(chain.eg_l1eg_et)):
#
#            if abs(chain.eg_l1eg_eta[jj]) > 1.2: continue
#
#            if jj==match_index_l1[0]: continue
#
#            dr = deltaR(chain.eg_gen_eta[gindex], chain.eg_gen_phi[gindex],
#                        chain.eg_l1eg_eta[jj], chain.eg_l1eg_phi[jj])
#
#            if maxl1_dr > dr:
#
#                maxl1_dr = dr
#                l1index = jj
#                
#
#
#
#        getattr(out, 'gen_e' + str(gindex+1) + '_l1_exclusive_dr')[0] = maxl1_dr
#
#
#        hltindex = -1
#        maxhlt_dr = 99
#
#        for jj in range(len(chain.eg_et)):
#
#            if abs(chain.eg_eta[jj]) > 1.2: continue
#            if not hlt_criteria(chain, jj): continue
#
#            if jj==match_index_hlt[0]: continue
#
#            dr = deltaR(chain.eg_gen_eta[gindex], chain.eg_gen_phi[gindex],
#                        chain.eg_eta[jj], chain.eg_phi[jj]) 
#            
#           
#            if maxhlt_dr > dr:
#                
#                maxhlt_dr = dr
#                hltindex = jj
#
#
#        getattr(out, 'gen_e' + str(gindex+1) + '_hlt_exclusive_dr')[0] = maxhlt_dr




        

    #out.isgjson[0] = chain.isgjson
    #out.instL[0] = chain.instL
    #out.npu[0] = chain.npu

    out.tree.Fill()


    evtid += 1

    
print(Nevt, 'evt processed.', evtid, 'evt has matching')

out.endJob()
