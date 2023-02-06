import ROOT
import math 
from TreeProducerCommon import *


class TreeProducer_eff(TreeProducerCommon):
    """Class to create a custom output file & tree; as well as create and contain branches."""

    def __init__(self, name, **kwargs):
        print('TreeProducer is called for', name)
        super(TreeProducer_eff, self).__init__(name,**kwargs)


        for name in ['e1', 'e2']:
            self.addBranch('gen_' + name + '_pt',                  'f')
            self.addBranch('gen_' + name + '_eta',                  'f')
            self.addBranch('gen_' + name + '_phi',                  'f')
            self.addBranch('gen_' + name + '_vz',                  'f')

            self.addBranch('gen_' + name + '_l1_dr',                  'f')
            self.addBranch('gen_' + name + '_l1_dphi',                  'f')
            self.addBranch('gen_' + name + '_l1_deta',                  'f')
            self.addBranch('gen_' + name + '_hlt_dr',                  'f')
            self.addBranch('gen_' + name + '_hlt_dphi',                  'f')
            self.addBranch('gen_' + name + '_hlt_deta',                  'f')


            self.addBranch(name + '_l1_pt',                  'f')
            self.addBranch(name + '_l1_eta',                  'f')
            self.addBranch(name + '_l1_phi',                  'f')

            self.addBranch(name + '_hlt_pt',                  'f')
            self.addBranch(name + '_hlt_eta',                  'f')
            self.addBranch(name + '_hlt_phi',                  'f')
            self.addBranch(name + '_hlt_energy', 'f')
            self.addBranch(name + '_hlt_rawEnergy', 'f')
            self.addBranch(name + '_hlt_phiWidth', 'f')
            self.addBranch(name + '_hlt_nrClus', 'i')
            self.addBranch(name + '_hlt_seedId', 'f')
            self.addBranch(name + '_hlt_seedDet', 'i')
            self.addBranch(name + '_hlt_seednrCrystals', 'i')
            self.addBranch(name + '_hlt_sigmaIEtaIEta', 'f')
            self.addBranch(name + '_hlt_sigmaIEtaIEtaNoise', 'f')
            self.addBranch(name + '_hlt_ecalPFIsol', 'f')
            self.addBranch(name + '_hlt_hcalPFIsol', 'f')
            self.addBranch(name + '_hlt_trkIsol', 'f')
            self.addBranch(name + '_hlt_trkChi2', 'f')
            self.addBranch(name + '_hlt_trkMissHits', 'f')
            self.addBranch(name + '_hlt_trkNrLayerIT', 'i')
            self.addBranch(name + '_hlt_trkValidHits', 'i')
            self.addBranch(name + '_hlt_hcalHForHoverE', 'f')
            self.addBranch(name + '_hlt_invESeedInvP', 'f')
            self.addBranch(name + '_hlt_invEInvP', 'f')
            self.addBranch(name + '_hlt_trkDEta', 'f')
            self.addBranch(name + '_hlt_trkDEtaSeed', 'f')
            self.addBranch(name + '_hlt_trkDPhi', 'f')
            self.addBranch(name + '_hlt_pms2', 'f')
            self.addBranch(name + '_hlt_dxy', 'f')
            self.addBranch(name + '_hlt_vz', 'f')
            

        
        self.addBranch('l1_eedr',                  'f')
        self.addBranch('hlt_eedr',                  'f')
        self.addBranch('l1_mee',                  'f')
        self.addBranch('hlt_mee',                  'f')

        self.addBranch('gen_mass',                  'f')
        self.addBranch('gen_dr',                  'f')

        self.addBranch('isgjson', '?')
        self.addBranch('instL', 'f')
        self.addBranch('npu', 'i')

