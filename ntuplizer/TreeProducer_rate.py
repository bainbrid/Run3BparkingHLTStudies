import ROOT
import math 
from TreeProducerCommon import *


#ptrange = np.arange(5, 12, 0.5).tolist() 


class TreeProducer_rate(TreeProducerCommon):
    """Class to create a custom output file & tree; as well as create and contain branches."""

    def __init__(self, name, **kwargs):
        print('TreeProducer_rate is called for', name)
        super(TreeProducer_rate, self).__init__(name,**kwargs)



        for pt in l1_ptrange:
#            for dr in drrange:
            #            self.addBranch('doubleE' + str(pt) + '_nomatch',                  '?')
#            self.addBranch('l1_doubleE' + '{0:.1f}'.format(pt).replace('.','p') + '_dR' + '{0:.1f}'.format(dr).replace('.','p'),                  '?')
            self.addBranch('l1_doubleE' + '{0:.1f}'.format(pt).replace('.','p'),                  '?')
                


        for pt in hlt_ptrange:
#            for dr in drrange:
                #            self.addBranch('doubleE' + str(pt),                  '?')
#            self.addBranch('doubleE' + '{0:.1f}'.format(pt).replace('.','p') + '_dR' + '{0:.1f}'.format(dr).replace('.','p'),                  '?')
            self.addBranch('doubleE' + '{0:.1f}'.format(pt).replace('.','p'),                  '?')

            self.addBranch('mass_doubleE' + '{0:.1f}'.format(pt).replace('.','p'),                  '?')

            setattr(self, 'dr_' + str(pt), ROOT.TH1F('dr_' + str(pt), 'dr_' + str(pt), 1000,0,2*math.pi))
            setattr(self, 'mass_' + str(pt), ROOT.TH1F('mass_' + str(pt), 'mass_' + str(pt), 1000,0,20))


        self.addBranch('isgjson', '?')
        self.addBranch('instL', 'f')
        self.addBranch('npu', 'i')
        self.addBranch('run', 'i')
        self.addBranch('ls', 'i')
        self.addBranch('evt', 'i')
        self.addBranch('bx', 'i')
