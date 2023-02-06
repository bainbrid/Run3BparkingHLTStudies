#import numpy as num
import os, math, sys
#from ROOT import TFile, TTree, TH1F, gROOT, TTree, Double, TChain, TLorentzVector, TVector3
import numpy as num
from array import array
import ROOT

root_dtype = {
  float: 'D',  int: 'I',  bool: 'O',
  'f':   'D',  'i': 'I',  '?':  'O',  'b': 'b', 'v':'D', 's':'s'
}

num_dtype = {
  'D':   'f',  'I': 'i',  'O':  '?',  'b': 'b'
}


#ptrange = num.arange(5, 15, 1).tolist() 
l1_ptrange = num.arange(4, 12, 0.5).tolist() 
hlt_ptrange = num.arange(4, 12, 0.5).tolist() 
#drrange = num.arange(0.3, 1.1, 0.1).tolist() 

class TreeProducerCommon(object):
    """Class to create a custom output file & tree; as well as create and contain branches."""
    
    def __init__(self, name, **kwargs):
        print('TreeProducerCommon is called for', name)
        
        self.name       = name
        
        # TREE
        self.outputfile = ROOT.TFile(name, 'RECREATE')
        self.tree       = ROOT.TTree('tree','tree')
        self.v = ROOT.std.vector( ROOT.std.string )()
        self.tree._v = self.v


        ###################################
#        self.hist = ROOT.TH1F('cutflow', 'cutflow', 15,0,15)

#        self.addBranch('b_eta',                  'f')


    def addBranch(self, name, dtype='f', default=None):
        """Add branch with a given name, and create an array of the same name as address."""
        if hasattr(self,name):
          print("ERROR! TreeProducerCommon.addBranch: Branch of name '%s' already exists!"%(name))
          exit(1)
        if isinstance(dtype,str):
          if dtype.lower()=='f': # 'f' is only a 'float32', and 'F' is a 'complex64', which do not work for filling float branches
            dtype = float        # float is a 'float64' ('f8')
          elif dtype.lower()=='i': # 'i' is only a 'int32'
            dtype = int            # int is a 'int64' ('i8')


#        if dtype=='v':
#            setattr(self,name,array('d', 1000*[0]))
#            setattr(self,name,array('d', 1000,dtype='f'))
#            getattr(self, treename).Branch(name, getattr(self,name), '%s[1000]/%s'%(name,root_dtype[dtype]))
#        elif dtype=='s':

#            self.tree.Branch( name, self.v )

#            setattr(self,name,None)

#            self.tree.Branch(name, getattr(self,name), '%s'%(name))
#        else:
        setattr(self,name,num.zeros(1,dtype=dtype))
        self.tree.Branch(name, getattr(self,name), '%s/%s'%(name,root_dtype[dtype]))

        if default!=None:
          getattr(self,name)[0] = default
        
    
    def endJob(self):
        """Write and close files after the job ends."""
#        self.hist.Write()
#        self.hammer_hist.Write()
        self.outputfile.Write()
        self.outputfile.Close()
        

