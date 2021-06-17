import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class ScaleUncertainty(Module):
    def __init__(
        self,
        xsecs,
        storeVariables = []
    ):
        self.xsecs = xsecs
        self.storeVariables = storeVariables
        
        
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        processName = None
        for process in self.xsecs.keys():
            if inputFile.GetName().find(process)>=0:
                if processName != None:
                    raise Exception("Process duplications: '%s' vs '%s' for file '%s'"%(processName,process,inputFile.GetName()))
                processName = process
        if processName == None:
            raise Exception("Process xsec not found for file '%s'"%(inputFile.GetName()))

        xsec = self.xsecs[processName]['weights']['1']['xsec']
        self.inclUp = xsec['up']/xsec['nominal']
        self.inclDown = xsec['down']/xsec['nominal']
         
        #total scale variation (NB: normalization is also part of theo. xsec so should not be considered again)   
        self.out.branch("scale_incl_nominal","F")
        self.out.branch("scale_incl_up","F")
        self.out.branch("scale_incl_down","F")
        
        #only residual shape variation wrt. uncertainty on normalization
        self.out.branch("scale_shapeonly_nominal","F")
        self.out.branch("scale_shapeonly_up","F")
        self.out.branch("scale_shapeonly_down","F")
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        nominal = getattr(event,"LHEWeights_murNominal_mufNominal_1")
        up = getattr(event,"LHEWeights_murUp_mufUp_1")
        down = getattr(event,"LHEWeights_murDown_mufDown_1")
    
        self.out.fillBranch("scale_incl_nominal",1)
        self.out.fillBranch("scale_incl_up",up)
        self.out.fillBranch("scale_incl_down",down)
        
        self.out.fillBranch("scale_shapeonly_nominal",1)
        self.out.fillBranch("scale_shapeonly_up",up/self.inclUp)
        self.out.fillBranch("scale_shapeonly_down",down/self.inclDown)
        
        return True
        
