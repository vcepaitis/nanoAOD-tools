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
        isSignal = True,
    ):
        self.xsecs = xsecs
        self.isSignal = isSignal
        
        
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.processName = None
        
        if self.isSignal:
            for process in self.xsecs.keys():
                if inputFile.GetName().replace("pt20", "all").find(process)>=0:
                    if self.processName != None:
                        raise Exception("Process duplications: '%s' vs '%s' for file '%s'"%(self.processName,process,inputFile.GetName()))
                    self.processName = process
            if self.processName == None:
                raise Exception("Process xsec not found for file '%s'"%(inputFile.GetName()))
            processNameHack = self.processName.replace("pt20", "all")

            xsec = self.xsecs[processNameHack]['weights']['1']['xsec']
            self.inclUp = xsec['up']/xsec['nominal']
            self.inclDown = xsec['down']/xsec['nominal']
         
            #only residual shape variation wrt. uncertainty on normalization
            self.out.branch("scale_shapeonly_nominal","F")
            self.out.branch("scale_shapeonly_up","F")
            self.out.branch("scale_shapeonly_down","F")
            
        
        #total scale variation (NB: normalization is also part of theo. xsec so should not be considered again)   
        self.out.branch("scale_incl_nominal","F")
        self.out.branch("scale_incl_up","F")
        self.out.branch("scale_incl_down","F")
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        if self.isSignal:
            up = getattr(event,"LHEWeights_murUp_mufUp_1")
            down = getattr(event,"LHEWeights_murDown_mufDown_1")

        elif not hasattr(event, "nLHEScaleWeight"):
            down = 1
            up = 1
            
        elif event.nLHEScaleWeight==9:
            down  = getattr(event,"LHEScaleWeight")[0]
            up  = getattr(event,"LHEScaleWeight")[8]
            
        else:
            down = 1
            up = 1
            
        #sanity check
        if math.fabs(up-1)>1 or math.fabs(down-1)>1:
            up = 1
            down =1 
            
        if self.processName!=None:
            self.out.fillBranch("scale_shapeonly_nominal",1)
            self.out.fillBranch("scale_shapeonly_up",up/self.inclUp)
            self.out.fillBranch("scale_shapeonly_down",down/self.inclDown)
        
            
        self.out.fillBranch("scale_incl_nominal",1)
        self.out.fillBranch("scale_incl_up",up)
        self.out.fillBranch("scale_incl_down",down)

        return True
        
