import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaPhi

class WbosonReconstruction(Module):
    def __init__(self,
        leptonObject = lambda event: Collection(event, "Muon")[0],
        metObject =lambda event: Object(event, "MET"),
        globalOptions={"isData":False}, 
        outputName='nominal'
    ):
        self.leptonObject = leptonObject
        self.metObject = metObject
        self.globalOptions=globalOptions
        self.outputName=outputName
        
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.outputName+"_mtw", "F")
        self.out.branch(self.outputName+"_deltaPhi", "F")
            
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        lepton = self.leptonObject(event)
        met = self.metObject(event)
        
        dPhi = deltaPhi(lepton,met)
        mtw = math.sqrt(2*lepton.pt*met.pt*(1-math.cos(dPhi)))

        self.out.fillBranch(self.outputName+"_mtw", mtw)
        self.out.fillBranch(self.outputName+"_deltaPhi", dPhi)
        
        return True
            
            
