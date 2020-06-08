import os
import sys
import math
import json
import ROOT
import random

from utils import deltaPhi, deltaR

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class InvariantSystem(Module):
    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "Muon"),
        outputName="sys",
        globalOptions={"isData": False}
    ):
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.globalOptions = globalOptions

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.outputName+"_mass", "F")
        self.out.branch(self.outputName+"_pt", "F")
        self.out.branch(self.outputName+"_eta", "F")
        self.out.branch(self.outputName+"_minDeltaR", "F")
        self.out.branch(self.outputName+"_minDeltaPhi", "F")
        self.out.branch(self.outputName+"_maxDeltaR", "F")
        self.out.branch(self.outputName+"_maxDeltaPhi", "F")
        self.out.branch(self.outputName+"_charge", "I")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        collection = self.inputCollection(event)
        vec = ROOT.TLorentzVector()
        charge = 1
        
        minDeltaR = 100.
        maxDeltaR = 0.
        minDeltaPhi = 100.
        maxDeltaPhi = 0.
        
        for i,obj in enumerate(collection):
            vec += obj.p4()
            charge *= obj.charge
            
            for j,obj2 in enumerate(collection):
                if j>=i:
                    continue
                minDeltaR = min(minDeltaR,deltaR(obj,obj2))
                maxDeltaR = max(maxDeltaR,deltaR(obj,obj2))
                minDeltaPhi = min(minDeltaPhi,deltaPhi(obj,obj2))
                maxDeltaPhi = max(maxDeltaPhi,deltaPhi(obj,obj2))
    
        self.out.fillBranch(self.outputName+"_mass", vec.M())
        self.out.fillBranch(self.outputName+"_pt", vec.Pt())
        self.out.fillBranch(self.outputName+"_eta", vec.Eta())
        self.out.fillBranch(self.outputName+"_minDeltaR",  minDeltaR)
        self.out.fillBranch(self.outputName+"_minDeltaPhi",  minDeltaPhi)
        self.out.fillBranch(self.outputName+"_maxDeltaR",  maxDeltaR)
        self.out.fillBranch(self.outputName+"_maxDeltaPhi",  maxDeltaPhi)
        self.out.fillBranch(self.outputName+"_charge", charge)

        return True
