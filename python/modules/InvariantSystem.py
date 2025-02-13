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
        self.out.branch(self.outputName+"_charge", "I")
        
        self.out.branch(self.outputName+"_dRmin", "F")
        self.out.branch(self.outputName+"_dPhimin", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        collection = self.inputCollection(event)
        vec = ROOT.TLorentzVector()
        charge = 1
        
        dRmin = 100.
        dPhimin = 100.
        for i,obj in enumerate(collection):
            vec += obj.p4()
            charge *= obj.charge
            for j,obj2 in enumerate(collection):
                if i>=j:
                    continue
                dRmin = min(dRmin,deltaR(obj,obj2))
                dPhimin = min(dPhimin,math.fabs(deltaPhi(obj,obj2)))
            
        self.out.fillBranch(self.outputName+"_mass", vec.M())
        self.out.fillBranch(self.outputName+"_pt", vec.Pt())
        self.out.fillBranch(self.outputName+"_eta", vec.Eta())
        self.out.fillBranch(self.outputName+"_charge", charge)
        
        self.out.fillBranch(self.outputName+"_dRmin", dRmin)
        self.out.fillBranch(self.outputName+"_dPhimin", dPhimin)

        return True
        
