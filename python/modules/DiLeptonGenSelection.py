import os
import sys
import math
import json
import ROOT
import random

from utils import deltaPhi, deltaR

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class DiLeptonGenSelection(Module):
    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "GenPart"),
        outputName="Gen_",
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


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        collection = self.inputCollection(event)
        index = -1
        muons = []
        electrons = []
        taus = []
        for obj in collection:
            if abs(obj.pdgId) == 11:
                electrons.append(obj)
            elif abs(obj.pdgId) == 13:
                muons.append(obj)
            elif abs(obj.pdgId) == 15:
                taus.append(obj)

        leptons = None

        if len(electrons) < 2 and len(muons) < 2:
            leptons = taus
        elif len(taus) < 2 and len(muons) < 2:
            leptons = electrons
        elif len(taus) < 2 and len(electrons) < 2:
            leptons = muons

        print(leptons)
        for lep in leptons:
            print lep.pt, lep.pdgId
            #if obj.pdgId, obj.genPartIdxMother

        #self.out.fillBranch(self.outputName+"_mass", vec.M())
        return True
