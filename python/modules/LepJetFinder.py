import os
import sys
import math
import json
import ROOT
import random

import utils 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class LepJetFinder(Module):
    def __init__(
        self,
        jetCollection,
        leptonCollection,
        outputName="lepJet",
        storeKinematics=['pt', 'ptLeptonSubtracted', 'eta', 'phi', 'jetId', 'deltaR', 'nConstituents'],
    ):
        self.jetCollection = jetCollection
        self.leptonCollection = leptonCollection
        self.outputName = outputName
        self.storeKinematics = storeKinematics

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable, "F", lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        jetCollection = self.jetCollection(event)
        leptonCollection = self.leptonCollection(event)

        lepJets = []

        if len(jetCollection) == 0:
            setattr(event, self.outputName, lepJets)
            return True

        for lepton in leptonCollection:
            jet = jetCollection[0]
            deltaR = lepton.p4().DeltaR(jet.p4())
            for _jet in jetCollection:
                _deltaR = lepton.p4().DeltaR(_jet.p4())
                if _deltaR < deltaR:
                    jet = _jet
                    deltaR = _deltaR

            setattr(jet, "deltaR", deltaR)

            lepJets.append(jet)

        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName+"_"+variable,
                                map(lambda jet: getattr(jet, variable), lepJets))
        setattr(event, self.outputName, lepJets)

        return True
