import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class LepJetFinder(Module):
    def __init__(
        self,
        jetCollection,
        leptonCollection,
        outputName="lepJet",
    ):
        self.jetCollection = jetCollection
        self.leptonCollection = leptonCollection
        self.outputName = outputName

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.outputName+"_jetIdx",
                        "I", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_pt",
                        "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_eta",
                        "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_phi",
                        "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_deltaR",
                        "F", lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        jetCollection = self.jetCollection(event)
        leptonCollection = self.leptonCollection(event)
        jet_idxs = []
        jet_pts = []
        jet_etas = []
        jet_phis = []
        jet_deltaRs = []
        lepJets = []

        for lepton in leptonCollection:
            jet = jetCollection[0]
            deltaR = lepton.p4().DeltaR(jet.p4())
            for _jet in jetCollection:
                _deltaR = lepton.p4().DeltaR(_jet.p4())
                if _deltaR < deltaR:
                    jet = _jet
                    deltaR = _deltaR

            lepJets.append(jet)
            jet_deltaRs.append(lepton.p4().DeltaR(jet.p4()))
            jet_pts.append(jet.pt)
            jet_etas.append(jet.eta)
            jet_phis.append(jet.phi)
            jet_idxs.append(jet._index)

        self.out.fillBranch(self.outputName+"_jetIdx", jet_idxs)
        self.out.fillBranch(self.outputName+"_pt", jet_pts)
        self.out.fillBranch(self.outputName+"_eta", jet_etas)
        self.out.fillBranch(self.outputName+"_phi", jet_phis)
        self.out.fillBranch(self.outputName+"_deltaR", jet_deltaRs)
        setattr(event, self.outputName, lepJets)

        return True
