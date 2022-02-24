import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from utils import deltaR, deltaPhi


class MassReconstruction_electron_sf(Module):
    def __init__(
        self,
        globalOptions={"isData":False, "isSignal":False},
        outputName="category",
        tightLeptons=None,
        looseLeptons=None,

    ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.tightLeptons = tightLeptons
        self.looseLeptons = looseLeptons

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self.outputName+"_m_lll", "F")
        self.out.branch(self.outputName+"_minDeltaPhi_l1l3", "F")
        self.out.branch(self.outputName+"_minDeltaPhi_l2l3", "F")
        self.out.branch(self.outputName+"_minDeltaR_l1l3", "F")
        self.out.branch(self.outputName+"_minDeltaR_l2l3", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

                tightLeptons = []
        looseLeptons = []

        tightLeptons = tightMuons+tightElectrons
        looseLeptons = looseMuons+looseElectrons

        tightLeptons = sorted(tightLeptons, key=lambda x: x.pt, reverse=True)
        
        # select leading only, move subleading to "loose"
        looseLeptons.extend(tightLeptons[2:])
        tightLeptons = tightLeptons[:2]
        looseLeptons = sorted(looseLeptons, key=lambda x: x.pt, reverse=True)

        tightLeptons = self.tightLeptons(event)
        looseLeptons = self.looseLeptons(event)

        jets = self.jets(event)

        deltaR_l1l3 = -1.
        deltaR_l2l3 = -1.
        deltaPhi_l1l3 = -1.
        deltaPhi_l2l3 = -1.

        if len(tightLeptons) > 1 and len(looseLeptons) > 0:
            deltaPhi_l1l3 = math.fabs(deltaPhi(tightLeptons[0], looseLeptons[0]))
            deltaPhi_l2l3 = math.fabs(deltaPhi(tightLeptons[1], looseLeptons[0]))
            deltaR_l1l3 = deltaR(tightLeptons[0], looseLeptons[0])
            deltaR_l2l3 = deltaR(tightLeptons[1], looseLeptons[0])

        WCandidateLorentzVector = ROOT.TLorentzVector(0,0,0,0)
        WCandidateLorentzVector = tightLeptons[0].p4() + tightLeptons[1].p4() + looseLeptons[0].p4()

        self.out.fillBranch(self.outputName+"_m_lll", WCandidateLorentzVector.M())
        self.out.fillBranch(self.outputName+"_minDeltaPhi_l1l3", deltaPhi_l1l3)
        self.out.fillBranch(self.outputName+"_minDeltaPhi_l2l3", deltaPhi_l2l3)
        self.out.fillBranch(self.outputName+"_minDeltaR_l1l3", deltaR_l1l3)
        self.out.fillBranch(self.outputName+"_minDeltaR_l2l3", deltaR_l2l3)

        return True
