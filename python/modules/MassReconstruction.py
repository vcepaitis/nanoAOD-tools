import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from utils import deltaR, deltaPhi


class MassReconstruction(Module):
    def __init__(
        self,
        globalOptions={"isData":False, "isSignal":False},
        outputName="category",
        tightLeptons=None,
        looseLeptons=None,
        jets=lambda event: Collection(event, "Jet"),

    ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.tightLeptons = tightLeptons
        self.looseLeptons = looseLeptons
        self.jets = jets

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self.outputName+"_m_llj", "F")
        self.out.branch(self.outputName+"_deltaPhil1j", "F")
        self.out.branch(self.outputName+"_deltaRl2j", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def analyze(self, event):

        tightLeptons = self.tightLeptons(event)
        looseLeptons = self.looseLeptons(event)

        jets = self.jets(event)

        deltaR_l2j = -1.
        deltaPhi_l1j = -1.

        if len(tightLeptons) > 0 and len(jets) > 0:
            jets = sorted(jets, key=lambda jet: math.fabs(deltaPhi(jet, tightLeptons[0])), reverse=True)
            deltaPhi_l1j = math.fabs(deltaPhi(tightLeptons[0], jets[0]))

        if len(looseLeptons) > 0 and len(jets) > 0:
            jets = sorted(jets, key=lambda jet: deltaR(jet, looseLeptons[0]), reverse=False)
            deltaR_l2j = deltaR(looseLeptons[0], jets[0])

        jetHNLCandidates = []
        if len(jets) > 0:
            jetHNLCandidates.append(jets[0]) 
            # Option to add second candidate

        WCandidateLorentzVector = ROOT.TLorentzVector(0,0,0,0)

        for lepton in tightLeptons:
            WCandidateLorentzVector += lepton.p4()

        for jet in jetHNLCandidates:
            WCandidateLorentzVector += jet.p4()
            for lepton in looseLeptons:
                if deltaR(jet, lepton) < 0.4:
                    WCandidateLorentzVector += lepton.p4()

        self.out.fillBranch(self.outputName+"_m_llj", WCandidateLorentzVector.M())
        self.out.fillBranch(self.outputName+"_deltaPhil1j", deltaPhi_l1j)
        self.out.fillBranch(self.outputName+"_deltaRl2j", deltaR_l2j)
        setattr(event, self.outputName+"_jets", jetHNLCandidates)

        return True
