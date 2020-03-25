import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi


class EventObservables(Module):

    def __init__(
        self,
        jetCollection=lambda event: Collection(event, "Jet"),
        leptonCollection=lambda event: Collection(event, "Muon"),
        metInput=lambda event: Object(event, "MET"),
        outputName="EventObservables",
        globalOptions={"isData": False}
    ):
        self.globalOptions = globalOptions
        self.jetCollection = jetCollection
        self.leptonCollection = leptonCollection
        self.metInput = metInput
        self.outputName = outputName

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.outputName+"_ht", "F")
        self.out.branch(self.outputName+"_mht", "F")
        self.out.branch(self.outputName+"_mht_phi", "F")
        self.out.branch(self.outputName+"_mass", "F")
        self.out.branch(self.outputName+"_minPhi", "F")

        if self.leptonCollection is not None:
            self.out.branch(self.outputName+"_MET_NoMu_pt", "F")
            self.out.branch(self.outputName+"_MET_NoMu_phi", "F")
            self.out.branch(self.outputName+"_MET_Mu_mT", "F")
            self.out.branch(self.outputName+"_mht_NoMu", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        objs = self.jetCollection(event)
        met = self.metInput(event)
        vectorSum = ROOT.TLorentzVector()
        scalarPtSum = 0.0
        for obj in objs:
            vectorSum += obj.p4()
            scalarPtSum += obj.pt

        self.out.fillBranch(self.outputName+"_ht", scalarPtSum)
        setattr(event, self.outputName+"_ht", scalarPtSum)
        self.out.fillBranch(self.outputName+"_mht", vectorSum.Pt())
        setattr(event, self.outputName+"_mht", vectorSum.Pt())
        self.out.fillBranch(self.outputName+"_mht_phi", vectorSum.Phi())
        setattr(event, self.outputName+"_mht_phi", vectorSum.Phi())
        self.out.fillBranch(self.outputName+"_mass", vectorSum.M())
        setattr(event, self.outputName+"_mass", vectorSum.M())

        if self.leptonCollection is not None:
            met_px_lc = met.pt*math.sin(met.phi)
            met_py_lc = met.pt*math.cos(met.phi)

        lepton = self.leptonCollection(event)
            met_px_lc += lepton.p4().Px()
            met_py_lc += lepton.p4().Py()
            met_lc = math.sqrt(met_px_lc**2+met_py_lc**2)
            self.out.fillBranch(self.outputName+"_MET_NoMu_pt", met_lc)
            setattr(event, self.outputName+"_MET_NoMu_pt", met_lc)

        met_l_mT = math.sqrt(2.*lepton.pt*met.pt*(1. - math.cos(lepton.phi - met.phi)))

        self.out.fillBranch(self.outputName+"_MET_Mu_mT", met_l_mT)
            setattr(event, self.outputName+"_MET_Mu_mT", met_l_mT)

        self.out.fillBranch(self.outputName+"_mht_NoMu", vectorSum.Pt() + lepton.p4().Pt())
            setattr(event, self.outputName+"_mht_NoMu", vectorSum.Pt() + lepton.p4().Pt())

        minPhi = math.pi
        for obj in objs:
            negSum = -(vectorSum-obj.p4())
            minPhi = min(minPhi, math.fabs(deltaPhi(negSum.Phi(), obj.phi)))
        self.out.fillBranch(self.outputName+"_minPhi", minPhi)
        setattr(event, self.outputName+"_minPhi", minPhi)
        return True
