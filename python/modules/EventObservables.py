import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel \
    import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi


class EventObservables(Module):

    def __init__(
        self,
        jetCollection=lambda event: Collection(event, "Jet"),
        metInput=lambda event: Object(event, "MET"),
        outputName="EventObservables",
        globalOptions={"isData": False}
    ):
        self.jetCollection = jetCollection
        self.metInput = metInput
        self.outputName = outputName
        self.globalOptions = globalOptions

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self.outputName+"_met", "F")
        self.out.branch(self.outputName+"_met_phi", "F")
        self.out.branch(self.outputName+"_ht", "F")
        self.out.branch(self.outputName+"_hmass", "F")
        self.out.branch(self.outputName+"_mht", "F")
        self.out.branch(self.outputName+"_mht_phi", "F")
        self.out.branch(self.outputName+"_mht_met_dphi", "F")
        self.out.branch(self.outputName+"_minPhiStar", "F")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        jets = self.jetCollection(event)

        met = self.metInput(event)
        self.out.fillBranch(self.outputName+"_met", met.pt)
        self.out.fillBranch(self.outputName+"_met_phi", met.phi)

        # ht and mht
        vectorSum = ROOT.TLorentzVector()
        scalarPtSum = 0.0

        for jet in jets:
            vectorSum += jet.p4Subtracted
            scalarPtSum += jet.ptSubtracted
            
        mht_met_dphi = math.fabs(deltaPhi(vectorSum.Phi(),met.phi))

        # minPhiStar
        minPhiStar = math.pi
        for jet in jets:
            negSum = -(vectorSum-jet.p4Subtracted)
            minPhiStar = min(minPhiStar, math.fabs(deltaPhi(negSum.Phi(), jet.phi)))

        self.out.fillBranch(self.outputName+"_minPhiStar", minPhiStar)
        self.out.fillBranch(self.outputName+"_ht", scalarPtSum)
        self.out.fillBranch(self.outputName+"_hmass", vectorSum.M())
        self.out.fillBranch(self.outputName+"_mht", vectorSum.Pt())
        self.out.fillBranch(self.outputName+"_mht_phi", vectorSum.Phi())
        self.out.fillBranch(self.outputName+"_mht_met_dphi", mht_met_dphi)
        
        setattr(event, self.outputName+"_minPhiStar", minPhiStar)
        setattr(event, self.outputName+"_ht", scalarPtSum)
        setattr(event, self.outputName+"_hmass", vectorSum.M())
        setattr(event, self.outputName+"_mht", vectorSum.Pt())
        setattr(event, self.outputName+"_mht_phi", vectorSum.Phi())
        setattr(event, self.outputName+"_mht_met_dphi", mht_met_dphi)


        return True
