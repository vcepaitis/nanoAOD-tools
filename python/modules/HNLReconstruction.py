import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from utils import deltaR, deltaPhi


class HNLReconstruction(Module):
    def __init__(
        self,
        lepton1Object = lambda event: event.leadingLeptons[0],
        lepton2Object = None,
        jetCollection=lambda event: Collection(event, "Jet"),
        outputName="nominal",
        globalOptions={"isData":False, "isSignal":False},

    ):
        self.lepton1Object = lepton1Object
        self.lepton2Object = lepton2Object
        self.jetCollection = jetCollection

        self.outputName = outputName
        self.globalOptions = globalOptions

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        if self.lepton2Object==None:
            self.out.branch(self.outputName+"_m_l1j", "F")
            self.out.branch(self.outputName+"_deltaPt_l1j", "F")
            self.out.branch(self.outputName+"_deltaPhi_l1j", "F")
            self.out.branch(self.outputName+"_deltaR_l1j", "F")
            self.out.branch(self.outputName+"_deltaEta_l1j", "F")
        else:
            self.out.branch(self.outputName+"_m_llj", "F")
            self.out.branch(self.outputName+"_deltaPt_llj", "F")
            self.out.branch(self.outputName+"_deltaPhi_l1j", "F")
            self.out.branch(self.outputName+"_deltaEta_l1j", "F")
            self.out.branch(self.outputName+"_deltaR_l2j", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def analyze(self, event):

        lepton1 = self.lepton1Object(event)
        lepton2 = None if self.lepton2Object==None else self.lepton2Object(event)
          
        jets = self.jetCollection(event)

        if lepton2==None:
            if len(jets) > 0:
                #take jet opposite of l1 in transverse plane
                hnlJet = sorted(jets, key=lambda jet: math.fabs(deltaPhi(jet, lepton1)), reverse=True)[0]
                self.out.branch(self.outputName+"_m_l1j", (lepton1.p4()+hnlJet.p4()).M())
                self.out.branch(self.outputName+"_deltaPt_l1j", (lepton1.p4()-hnlJet.p4()).Pt())
                self.out.branch(self.outputName+"_deltaPhi_l1j", math.fabs(deltaPhi(lepton1, hnlJet)))
                self.out.branch(self.outputName+"_deltaR_l1j", deltaR(lepton1, hnlJet))
                self.out.branch(self.outputName+"_deltaEta_l1j", math.fabs(lepton1.eta-hnlJet.eta))
                setattr(event, self.outputName+"_hnlJets", [hnlJet])
            else:
                self.out.branch(self.outputName+"_m_lj", 0)
                self.out.branch(self.outputName+"_deltaPt_lj", 0)
                self.out.branch(self.outputName+"_deltaPhi_l1j", 0)
                self.out.branch(self.outputName+"_deltaR_l1j", 0)
                self.out.branch(self.outputName+"_deltaEta_l1j", 0)
                setattr(event, self.outputName+"_hnlJets", [])
        else:
            if len(jets) > 0:
                #take jet closest to l2
                hnlJet = sorted(jets, key=lambda jet: deltaR(jet, lepton2), reverse=False)[0]
                self.out.branch(self.outputName+"_m_llj", (lepton1.p4()+lepton2.p4()+hnlJet.p4()).M())
                self.out.branch(self.outputName+"_deltaPt_llj", (lepton1.p4()-lepton2.p4()-hnlJet.p4()).Pt())
                self.out.branch(self.outputName+"_deltaPhi_l1j", math.fabs(deltaPhi(lepton1,hnlJet)))
                self.out.branch(self.outputName+"_deltaEta_l1j", math.fabs(lepton1.eta-hnlJet.eta))
                self.out.branch(self.outputName+"_deltaR_l2j", deltaR(lepton2,hnlJet))
                setattr(event, self.outputName+"_hnlJets", [hnlJet])
                
            else:
                self.out.branch(self.outputName+"_m_llj", 0)
                self.out.branch(self.outputName+"_deltaPt_llj", 0)
                self.out.branch(self.outputName+"_deltaPhi_l1j", 0)
                self.out.branch(self.outputName+"_deltaEta_l1j", 0)
                self.out.branch(self.outputName+"_deltaR_l2j", 0)
                setattr(event, self.outputName+"_hnlJets", [])
        
        return True
