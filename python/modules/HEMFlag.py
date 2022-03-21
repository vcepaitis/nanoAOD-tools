import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class HEMFlag(Module):

    def __init__(
        self,
        inputDict={"nominal": lambda event: event.hnlJet_nominal},
        leadingLeptons=lambda event: event.leadingLeptons,
        subleadingLeptons=lambda event: event.subleadingleptons
    ):
        self.inputDict = inputDict
        self.leadingLeptons = leadingLeptons
        self.subleadingLeptons = subleadingLeptons
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for key in self.inputDict.keys():
            self.out.branch("HEMFlag_"+key+"_jet","I")
            self.out.branch("HEMFlag_"+key+"_l1","I")
            self.out.branch("HEMFlag_"+key+"_l2","I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        leadingLeptons = self.leadingLeptons(event)
        subleadingLeptons = self.subleadingLeptons(event)

        for key, value in self.inputDict.iteritems():
            hnlJets = value(event)

            if len(leadingLeptons) > 0:
                leadingLepton = leadingLeptons[0]

                if -3.0 <leadingLepton.eta< -1.3 and -1.57 < leadingLepton.phi < -0.87:
                    self.out.fillBranch("HEMFlag_"+key+"_l1", 1)
                else:
                    self.out.fillBranch("HEMFlag_"+key+"_l1", 0)
            else:
                self.out.fillBranch("HEMFlag_"+key+"_l1", 0)

            if len(subleadingLeptons) > 0:
                subleadingLepton = subleadingLeptons[0]

                if -3.0 <subleadingLepton.eta< -1.3 and -1.57 < subleadingLepton.phi < -0.87:
                    self.out.fillBranch("HEMFlag_"+key+"_l2", 1)
                else:
                    self.out.fillBranch("HEMFlag_"+key+"_l2", 0)

            else:
                self.out.fillBranch("HEMFlag_"+key+"_l2", 0)

            if len(hnlJets) > 0:
                hnlJet = hnlJets[0]
                if -3.0 <hnlJet.eta< -1.3 and -1.57 < hnlJet.phi < -0.87:
                    self.out.fillBranch("HEMFlag_"+key+"_jet", 1)
                else:
                    self.out.fillBranch("HEMFlag_"+key+"_jet", 0)      
            else:
                self.out.fillBranch("HEMFlag_"+key+"_jet", 0)      

        return True
        
