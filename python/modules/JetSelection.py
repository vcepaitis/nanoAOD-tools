import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR

class JetSelection(Module):

    def __init__(
         self,
         inputCollection = lambda event: Collection(event, "Jet"),
         leptonCollection = lambda event: [],
         outputName = "selectedJets",
         jetMinPt = 15.,
         jetMaxPt = 100.,
         jetMaxEta = 2.4,
         dRCleaning = 0.4,
         minNconstituents = 3,
         genJetMinPt = 0.,
         minRatio = 0.,
         flagDA = False,
         addSize = True,
         storeKinematics=['pt','eta'],
         globalOptions={"isData":False},
         jetId = -1
         ):
        self.globalOptions = globalOptions
        self.inputCollection = inputCollection
        self.leptonCollection = leptonCollection
        self.outputName = outputName
        self.jetMinPt = jetMinPt
        self.jetMaxPt = jetMaxPt
        self.jetMaxEta = jetMaxEta
        self.dRCleaning = dRCleaning
        self.minNconstituents = minNconstituents
        self.genJetMinPt = genJetMinPt
        self.minRatio = minRatio
        self.flagDA = flagDA
        self.addSize = addSize
        self.storeKinematics = storeKinematics
        self.jetId = jetId

        
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.addSize:
            self.out.branch("n"+self.outputName,"I")
            self.out.branch("nfailedId"+self.outputName,"I")
        if self.flagDA:
            self.out.branch(self.outputName+"_forDA","F",lenVar="nJet")
        
        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
                
                
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets = self.inputCollection(event)
        origJets = Collection(event,"Jet")
        if not self.globalOptions["isData"]:
            genJets = Collection(event,"GenJet")
        selectedJets = []
        unselectedJets = []
        
        if self.flagDA:
            flagsDA = [0.]*event.nJet
        
        failedId = 0
        
        for jet in jets:
            if jet.pt>self.jetMinPt and jet.pt<self.jetMaxPt and math.fabs(jet.eta)<self.jetMaxEta and (jet.jetId>self.jetId) and (jet.nConstituents>self.minNconstituents):
                leptons = self.leptonCollection(event)
                if self.dRCleaning>0. and leptons!=None and len(leptons)>0:
                    mindr = min(map(lambda lepton: deltaR(lepton,jet),leptons))
                    if mindr<self.dRCleaning:
                        unselectedJets.append(jet)
                        continue
            else:
                unselectedJets.append(jet)
                continue
                        
            if self.flagDA:
                flagsDA[jet._index]=1.
            if not self.globalOptions["isData"]:
                if jet.genJetIdx == -1:
                    selectedJets.append(jet)
                    continue
                elif jet.genJetIdx >= len(genJets):
                    unselectedJets.append(jet)
                    continue
                else:
                    genJet = genJets[jet.genJetIdx]
                if genJet.pt>self.genJetMinPt and jet.pt/genJet.pt>self.minRatio:
                    selectedJets.append(jet)
                else:
                    unselectedJets.append(jet)
            else:
                selectedJets.append(jet)

        if self.addSize:
            self.out.fillBranch("n"+self.outputName,len(selectedJets))
            self.out.fillBranch("nfailedId"+self.outputName,failedId)
        if self.flagDA:
            self.out.fillBranch(self.outputName+"_forDA",flagsDA)
        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName+"_"+variable,map(lambda jet: getattr(jet,variable),selectedJets))
        
        setattr(event,self.outputName,selectedJets)
        setattr(event,self.outputName+"_unselected",unselectedJets)
        
        return True
    
