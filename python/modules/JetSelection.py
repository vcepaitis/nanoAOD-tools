import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi


class JetSelection(Module):

    LOOSE = 0
    TIGHT = 1
    TIGHTLEPVETO = 2
    
    def __init__(
         self,
         inputCollection=lambda event: Collection(event, "Jet"),
         leptonCollectionDRCleaning=lambda event: [],
         leptonCollectionP4Subraction=lambda event: [],
         outputName="selectedJets",
         jetMinPt=15.,
         jetMinEta=-1.,
         jetMaxEta=2.4,
         jetMinNConstituents=-1,
         dRCleaning=0.4,
         dRP4Subtraction=0.4,
         flagDA=False,
         storeKinematics=['pt', 'eta', 'phi', 'minDeltaRSubtraction', 'ptLepton', 'ptSubtracted'],
         globalFeatures = [],
         globalOptions={"isData": False, "year": 2016},
         jetId=TIGHT
     ):
        self.globalOptions = globalOptions

        self.inputCollection = inputCollection
        self.leptonCollectionDRCleaning = leptonCollectionDRCleaning
        self.leptonCollectionP4Subraction = leptonCollectionP4Subraction
        self.outputName = outputName
        self.jetMinPt = jetMinPt
        self.jetMinEta = jetMinEta
        self.jetMaxEta = jetMaxEta
        self.jetMinNConstituents = jetMinNConstituents
        self.dRCleaning = dRCleaning
        self.dRP4Subtraction = dRP4Subtraction
        self.flagDA = flagDA
        self.storeKinematics = storeKinematics
        self.globalFeatures = globalFeatures
        if jetId==JetSelection.LOOSE and (globalOptions["year"] == 2017 or globalOptions["year"] == 2018):
            self.jetId = JetSelection.TIGHT
        else:
            self.jetId = jetId

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.flagDA:
            self.out.branch("Jet_forDA", "F", lenVar="nJet")

        self.out.branch("n"+self.outputName, "I")
        for variable in self.storeKinematics+self.globalFeatures:
            self.out.branch(self.outputName+"_"+variable, "F", lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)
        
        
        jetglobal = Collection(event, "global")
        jetglobal_indices = [global_jet.jetIdx for global_jet in jetglobal]
        
        selectedJets = []
        unselectedJets = []

        leptonsForDRCleaning = self.leptonCollectionDRCleaning(event)
        leptonsForP4Subtraction = self.leptonCollectionP4Subraction(event)

        if self.flagDA:
            flagsDA = [0.]*event.nJet

        for jet in jets:            
            #find global jet to access more properties
            global_jet = None
            if len(self.globalFeatures)>0 and math.fabs(jet.eta)<2.4:
                try:
                    global_jet_index = jetglobal_indices.index(jet._index)
                    global_jet = jetglobal[global_jet_index]
                    if abs(jet.eta - global_jet.eta) > 0.01 or \
                       abs(jet.phi - global_jet.phi) > 0.01:
                           print "Warning ->> jet might be mismatched!"
                except ValueError:
                    print "WARNING: jet (pt: %s, eta: %s) does not have a matching global jet" % (jet.pt, jet.eta)
                    
            if global_jet:
                for feature in self.globalFeatures:
                    setattr(jet,feature,getattr(global_jet,feature))
            else:
                for feature in self.globalFeatures:
                    setattr(jet,feature,0.0)
  
            if jet.pt<self.jetMinPt:
                unselectedJets.append(jet)
                continue   

            if math.fabs(jet.eta) > self.jetMaxEta:
                unselectedJets.append(jet)
                continue
                
            if (self.jetMinEta>0.) and (math.fabs(jet.eta) < self.jetMinEta):
                unselectedJets.append(jet)
                continue
                
            if (jet.jetId & (1 << self.jetId)) == 0:
                unselectedJets.append(jet)
                continue
                
            if self.jetMinNConstituents > 0 and jet.nConstituents < self.jetMinNConstituents:
                unselectedJets.append(jet.nConstituents)
                continue

            minDeltaRSubtraction = 999.

            leptonP4 = ROOT.TLorentzVector(0,0,0,0)
            if len(leptonsForP4Subtraction) > 0:
                for lepton in leptonsForP4Subtraction:
                    minDeltaRSubtraction = min(minDeltaRSubtraction, deltaR(lepton, jet))
                    if deltaR(lepton,jet)<self.dRP4Subtraction:
                        leptonP4 += lepton.p4()
         
            setattr(jet,"minDeltaRSubtraction", minDeltaRSubtraction)
            leptonPt = leptonP4.Pt()
            setattr(jet,"ptLepton", leptonPt)
            jetP4LeptonSubtracted = jet.p4()-leptonP4 
            setattr(jet,"ptSubtracted", jetP4LeptonSubtracted.Pt())

            if len(leptonsForDRCleaning) > 0:
                mindphi = min(map(lambda lepton: math.fabs(deltaPhi(lepton, jet)), leptonsForDRCleaning))
                mindr = min(map(lambda lepton: deltaR(lepton, jet), leptonsForDRCleaning))
                
                if mindr < self.dRCleaning:
                    unselectedJets.append(jet)
                    continue
                    
                setattr(jet,"minDPhiClean",mindphi)
                setattr(jet,"minDRClean",mindr)
            else:
                setattr(jet,"minDPhiClean",100)
                setattr(jet,"minDRClean",100)

            selectedJets.append(jet)

            if self.flagDA:
                flagsDA[jet._index] = 1.

        if self.flagDA:
            self.out.fillBranch("Jet_forDA", flagsDA)

        self.out.fillBranch("n"+self.outputName, len(selectedJets))
        for variable in self.storeKinematics+self.globalFeatures:
            self.out.fillBranch(self.outputName+"_"+variable,
                                map(lambda jet: getattr(jet, variable), selectedJets))

        setattr(event, self.outputName, selectedJets)
        setattr(event, self.outputName+"_unselected", unselectedJets)

        return True

