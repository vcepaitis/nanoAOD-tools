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
         jetMinPtMerged=15.,
         jetMinPt=15.,
         jetMinEta=-1.,
         jetMaxEta=2.4,
         jetMinNConstituents=-1,
         dRCleaning=0.4,
         dRP4Subtraction=0.4,
         flagDA=False,
         globalFeatures = [],
         storeKinematics=['pt', 'eta', 'phi', 'minDeltaRSubtraction', 'ptLepton', 'ptOriginal', 'ptSubtracted', 'rawFactor', 'ptRaw'],
         globalOptions={"isData": False, "year": 2016},
         jetId=TIGHT
     ):
        self.globalOptions = globalOptions

        self.inputCollection = inputCollection
        self.leptonCollectionDRCleaning = leptonCollectionDRCleaning
        self.leptonCollectionP4Subraction = leptonCollectionP4Subraction
        self.outputName = outputName
        self.jetMinPtMerged = jetMinPtMerged
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
            
        #reset jet pt in case modules are chained
        for jet in jets:
            if hasattr(jet,"p4Original"):
                jet.pt = jet.p4Original.Pt()
                jet.eta = jet.p4Original.Eta()
                jet.phi = jet.p4Original.Phi()
                jet.mass = jet.p4Original.M()

        for jet in jets:            
            #find global jet to access more properties
            global_jet = None
            if len(self.globalFeatures)>0 and math.fabs(jet.eta)<2.4:
                try:
                    global_jet_index = jetglobal_indices.index(jet._index)
                    global_jet = jetglobal[global_jet_index]
                    if math.fabs(jet.eta - global_jet.eta) > 0.01 or \
                       math.fabs(deltaPhi(jet.phi,global_jet.phi)) > 0.01:
                           print "Warning ->> jet might be mismatched! (phi: %.2f, eta: %.2f) != (phi: %.2f, eta: %.2f)"%(jet.phi, jet.eta, global_jet.phi, global_jet.eta)
                except ValueError:
                    print "WARNING: jet (pt: %s, eta: %s) does not have a matching global jet" % (jet.pt, jet.eta)
                    
            if global_jet:
                for feature in self.globalFeatures:
                    setattr(jet,feature,getattr(global_jet,feature))
            else:
                for feature in self.globalFeatures:
                    setattr(jet,feature,0.0)
  
        for jet in jets:    
            #undo jet uncertainty including nominal JER smearing
            uncFactor = jet.uncFactor if hasattr(jet,"uncFactor") else 1.
            jet.ptRaw = jet.pt*(1. - jet.rawFactor)/uncFactor
            jet.massRaw = jet.mass * (1 - jet.rawFactor) #mass not modified through uncertainties
            jet.p4Raw = ROOT.TLorentzVector()
            jet.p4Raw.SetPtEtaPhiM(jet.ptRaw,jet.eta,jet.phi,jet.massRaw)
            jet.p4Original = jet.p4()
            jet.ptOriginal = jet.pt
        
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
                unselectedJets.append(jet)
                continue

            minDeltaRSubtraction = 999.

            leptonP4 = ROOT.TLorentzVector(0,0,0,0)
            if len(leptonsForP4Subtraction) > 0:
                for lepton in leptonsForP4Subtraction:
                    minDeltaRSubtraction = min(minDeltaRSubtraction, deltaR(lepton, jet))
                    if deltaR(lepton,jet)<self.dRP4Subtraction:
                        leptonP4 += lepton.p4()
                      
            jet.p4Lepton = leptonP4    
            jet.p4Subtracted = jet.p4()-leptonP4 
            jet.p4OriginalSubtracted = jet.p4Original-leptonP4 
            jet.p4RawSubtracted = jet.p4Raw-leptonP4 
            
            jet.minDeltaRSubtraction = minDeltaRSubtraction

            if jet.minDeltaRSubtraction < 0.4:
                if jet.pt < self.jetMinPtMerged:
                    unselectedJets.append(jet)
                    continue    
            
            # if leptonP4.Pt()>1e-3 and jet.p4Subtracted.Pt()<20.:
            #     #reset jet pt & recaculate p4Subtracted
            #     jet.pt = jet.ptRaw 
            #     jet.mass = jet.massRaw 
            #     jet.p4Subtracted = jet.p4()-leptonP4
            
            jet.ptLepton = jet.p4Lepton.Pt() 
            jet.ptSubtracted = jet.p4Subtracted.Pt()
            jet.ptOriginalSubtracted = jet.p4OriginalSubtracted.Pt()
            jet.ptRawSubtracted = jet.p4RawSubtracted.Pt()

            if len(leptonsForDRCleaning) > 0:
                mindphi = min(map(lambda lepton: math.fabs(deltaPhi(lepton, jet)), leptonsForDRCleaning))
                mindr = min(map(lambda lepton: deltaR(lepton, jet), leptonsForDRCleaning))
                
                setattr(jet,"minDPhiClean",mindphi)
                setattr(jet,"minDRClean",mindr)
                
                if mindr < self.dRCleaning:
                    unselectedJets.append(jet)
                    continue
                    
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

