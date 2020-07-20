import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from utils import deltaPhi, deltaR

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

def getDisplacement(p1, p2):
    x1 = p1.vertex_x
    x2 = p2.vertex_x
    y1 = p1.vertex_y
    y2 = p2.vertex_y
    z1 = p1.vertex_z
    z2 = p2.vertex_z

    dy = math.sqrt((x1-x2) ** 2 + (y1-y2) ** 2)
    d3 = math.sqrt((x1-x2) ** 2 + (y1-y2) ** 2 + (z1-z2) ** 2)
    return d3, dy

def matchLepton(genLepton, recoLeptons, minDeltaRThreshold=0.4):
    if len(recoLeptons) > 0:
        minDeltaR = 1000.
        matchedLepton = recoLeptons[0]
        for recoLepton in recoLeptons:
            if deltaR(recoLepton, genLepton) < minDeltaR:
                matchedLepton = recoLepton
                minDeltaR = deltaR(recoLepton, genLepton) 
        if minDeltaR < minDeltaRThreshold:
            return matchedLepton
        else:
            return None
    else:
        return None

class LeptonGenEfficiency(Module):
    def __init__(
        self,
        genInputCollection=lambda event: Collection(event, "GenPart"),
        electronCollection=lambda event: Collection(event, "Electron"),
        muonCollection=lambda event: Collection(event, "Muon"),
        jetCollection=lambda event: Collection(event, "Jet"),
        globalOptions={"isData": False}
    ):
        self.genInputCollection = genInputCollection
        self.electronCollection = electronCollection
        self.muonCollection = muonCollection
        self.jetCollection = jetCollection
        self.globalOptions = globalOptions
        features = ["pt", "eta", "phi", "reco_pt", "reco_eta", "reco_phi", "reco_deltaR", "is_matched", "d3", "dxy"]
        self.electron_features = features + ["mvaFall17V2noIso_WPL", "customID"]
        self.muon_features = features + ["looseId", "mediumId", "tightId"]
        self.tau_features = features + ["looseId"]

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nmuon", "I")
        self.out.branch("nelectron", "I")
        self.out.branch("ntau", "I")

        for feature in self.muon_features:
            self.out.branch("muon_"+feature, "F", lenVar="nmuon")
        for feature in self.electron_features:
            self.out.branch("electron_"+feature, "F", lenVar="nelectron")
        for feature in self.tau_features:
            self.out.branch("tau_"+feature, "F", lenVar="ntau")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        genParticles = self.genInputCollection(event)
        electrons = self.electronCollection(event)
        muons = self.muonCollection(event)
        jets = self.jetCollection(event)

        genElectrons = list(filter(lambda genParticle: abs(genParticle.pdgId) == 11 and genParticle.genPartIdxMother != -1, genParticles))
        genMuons = list(filter(lambda genParticle: abs(genParticle.pdgId) == 13 and genParticle.genPartIdxMother != -1, genParticles))
        genTaus = list(filter(lambda genParticle: abs(genParticle.pdgId) == 15 and genParticle.genPartIdxMother != -1, genParticles))

        electronsFromHNL = []
        muonsFromHNL = []
        tausFromHNL = []

        for lepton in genElectrons+genMuons+genTaus:
            if abs(genParticles[lepton.genPartIdxMother].pdgId) == 9990012:
                if abs(lepton.pdgId) == 11: 
                    electronsFromHNL.append(lepton)
                elif abs(lepton.pdgId) == 13: 
                    muonsFromHNL.append(lepton)
                elif abs(lepton.pdgId) == 15:
                    decaysToE = 0
                    decaysToMu = 0
                    for lep in genElectrons+genMuons:
                        if lepton._index == lep.genPartIdxMother:
                            if abs(lep.pdgId) == 11:
                                decaysToE = True
                                electronsFromHNL.append(lep)
                                break
                            elif abs(lep.pdgId) == 13:
                                decaysToMu = True
                                muonsFromHNL.append(lep)
                                break
                    if decaysToE+decaysToMu == 0:
                        tausFromHNL.append(lepton)
                    #print("bitmap: {0:b} ".format(lepton.statusFlags))

        self.out.fillBranch("nelectron", len(electronsFromHNL)) 
        self.out.fillBranch("nmuon", len(muonsFromHNL))
        self.out.fillBranch("ntau", len(tausFromHNL))

        #for i, part in enumerate(genParticles):
            #print i, part.pdgId, part.genPartIdxMother

        for lepton in electronsFromHNL+muonsFromHNL+tausFromHNL:
            if abs(lepton.pdgId) == 11:
                matchedLepton = matchLepton(lepton, electrons)        
            elif abs(lepton.pdgId) == 13:
                matchedLepton = matchLepton(lepton, muons)   
            elif abs(lepton.pdgId) == 15:
                matchedLepton = matchLepton(lepton, jets)     

            d3, dxy = getDisplacement(lepton, genParticles[lepton.genPartIdxMother])
            lepton.d3 = d3
            lepton.dxy = dxy

            if matchedLepton:
                lepton.reco_pt = matchedLepton.pt
                lepton.reco_eta = matchedLepton.eta
                lepton.reco_phi = matchedLepton.phi
                lepton.reco_deltaR = deltaR(matchedLepton, lepton)
                lepton.is_matched = 1
                lepton.matchedLepton = matchedLepton
            else:
                lepton.reco_pt = -999.
                lepton.reco_eta = -999.
                lepton.reco_phi = -999.
                lepton.reco_deltaR = -999.
                lepton.is_matched = 0

        for electron in electronsFromHNL:
            electron.d3 = d3
            electron.dxy = dxy
            if electron.is_matched:
                electron.mvaFall17V2noIso_WPL = electron.matchedLepton.mvaFall17V2noIso_WPL
                bitmap = electron.matchedLepton.vidNestedWPBitmap
                # decision for each cut represented by 3 bits (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
                # Electron_vidNestedWPBitmap
                cuts = np.empty(10)

                for i in range(10):
                    cuts[i] = (bitmap >> i*3) & 0x7

                electron.MinPtCut = cuts[0]
                electron.GsfEleSCEtaMultiRangeCut = cuts[1]
                electron.GsfEleDEtaInSeedCut = cuts[2]
                electron.GsfEleDPhiInCut = cuts[3]
                electron.GsfEleFull5x5SigmaIEtaIEtaCut = cuts[4]
                electron.GsfEleHadronicOverEMEnergyScaledCut = cuts[5]
                electron.GsfEleEInverseMinusPInverseCut = cuts[6]
                electron.GsfEleRelPFIsoScaledCut = cuts[7]
                electron.GsfEleConversionVetoCut = cuts[8]
                electron.GsfEleMissingHitsCut = cuts[9]

                electron.customID = electron.GsfEleSCEtaMultiRangeCut>0 and \
                                    electron.GsfEleDEtaInSeedCut>0 and \
                                    electron.GsfEleDPhiInCut>0 and \
                                    electron.GsfEleFull5x5SigmaIEtaIEtaCut>0 and \
                                    electron.GsfEleEInverseMinusPInverseCut>0 and \
                                    electron.GsfEleConversionVetoCut>0
            else:
                electron.mvaFall17V2noIso_WPL = -1.
                electron.customID = -1.

        for muon in muonsFromHNL:
            muon.d3 = d3
            muon.dxy = dxy
            if muon.is_matched:
                muon.looseId = muon.matchedLepton.looseId
                muon.mediumId = muon.matchedLepton.mediumId
                muon.tightId = muon.matchedLepton.tightId
            else:
                muon.looseId = -1.
                muon.mediumId = -1.
                muon.tightId = -1.

        for tau in tausFromHNL:
            tau.d3 = d3
            tau.dxy = dxy
            if tau.is_matched:
                tau.looseId = tau.matchedLepton.jetId>0
            else:
                tau.looseId = -1.


        for feature in self.electron_features:
            self.out.fillBranch("electron_"+feature, [getattr(lepton, feature) for lepton in electronsFromHNL])
        for feature in self.muon_features:
            self.out.fillBranch("muon_"+feature, [getattr(lepton, feature) for lepton in muonsFromHNL])
        for feature in self.tau_features:
            self.out.fillBranch("tau_"+feature, [getattr(lepton, feature) for lepton in tausFromHNL])

        return True
