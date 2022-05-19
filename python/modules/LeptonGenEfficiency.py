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


# Calculate the gen-level displacement between two vertices
def getDisplacement(p1, p2):
    x1 = p1.vertex_x
    x2 = p2.vertex_x
    y1 = p1.vertex_y
    y2 = p2.vertex_y
    z1 = p1.vertex_z
    z2 = p2.vertex_z

    Lxy = math.sqrt((x1-x2) ** 2 + (y1-y2) ** 2)
    Lz = math.sqrt((z1-z2) ** 2)
    return Lxy, Lz

# Returned truth-matched reco-level leptons with a given deltaR threshold
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

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        self.features = ["pt", "eta", "phi", "RecoPt", "RecoEta", "RecoPhi", "RecoDeltaR", "isRecoMatched", "Lz", "Lxy"]
        self.electron_features = self.features + ["mvaFall17V2noIso_WPL", "customID"]
        self.muon_features  = self.features + ["looseId", "mediumId", "tightId"]
        self.jet_features = ["pt", "eta", "phi", "is_matched", "dxy", "dz", "reco_pt", "reco_eta", "reco_phi", "tightId", "cpf_sum", "npf_sum"]
        self.hnl_features = ["pt", "eta", "phi", "mass", "boost", "Lxy"]
        self.tau_features = self.features

        self.out = wrappedOutputTree
        self.out.branch("nGenMuons", "I")
        self.out.branch("nGenElectrons", "I")
        self.out.branch("nGenHadTaus", "I")
        self.out.branch("nLeptons", "I")
        self.out.branch("nGenJets", "I")
        self.out.branch("nGenLeptonsFromV", "I")
        self.out.branch("nGenLeptonsFromTauFromV", "I")

        for feature in self.muon_features:
            self.out.branch("GenMuon_"+feature, "F", lenVar="nGenMuons")
        for feature in self.electron_features:
            self.out.branch("GenElectron_"+feature, "F", lenVar="nGenElectrons")
        for feature in self.tau_features:
            self.out.branch("GenHadTau_"+feature, "F", lenVar="nGenHadTaus")
        for feature in self.features:
            self.out.branch("GenLepton_"+feature, "F", lenVar="nGenLeptons")
        for feature in self.jet_features:
             self.out.branch("GenJet_"+feature, "F", lenVar="nGenJets")
        for feature in self.hnl_features:
            self.out.branch("HNL_"+feature, "F")
    
        self.out.branch("GenLeptonsFromV_pt", "F", lenVar="nGenLeptonsFromV")
        self.out.branch("GenLeptonsFromV_eta", "F", lenVar="nGenLeptonsFromV")
        self.out.branch("GenLeptonsFromV_pdgId", "F", lenVar="nGenLeptonsFromV")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        # Gen-level 
        genParticles = self.genInputCollection(event)
        genJets = Collection(event, "GenJet")

        # Reco-level
        electrons = self.electronCollection(event)
        muons = self.muonCollection(event)
        jets = self.jetCollection(event)

        # Extra reco-level
        cpfs = Collection(event, "cpf")
        npfs = Collection(event, "npf")

        genLeptonsFromV = []
        HNLHadronicDaughters = []
        hnl = None

        for genParticle in genParticles:
            # Check if has a mother, otherwise skip
            if genParticle.genPartIdxMother == -1:
                continue
            if abs(genParticle.pdgId) in [9900012, 9990012]:
                continue
            # Check if from (anti)HNL decay
            if abs(genParticles[genParticle.genPartIdxMother].pdgId) in [9900012, 9990012]:
                hnl = genParticles[genParticle.genPartIdxMother]
                hnl_daughter = genParticle
                # Check if not muon, electron, tau
                if abs(hnl_daughter.pdgId) not in [11, 12, 13, 14, 15, 16]:
                    HNLHadronicDaughters.append(hnl_daughter)
        if hnl is None:
            print("No HNL found -- skipping event!")
            return False

        # Calculate some features
        hnl.boost = hnl.pt/hnl.mass
        hnl.Lxy, hnl.Lz = getDisplacement(hnl, hnl_daughter)

        for feature in self.hnl_features:
            self.out.fillBranch("HNL_"+feature, getattr(hnl, feature))

        for genParticle in genParticles:
            # Proceed if particle has a mother and the particle is a muon or electron
            if genParticle.genPartIdxMother != -1 and abs(genParticle.pdgId) in [11, 13]:
                mother = genParticles[genParticle.genPartIdxMother]
                # W or Z boson
                if abs(mother.pdgId) in [23, 24]:
                    genLeptonsFromV.append(genParticle)
                # Special case: W/Z -> tau -> tau_l
                elif abs(mother.pdgId) == 15:
                    if mother.genPartIdxMother != -1:
                        grandmother = genParticles[mother.genPartIdxMother]
                        if abs(grandmother.pdgId) in [23, 24]:
                            genLeptonsFromV.append(genParticle)

        # Need to sort gen leptons
        genLeptonsFromV = sorted(genLeptonsFromV, key=lambda x: x.pt, reverse=True)

        self.out.fillBranch("nGenLeptonsFromV", len(genLeptonsFromV)) 
        self.out.fillBranch("GenLeptonsFromV_pt", [getattr(lepton, "pt") for lepton in genLeptonsFromV])
        self.out.fillBranch("GenLeptonsFromV_eta", [getattr(lepton, "eta") for lepton in genLeptonsFromV])
        self.out.fillBranch("GenLeptonsFromV_pdgId", [getattr(lepton, "pdgId") for lepton in genLeptonsFromV])

        genElectrons = list(filter(lambda genParticle: abs(genParticle.pdgId) == 11 and genParticle.genPartIdxMother != -1, genParticles))
        genMuons = list(filter(lambda genParticle: abs(genParticle.pdgId) == 13 and genParticle.genPartIdxMother != -1, genParticles))
        genTaus = list(filter(lambda genParticle: abs(genParticle.pdgId) == 15 and genParticle.genPartIdxMother != -1, genParticles))

        electronsFromHNL = []
        muonsFromHNL = []
        hadTausFromHNL = []

        for lepton in genElectrons+genMuons+genTaus:
            # Check if parent is HNL and categorise
            if abs(genParticles[lepton.genPartIdxMother].pdgId) in [9900012, 9990012]:
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
                        hadTausFromHNL.append(lepton)

        electronsFromHNL = sorted(electronsFromHNL, key=lambda x: x.pt, reverse=True)
        muonsFromHNL = sorted(muonsFromHNL, key=lambda x: x.pt, reverse=True)
        hadTausFromHNL = sorted(hadTausFromHNL, key=lambda x: x.pt, reverse=True)

        self.out.fillBranch("nGenElectrons", len(electronsFromHNL)) 
        self.out.fillBranch("nGenMuons", len(muonsFromHNL))
        self.out.fillBranch("nGenHadTaus", len(hadTausFromHNL))

        # Match to reco leptons
        # Features common for all leptons
        for lepton in electronsFromHNL+muonsFromHNL+hadTausFromHNL:
            if abs(lepton.pdgId) == 11:
                matchedLepton = matchLepton(lepton, electrons)        
            elif abs(lepton.pdgId) == 13:
                matchedLepton = matchLepton(lepton, muons)   
            elif abs(lepton.pdgId) == 15:
                matchedLepton = matchLepton(lepton, jets)     

            # Calculate some properties
            Lxy, Lz = getDisplacement(lepton, genParticles[lepton.genPartIdxMother])
            lepton.Lz = Lz
            lepton.Lxy = Lxy

            lepton.RecoPt = -999.
            lepton.RecoEta = -999.
            lepton.RecoPhi = -999.
            lepton.RecoDeltaR = -999.
            lepton.isRecoMatched = 0

            if matchedLepton:
                lepton.RecoPt = matchedLepton.pt
                lepton.RecoEta = matchedLepton.eta
                lepton.RecoPhi = matchedLepton.phi
                lepton.RecoDeltaR = deltaR(matchedLepton, lepton)
                lepton.isRecoMatched = 1
                lepton.matchedLepton = matchedLepton

        # Lepton-specific extra features
        for electron in electronsFromHNL:
            # Check if matched reco-electron passes ID
            if electron.isRecoMatched:
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

                electron.customID = electron.GsfEleSCEtaMultiRangeCut>1 and \
                                    electron.GsfEleDEtaInSeedCut>1 and \
                                    electron.GsfEleDPhiInCut>1 and \
                                    electron.GsfEleFull5x5SigmaIEtaIEtaCut>1 and \
                                    electron.GsfEleEInverseMinusPInverseCut>1 and \
                                    electron.GsfEleConversionVetoCut>1
            else:
                electron.mvaFall17V2noIso_WPL = -1.
                electron.customID = -1.

        for muon in muonsFromHNL:
            # Check if matched reco-muon passes ID

            if muon.isRecoMatched:
                muon.looseId = muon.matchedLepton.looseId
                muon.mediumId = muon.matchedLepton.mediumId
                muon.tightId = muon.matchedLepton.tightId
            else:
                muon.looseId = -1.
                muon.mediumId = -1.
                muon.tightId = -1.

        leptonsFromHNL = electronsFromHNL+muonsFromHNL
        leptonsFromHNL = sorted(leptonsFromHNL, key=lambda x: x.pt, reverse=True)

        self.out.fillBranch("nGenLeptons", len(leptonsFromHNL))    

        for feature in self.electron_features:
            self.out.fillBranch("GenElectron_"+feature, [getattr(lepton, feature) for lepton in electronsFromHNL])
        for feature in self.muon_features:
            self.out.fillBranch("GenMuon_"+feature, [getattr(lepton, feature) for lepton in muonsFromHNL])
        for feature in self.features:
            self.out.fillBranch("GenLepton_"+feature, [getattr(lepton, feature) for lepton in leptonsFromHNL])
        for feature in self.tau_features:
            self.out.fillBranch("GenHadTau_"+feature, [getattr(lepton, feature) for lepton in hadTausFromHNL])

        selectedgenjets = []

        for genjet in genJets:
            matched_jet = None
            genjet.is_matched = False
            genjet.reco_pt = -1.
            genjet.reco_eta = -1.
            genjet.reco_phi = -1.
            genjet.tightId = -1

            for jet in jets:
                if deltaR(jet, genjet) < 0.4:
                    genjet.is_matched = True
                    matched_jet = jet
                    genjet.reco_pt = jet.pt
                    genjet.reco_eta = jet.eta
                    genjet.reco_phi = jet.phi
                    genjet.tightId = jet.jetId > 1
                    break
            genjet.cpf_sum = -1
            genjet.npf_sum = -1
            for quark in HNLHadronicDaughters:
                if deltaR(quark, genjet) < 0.4:
                    if matched_jet:
                        cpf_sum = 0
                        npf_sum = 0
                        for cpf in cpfs:
                            if cpf.jetIdx == matched_jet._index:
                                cpf_sum += cpf.ptrel
                        for npf in npfs:
                            if npf.jetIdx == matched_jet._index:
                                npf_sum += npf.ptrel  
                        genjet.cpf_sum = cpf_sum 
                        genjet.npf_sum = npf_sum
                    selectedgenjets.append(genjet)
                    dxy, dz = getDisplacement(hnl, quark)
                    genjet.dxy = dxy
                    genjet.dz = dz
                    break
    
        self.out.fillBranch("nGenJets", len(selectedgenjets))   

        for feature in self.jet_features:
            self.out.fillBranch("GenJet_"+feature, [getattr(jet, feature) for jet in selectedgenjets]) 

        return True
