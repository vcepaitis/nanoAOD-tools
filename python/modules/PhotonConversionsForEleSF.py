import os
import sys
import math
import json
import ROOT
import random

from utils import deltaR, deltaPhi, splitNameType


from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

cpfFeatures = [
    "trackEtaRel",
    "trackPtRel",
    "trackPPar",
    "trackDeltaR",
    "trackPParRatio",
    "trackPtRatio",
    "trackSip2dVal",
    "trackSip2dSig",
    "trackSip3dVal",
    "trackSip3dSig",
    "trackJetDistVal",
    "trackJetDistSig",
    "drminsv",
    "vertex_association",
    "fromPV",
    "track_chi2",
    "track_quality",
    "track_numberOfValidPixelHits",
    "track_pixelLayersWithMeasurement",
    "track_numberOfValidStripHits",
    "track_stripLayersWithMeasurement", 
    "matchedMuon",
    "matchedElectron",
    "matchedSV",
    "matchedSV_adapted",
    "track_ndof"
]

class PhotonConversionsForEleSF(Module):

    def __init__(
        self,
        tightMuonsCollection = lambda event: Collection(event, "Muon"),
        tightElectronsCollection = lambda event: Collection(event, "Electron"),
        looseMuonCollection = lambda event: Collection(event, "Muon"),
        looseElectronCollection = lambda event: Collection(event, "Electron"),
        outputName = "Leptons",
        globalOptions={"isData": False, "year": 2016},
        storeLeadingKinematics=["pt", "eta", "phi", "charge", "isMuon", "isElectron", "relIso", "dxy", "dz", 'dxysig', 'dzsig', 'isTriggerMatched'],
        storeSubleadingKinematics=["pt", "eta", "phi", "charge", "isMuon", "isElectron", "relIso", "dxy", "dz", 'dxysig', 'dzsig', 'isTriggerMatched'],
        storeTrailingKinematics=["pt", "eta", "phi", "charge/I", "isMuon/I", "isElectron/I", "relIso", "dxy", "dz", 'dxysig', 'dzsig', 'looseId', 'tightId', 'isCustomID', 'isCustomNoConvID']
    ):

        self.globalOptions = globalOptions
        self.tightMuonsCollection = tightMuonsCollection
        self.tightElectronsCollection = tightElectronsCollection
        self.looseMuonCollection = looseMuonCollection
        self.looseElectronCollection = looseElectronCollection
        self.outputName = outputName
        self.storeLeadingKinematics = storeLeadingKinematics
        self.storeSubleadingKinematics = storeSubleadingKinematics
        self.storeTrailingKinematics = storeTrailingKinematics

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nleading"+self.outputName, "I")
        self.out.branch("nsubleading"+self.outputName, "I")
        self.out.branch("ntrailing"+self.outputName, "I")

        self.out.branch(self.outputName+"_muonmuonelectron", "I")
        self.out.branch(self.outputName+"_electronelectronelectron", "I")

        self.out.branch("trilepton_mass", "F")
        self.out.branch("trilepton_pt", "F")
        self.out.branch("trilepton_eta", "F")
        self.out.branch("trilepton_charge", "F")
        self.out.branch(self.outputName+"_dPhi_l1l3", "F")
        self.out.branch(self.outputName+"_dEta_l1l3", "F")
        self.out.branch(self.outputName+"_dR_l1l3", "F")
        self.out.branch(self.outputName+"_dPhi_l2l3", "F")
        self.out.branch(self.outputName+"_dEta_l2l3", "F")
        self.out.branch(self.outputName+"_dR_l2l3", "F")

        
        for i,variable in enumerate(self.storeLeadingKinematics):
            name,dtype = splitNameType(variable)
            self.storeLeadingKinematics[i] = name
            self.out.branch("leading"+self.outputName+"_"+name, "F")
        
        for i,variable in enumerate(self.storeSubleadingKinematics):
            name,dtype = splitNameType(variable)
            self.storeSubleadingKinematics[i] = name
            self.out.branch("subleading"+self.outputName+"_"+name, "F")

        for i,variable in enumerate(self.storeTrailingKinematics):
            name,dtype = splitNameType(variable)
            self.storeTrailingKinematics[i] = name
            self.out.branch("trailing"+self.outputName+"_"+name, "F")
            
        '''
        self.out.branch("subleading"+self.outputName+"_cpfMatch", "I", lenVar="nsubleading"+self.outputName)
        for cpfFeature in cpfFeatures:
            self.out.branch("subleading"+self.outputName+"_"+cpfFeature, "F", lenVar="nsubleading"+self.outputName)
        '''
        #for variable in self.storeKinematics:
            #self.out.branch(self. outputName+"_"+variable,"F",lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


        
    def triggerMatched(self, lepton, trigger_objects):
        if lepton.isElectron:
            trigger_objects = filter(lambda obj: abs(obj.id) == 11, trigger_objects)
        elif lepton.isMuon:
            trigger_objects = filter(lambda obj: abs(obj.id) == 13, trigger_objects)

        if len(trigger_objects) == 0:
            return False

        min_delta_R = min(map(lambda obj: deltaR(lepton, obj), trigger_objects))

        if lepton.isElectron:
            if min_delta_R < 0.3:
                return True
            else:
                return False

        elif lepton.isMuon:
            if min_delta_R < 0.1:
                return True
            else:
                return False


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
    
        triggerObjects = Collection(event, "TrigObj")

        tightMuons = self.tightMuonsCollection(event)
        tightElectrons = self.tightElectronsCollection(event)

        looseMuons = self.looseMuonCollection(event)
        looseElectrons = self.looseElectronCollection(event)
        
        cpfCandidates = Collection(event,"cpf")


        for lepton in tightMuons+looseMuons:
            lepton.isMuon = 1
            lepton.isElectron = 0
            lepton.relIso = lepton.pfRelIso04_all
            lepton.cpf_match = None
            lepton.cpf_match_dR = 100.
            for cpfCandidate in cpfCandidates:
                if cpfCandidate.matchedMuon>0.5:
                    pt = math.sqrt(cpfCandidate.px**2+cpfCandidate.py**2)
                    p = math.sqrt(cpfCandidate.px**2+cpfCandidate.py**2+cpfCandidate.pz**2)
                    eta = math.atanh(cpfCandidate.pz/p)
                    phi = math.atan2(cpfCandidate.py,cpfCandidate.px)
                    dR = math.sqrt((eta-lepton.eta)**2+deltaPhi(phi,lepton.phi)**2)
                    if math.fabs(pt/lepton.pt-1)<0.1 and dR<0.02 and dR<lepton.cpf_match_dR:
                        lepton.cpf_match = cpfCandidate
                        lepton.cpf_match_dR = dR

        for lepton in tightElectrons+looseElectrons:
            lepton.isMuon = 0
            lepton.isElectron = 1
            lepton.relIso = lepton.pfRelIso03_all
            lepton.cpf_match = None
            lepton.cpf_match_dR = 100.
            
            for cpfCandidate in cpfCandidates:
                if cpfCandidate.matchedElectron>0.5:
                    pt = math.sqrt(cpfCandidate.px**2+cpfCandidate.py**2)
                    p = math.sqrt(cpfCandidate.px**2+cpfCandidate.py**2+cpfCandidate.pz**2)
                    eta = math.atanh(cpfCandidate.pz/p)
                    phi = math.atan2(cpfCandidate.py,cpfCandidate.px)
                    dR = math.sqrt((eta-lepton.eta)**2+deltaPhi(phi,lepton.phi)**2)
                    if math.fabs(pt/lepton.pt-1)<0.1 and dR<0.02 and dR<lepton.cpf_match_dR:
                        lepton.cpf_match = cpfCandidate
                        lepton.cpf_match_dR = dR

        tightLeptons = []
        looseLeptons = []

        tightLeptons = tightMuons+tightElectrons
        looseLeptons = looseMuons+looseElectrons

        tightLeptons = sorted(tightLeptons, key=lambda x: x.pt, reverse=True)

        if len(tightLeptons)<2:
            print("Need at least two tight leptons")
            return False

        # if len(tightLeptons+looseLeptons)<3:
        #     print("Need at least three leptons")
        #     return False

        # select leading only, move subleading to "loose"
        looseLeptons.extend(tightLeptons[2:])
        tightLeptons = tightLeptons[:2]
        if len(looseLeptons)>0:
            looseLeptons = sorted(looseLeptons, key=lambda x: x.pt, reverse=True)

        muonmuonelectron = 0
        electronelectronelectron = 0

        ## flavour categorisation :
        if len(tightLeptons) > 0 and len(looseLeptons) > 0:
            # Ensure pt(l1) > pt(l3) and pt(l2) > pt(l3)
            if (tightLeptons[0].pt < looseLeptons[0].pt) or (tightLeptons[1].pt < looseLeptons[0].pt):
                return False
                
            if tightLeptons[0].isMuon and tightLeptons[1].isMuon and looseLeptons[0].isElectron:
                muonmuonelectron = 1

            elif tightLeptons[0].isElectron and tightLeptons[1].isElectron and looseLeptons[0].isElectron:
                electronelectronelectron = 1

        for lepton in tightLeptons+looseLeptons:
            if self.triggerMatched(lepton, triggerObjects):
                setattr(lepton,"isTriggerMatched",1)
            else:
                setattr(lepton,"isTriggerMatched",0)

        for lepton in tightLeptons+looseLeptons:
            if lepton.dxyErr < 1e-6:
                 setattr(lepton, "dxysig", -1.)
            else:
                 setattr(lepton, "dxysig", math.fabs(lepton.dxy)/math.fabs(lepton.dxyErr))

            if lepton.dzErr < 1e-6:
                 setattr(lepton, "dzsig", -1.)
            else:
                 setattr(lepton, "dzsig", math.fabs(lepton.dz)/math.fabs(lepton.dzErr))

        # self.out.fillBranch("nleading"+self.outputName, int(float(len(tightLeptons))/2.))
        # self.out.fillBranch("nsubleading"+self.outputName, int(float(len(tightLeptons))/2.))
        self.out.fillBranch("nleading"+self.outputName, 1)
        self.out.fillBranch("nsubleading"+self.outputName, 1)
        self.out.fillBranch("ntrailing"+self.outputName, len(looseLeptons))

        for variable in self.storeLeadingKinematics:
            #self.out.fillBranch("leading"+self.outputName+"_"+variable,map(lambda lepton: getattr(lepton,variable),tightLeptons))
            self.out.fillBranch("leading"+self.outputName+"_"+variable,getattr(tightLeptons[0], variable))
        for variable in self.storeSubleadingKinematics:
            #self.out.fillBranch("subleading"+self.outputName+"_"+variable,map(lambda lepton: getattr(lepton,variable),tightLeptons))
            self.out.fillBranch("subleading"+self.outputName+"_"+variable,getattr(tightLeptons[1], variable))
        if len(looseLeptons)>0 and looseLeptons[0].isElectron:
            for variable in self.storeTrailingKinematics:
                #self.out.fillBranch("trailing"+self.outputName+"_"+variable,map(lambda lepton: getattr(lepton,variable),looseLeptons))
                self.out.fillBranch("trailing"+self.outputName+"_"+variable,getattr(looseLeptons[0], variable))
        else:
            for variable in self.storeTrailingKinematics:
                self.out.fillBranch("trailing"+self.outputName+"_"+variable, -99)


        self.out.fillBranch(self.outputName+"_muonmuonelectron", muonmuonelectron)
        self.out.fillBranch(self.outputName+"_electronelectronelectron", electronelectronelectron)

        if len(looseLeptons)>0 and looseLeptons[0].isElectron:
            trilepton = ROOT.TLorentzVector()
            trilepton = tightLeptons[0].p4()+tightLeptons[1].p4()+looseLeptons[0].p4()
            trilepton_charge = tightLeptons[0].charge*tightLeptons[1].charge*looseLeptons[0].charge

            self.out.fillBranch("trilepton_mass", trilepton.M())
            self.out.fillBranch("trilepton_pt", trilepton.Pt())
            self.out.fillBranch("trilepton_eta", trilepton.Eta())
            self.out.fillBranch("trilepton_charge", trilepton_charge)
            self.out.fillBranch(self.outputName+"_dPhi_l1l3", math.fabs(deltaPhi(tightLeptons[0],looseLeptons[0])))
            self.out.fillBranch(self.outputName+"_dEta_l1l3", math.fabs(tightLeptons[0].eta-looseLeptons[0].eta))
            self.out.fillBranch(self.outputName+"_dR_l1l3", deltaR(tightLeptons[0],looseLeptons[0]))
            self.out.fillBranch(self.outputName+"_dPhi_l2l3", math.fabs(deltaPhi(tightLeptons[1],looseLeptons[0])))
            self.out.fillBranch(self.outputName+"_dEta_l2l3", math.fabs(tightLeptons[1].eta-looseLeptons[0].eta))
            self.out.fillBranch(self.outputName+"_dR_l2l3", deltaR(tightLeptons[1],looseLeptons[0]))
        else:
            self.out.fillBranch("trilepton_mass", -99)
            self.out.fillBranch("trilepton_pt", -99)
            self.out.fillBranch("trilepton_eta", -99)
            self.out.fillBranch("trilepton_charge", -99)
            self.out.fillBranch(self.outputName+"_dPhi_l1l3", -99)
            self.out.fillBranch(self.outputName+"_dEta_l1l3", -99)
            self.out.fillBranch(self.outputName+"_dR_l1l3", -99)
            self.out.fillBranch(self.outputName+"_dPhi_l2l3", -99)
            self.out.fillBranch(self.outputName+"_dEta_l2l3", -99)
            self.out.fillBranch(self.outputName+"_dR_l2l3", -99)

        setattr(event, "leading"+self.outputName, tightLeptons[0])
        setattr(event, "subleading"+self.outputName, tightLeptons[1])
        if len(looseLeptons)>0 and looseLeptons[0].isElectron:
            setattr(event, "trailing"+self.outputName, looseLeptons[0])
        else:
            setattr(event, "trailing"+self.outputName, None)

        return True
