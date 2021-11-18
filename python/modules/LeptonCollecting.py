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

class LeptonCollecting(Module):

    def __init__(
        self,
        tightMuonsCollection = lambda event: Collection(event, "Muon"),

        tightElectronsCollection = lambda event: Collection(event, "Electron"),

        looseMuonCollection = lambda event: Collection(event, "Muon"),
        looseElectronCollection = lambda event: Collection(event, "Electron"),

        outputName = "Leptons",
        globalOptions={"isData": False, "year": 2016},
        storeLeadingKinematics=["pt", "eta", "phi", "charge", "isMuon", "isElectron", "relIso", "dxy", "dz", 'dxysig', 'dzsig', 'isTriggerMatched'],
        storeSubleadingKinematics=["pt", "eta", "phi", "charge", "isMuon", "isElectron", "relIso", "dxy", "dz", 'dxysig', 'dzsig', 'looseId', 'tightId']
    ):

        self.globalOptions = globalOptions
        self.tightMuonsCollection = tightMuonsCollection
        self.tightElectronsCollection = tightElectronsCollection
        self.looseMuonCollection = looseMuonCollection
        self.looseElectronCollection = looseElectronCollection
        self.outputName = outputName
        self.storeLeadingKinematics = storeLeadingKinematics
        self.storeSubleadingKinematics = storeSubleadingKinematics

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nleading"+self.outputName, "I")
        self.out.branch("nsubleading"+self.outputName, "I")

        self.out.branch(self.outputName+"_muonmuon", "I")
        self.out.branch(self.outputName+"_electronelectron", "I")
        self.out.branch(self.outputName+"_muonelectron", "I")
        self.out.branch(self.outputName+"_electronmuon", "I")
        self.out.branch(self.outputName+"_muonjets", "I")
        self.out.branch(self.outputName+"_electronjets", "I")
        self.out.branch(self.outputName+"_subLeptonTightId", "I")
        self.out.branch(self.outputName+"_subLeptonChi2", "F")
        
        for i,variable in enumerate(self.storeLeadingKinematics):
            name,dtype = splitNameType(variable)
            self.storeLeadingKinematics[i] = name
            self.out.branch("leading"+self.outputName+"_"+name, dtype, lenVar="nleading"+self.outputName)
        
        for i,variable in enumerate(self.storeSubleadingKinematics):
            name,dtype = splitNameType(variable)
            self.storeSubleadingKinematics[i] = name
            self.out.branch("subleading"+self.outputName+"_"+name, dtype, lenVar="nsubleading"+self.outputName)
            
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
        muonCandidates = Collection(event,"muon")

        electronCandidates = Collection(event,"electron")



        for lepton in tightMuons+looseMuons:
            lepton.isMuon = 1
            lepton.isElectron = 0
            lepton.chi2 = -100.
            lepton.muon_match = None
            lepton.muon_match_dR = 100.
            for muonCandidate in muonCandidates:
               pt  = math.sqrt(muonCandidate.px**2+muonCandidate.py**2)
               p   = math.sqrt(muonCandidate.px**2+muonCandidate.py**2+muonCandidate.pz**2)
               eta = math.atanh(muonCandidate.pz/p)
               phi = math.atan2(muonCandidate.py,muonCandidate.px)
	       dphi =math.fabs( phi - lepton.phi)
               dpt = pt - lepton.pt
               dR =  math.sqrt((eta-lepton.eta)**2+ (phi - lepton.phi)**2)
               if  dphi < 0.1 and dR<0.1 :
               #if  dR<0.1 :
                        lepton.muon_match = muonCandidate
                        lepton.muon_match_dR = dR
			lepton.chi2 = muonCandidate.chi2 

        for lepton in tightElectrons+looseElectrons:
            lepton.isMuon = 0
            lepton.isElectron = 1
            lepton.relIso = lepton.pfRelIso03_all
            lepton.chi2 = -100.
            lepton.electron_match_dR = 100.
            lepton.electron_match = None
            for electronCandidate in electronCandidates:

               pt  = math.sqrt(electronCandidate.px**2+electronCandidate.py**2)
               p   = math.sqrt(electronCandidate.px**2+electronCandidate.py**2+electronCandidate.pz**2)
               eta = math.atanh(electronCandidate.pz/p)
               phi = math.atan2(electronCandidate.py,electronCandidate.px)

               dpt = pt - lepton.pt
               dphi =math.fabs( phi - lepton.phi)
               dR =  math.sqrt((eta-lepton.eta)**2 +(phi - lepton.phi)**2)
               if dR< 0.1 and dphi < 0.1:
               #if  dR<0.1 :
                        lepton.electron_match = electronCandidate
                        lepton.electron_match_dR = dR
			lepton.chi2 = electronCandidate.chi2
    
	########

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
        
        # select leading only, move subleading to "loose"
        looseLeptons.extend(tightLeptons[1:])
        tightLeptons = tightLeptons[:1]
        looseLeptons = sorted(looseLeptons, key=lambda x: x.pt, reverse=True)

        muonmuon = 0
        electronelectron = 0
        muonelectron = 0
        electronmuon = 0
        muonjets = 0
        electronjets = 0

        ## flavour categorisation :
        if len(tightLeptons) > 0 and len(looseLeptons) > 0:
            # Ensure pt(l1) > pt(l2)
            if tightLeptons[0].pt < looseLeptons[0].pt:
                return False
                
            if tightLeptons[0].isMuon and looseLeptons[0].isMuon:
                muonmuon = 1

            elif tightLeptons[0].isElectron and looseLeptons[0].isElectron:
                electronelectron= 1

            elif tightLeptons[0].isMuon and looseLeptons[0].isElectron:
                muonelectron = 1
               
            elif tightLeptons[0].isElectron and looseLeptons[0].isMuon:
                electronmuon = 1

        elif len(tightLeptons) > 0:
            if tightLeptons[0].isMuon:
                muonjets = 1
            elif tightLeptons[0].isElectron:
                electronjets = 1

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

        lepton_id = -1
        lepton_chi2 = -1. 
        if len(looseLeptons) > 0:
          if looseLeptons[0].isMuon: 
		lepton_id = looseLeptons[0].tightId
		lepton_chi2 = looseLeptons[0].chi2
          if looseLeptons[0].isElectron: 
                lepton_id = looseLeptons[0].mvaFall17V1Iso_WP90
		lepton_chi2 = looseLeptons[0].chi2

        self.out.fillBranch("nleading"+self.outputName, len(tightLeptons))
        self.out.fillBranch("nsubleading"+self.outputName, len(looseLeptons))

        for variable in self.storeLeadingKinematics:
            self.out.fillBranch("leading"+self.outputName+"_"+variable,map(lambda lepton: getattr(lepton,variable),tightLeptons))
        for variable in self.storeSubleadingKinematics:
            self.out.fillBranch("subleading"+self.outputName+"_"+variable,map(lambda lepton: getattr(lepton,variable),looseLeptons))

        ''' 
        self.out.fillBranch("subleading"+self.outputName+"_cpfMatch", [1 if lepton.cpf_match!=None else 0 for lepton in looseLeptons])
        for cpfFeature in cpfFeatures:
            arr = [getattr(lepton.cpf_match,cpfFeature) if lepton.cpf_match!=None else 0 for lepton in  looseLeptons]
            self.out.fillBranch("subleading"+self.outputName+"_"+cpfFeature,arr)
        '''
                       
          
        self.out.fillBranch(self.outputName+"_muonmuon", muonmuon)
        self.out.fillBranch(self.outputName+"_electronelectron", electronelectron)
        self.out.fillBranch(self.outputName+"_muonelectron", muonelectron)
        self.out.fillBranch(self.outputName+"_electronmuon", electronmuon)
        self.out.fillBranch(self.outputName+"_muonjets", muonjets)
        self.out.fillBranch(self.outputName+"_electronjets", electronjets)
        self.out.fillBranch(self.outputName+"_subLeptonTightId", lepton_id)
        self.out.fillBranch(self.outputName+"_subLeptonChi2", lepton_chi2)
   

        setattr(event, "leading"+self.outputName, tightLeptons)
        setattr(event, "subleading"+self.outputName, looseLeptons)

        return True
