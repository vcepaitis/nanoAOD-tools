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
        leadingLeptons=lambda event: [],
        subleadingLeptons=lambda event: [],
        metInput=lambda event: Object(event, "MET"),
        outputName="EventObservables",
        globalOptions={"isData": False}
    ):
        self.jetCollection = jetCollection
        self.leadingLeptons = leadingLeptons
        self.subleadingLeptons = subleadingLeptons
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
        
        self.out.branch(self.outputName+"_mJetsLepton", "F")
        self.out.branch(self.outputName+"_mMaxJetLepton", "F")
        self.out.branch(self.outputName+"_jetMaxMin_ptR", "F") #ratio between highest and lowest jet pt
        self.out.branch(self.outputName+"_htSJet_ptRMax", "F") #max ratio between ht and any jet pt 
        self.out.branch(self.outputName+"_htSJet_ptRMin", "F") #min ratio between ht and any jet pt 
        self.out.branch(self.outputName+"_htVecJet_ptRMax", "F") #max ratio between ht and any jet pt 
        self.out.branch(self.outputName+"_htVecJet_ptRMin", "F") #min ratio between ht and any jet pt 
        
        self.out.branch(self.outputName+"_htVecLepton_dphi", "F")
        self.out.branch(self.outputName+"_htVecLepton_ptR", "F")
        self.out.branch(self.outputName+"_htSLepton_ptR", "F")
        
        self.out.branch(self.outputName+"_minJetSubLepton_dphi", "F")
        self.out.branch(self.outputName+"_minJetSubLepton_dR", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        jets = self.jetCollection(event)
        leadingLepton = self.leadingLeptons(event)[0]
        looseLeptons = self.subleadingLeptons(event)

        met = self.metInput(event)
        self.out.fillBranch(self.outputName+"_met", met.pt)
        self.out.fillBranch(self.outputName+"_met_phi", met.phi)

        # ht and mht
        vectorSum = ROOT.TLorentzVector()
        scalarPtSum = 0.0

        for jet in jets:
            vectorSum += jet.p4()
            scalarPtSum += jet.pt
            
        mht_met_dphi = math.fabs(deltaPhi(vectorSum.Phi(),met.phi))

        # minPhiStar
        minPhiStar = math.pi
        for jet in jets:
            negSum = -(vectorSum-jet.p4())
            minPhiStar = min(minPhiStar, math.fabs(deltaPhi(negSum.Phi(), jet.phi)))
            
        if len(jets)>0:
            self.out.fillBranch(self.outputName+"_mJetsLepton", (vectorSum+leadingLepton.p4()).M())
            self.out.fillBranch(self.outputName+"_mMaxJetLepton", max(map(lambda j: (j.p4()+leadingLepton.p4()).M(),jets)))
            self.out.fillBranch(self.outputName+"_jetMaxMin_ptR", jets[0].pt/jets[-1].pt) #ratio between highest and lowest jet pt
            self.out.fillBranch(self.outputName+"_htSJet_ptRMax", max(map(lambda j: j.pt/scalarPtSum,jets))) #max ratio between ht and any jet pt 
            self.out.fillBranch(self.outputName+"_htSJet_ptRMin", min(map(lambda j: j.pt/scalarPtSum,jets))) #min ratio between ht and any jet pt 
            self.out.fillBranch(self.outputName+"_htVecJet_ptRMax", max(map(lambda j: j.pt/vectorSum.Pt(),jets))) #max ratio between ht and any jet pt 
            self.out.fillBranch(self.outputName+"_htVecJet_ptRMin", min(map(lambda j: j.pt/vectorSum.Pt(),jets))) #min ratio between ht and any jet pt 
            
            self.out.fillBranch(self.outputName+"_htVecLepton_dphi", math.fabs(deltaPhi(vectorSum.Phi(),leadingLepton.phi)))
            self.out.fillBranch(self.outputName+"_htVecLepton_ptR", vectorSum.Pt()/leadingLepton.pt)
            self.out.fillBranch(self.outputName+"_htSLepton_ptR", scalarPtSum/leadingLepton.pt)

            self.out.fillBranch(self.outputName+"_minPhiStar", minPhiStar)
            self.out.fillBranch(self.outputName+"_ht", scalarPtSum)
            self.out.fillBranch(self.outputName+"_hmass", vectorSum.M())
            self.out.fillBranch(self.outputName+"_mht", vectorSum.Pt())
            self.out.fillBranch(self.outputName+"_mht_phi", vectorSum.Phi())
            self.out.fillBranch(self.outputName+"_mht_met_dphi", mht_met_dphi)
        else:
            self.out.fillBranch(self.outputName+"_mJetsLepton", 0.)
            self.out.fillBranch(self.outputName+"_mMaxJetLepton", 0.)
            self.out.fillBranch(self.outputName+"_jetMaxMin_ptR", 0.) #ratio between highest and lowest jet pt
            self.out.fillBranch(self.outputName+"_htSJet_ptRMax", 0.) #max ratio between ht and any jet pt 
            self.out.fillBranch(self.outputName+"_htSJet_ptRMin", 0.) #min ratio between ht and any jet pt 
            self.out.fillBranch(self.outputName+"_htVecJet_ptRMax", 0.) #max ratio between ht and any jet pt 
            self.out.fillBranch(self.outputName+"_htVecJet_ptRMin", 0.) #min ratio between ht and any jet pt 
            
            self.out.fillBranch(self.outputName+"_htVecLepton_dphi", 0.)
            self.out.fillBranch(self.outputName+"_htVecLepton_ptR", 0.)
            self.out.fillBranch(self.outputName+"_htSLepton_ptR", 0.)

            self.out.fillBranch(self.outputName+"_minPhiStar", 0.)
            self.out.fillBranch(self.outputName+"_ht", 0.)
            self.out.fillBranch(self.outputName+"_hmass", 0.)
            self.out.fillBranch(self.outputName+"_mht", 0.)
            self.out.fillBranch(self.outputName+"_mht_phi", 0.)
            self.out.fillBranch(self.outputName+"_mht_met_dphi", 0.)

        
        if len(looseLeptons)>0 and len(jets)>0:
            looseLepton = looseLeptons[0]
            self.out.fillBranch(self.outputName+"_minJetSubLepton_dphi", min(map(lambda j: math.fabs(deltaPhi(looseLepton.phi,j.phi)), jets)))
            self.out.fillBranch(self.outputName+"_minJetSubLepton_dR", min(map(lambda j: deltaR(looseLepton,j), jets)))
        else:
            self.out.fillBranch(self.outputName+"_minJetSubLepton_dphi", 10.)
            self.out.fillBranch(self.outputName+"_minJetSubLepton_dR", 10.)
            


        return True
