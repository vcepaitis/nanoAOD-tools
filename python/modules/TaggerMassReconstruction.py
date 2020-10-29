import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class TaggerMassReconstruction(Module):
    def __init__(
        self,
        globalOptions={"isData":False, "isSignal":False},
        outputName="category_simplified",
        tightLeptons=None,
        looseLeptons=None,
        lepJets=None,
        resJets=None,
        taggerName="llpdnnx",
        jetLabels=['LLP_Q','LLP_QE','LLP_QMU']

    ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.tightLeptons = tightLeptons
        self.looseLeptons = looseLeptons
        self.lepJets = lepJets
        self.resJets = resJets
        self.taggerName = taggerName
        self.jetLabels = jetLabels

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self.outputName+"_"+self.taggerName+"_m_llj", "F")
        self.out.branch(self.outputName+"_"+self.taggerName+"_m_lljj", "F")
        #self.out.branch(self.outputName+"_m_HNL", "F")
        self.out.branch(self.outputName+"_"+self.taggerName+"_max", "F")
        self.out.branch(self.outputName+"_"+self.taggerName+"_max2nd", "F")

        '''
        for label in self.jetLabels:
            self.out.branch("{}_lepJet_{}_{}".format(self.outputName, self.taggerName, label), "F", lenVar="n"+self.outputName+"_lepJets")
            self.out.branch("{}_resJet_{}_{}".format(self.outputName, self.taggerName, label), "F", lenVar="n"+self.outputName+"_resJets")
        '''

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        tightLeptons = self.tightLeptons(event)
        looseLeptons = self.looseLeptons(event)

        lepJets = getattr(event, self.lepJets)
        resJets = getattr(event, self.resJets)

        jets = lepJets + resJets

        for jet in jets:
            taggerScore = getattr(jet, self.taggerName)
            for label in self.jetLabels:
                setattr(jet, "{}_{}".format(self.taggerName, label), taggerScore[label])

        if len(lepJets) > 0:
            lepJets = sorted(lepJets, key=lambda lepJet: \
                max(getattr(lepJet, self.taggerName+"_LLP_QMU"), getattr(lepJet, self.taggerName+"_LLP_QE")),
                reverse=True)
        if len(resJets) > 0:
            resJets = sorted(resJets, key=lambda resJet: getattr(resJet, self.taggerName+"_LLP_Q"), reverse=True)

        nlepJets = len(lepJets)
        nresJets = len(resJets)
        nJets = nlepJets+nresJets

        '''
        for label in self.jetLabels:
            self.out.fillBranch("{}_lepJet_{}_{}".format(self.outputName, self.taggerName, label), [getattr(jet, "{}_{}".format(self.taggerName, label)) for jet in lepJets])
            self.out.fillBranch("{}_resJet_{}_{}".format(self.outputName, self.taggerName, label), [getattr(jet, "{}_{}".format(self.taggerName, label)) for jet in resJets])
        '''

        WCandidateLorentzVector = ROOT.TLorentzVector(0,0,0,0)

        maxScore = -1
        max2ndScore = -1.

        for lepton in looseLeptons+tightLeptons:
            WCandidateLorentzVector += lepton.p4()

        if len(lepJets) > 0:
            WCandidateLorentzVector += lepJets[0].p4()
            maxScore = max(lepJets[0].llpdnnx_LLP_QMU, lepJets[0].llpdnnx_LLP_QE)

        elif len(resJets) > 0:
            WCandidateLorentzVector += resJets[0].p4()
            maxScore = resJets[0].llpdnnx_LLP_Q

        self.out.fillBranch(self.outputName+"_"+self.taggerName+"_m_llj", WCandidateLorentzVector.M())
        self.out.fillBranch(self.outputName+"_"+self.taggerName+"_max", maxScore)

        if len(lepJets) == 0 and len(resJets) > 1:
            WCandidateLorentzVector += resJets[1].p4()
            max2ndScore = resJets[1].llpdnnx_LLP_Q

        self.out.fillBranch(self.outputName+"_"+self.taggerName+"_m_lljj", WCandidateLorentzVector.M())
        self.out.fillBranch(self.outputName+"_"+self.taggerName+"_max2nd", max2ndScore)

        return True
