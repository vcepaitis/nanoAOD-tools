import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from utils import deltaR, getCtauLabel, deltaPhi

class SimplifiedEventCategorization(Module):
    def __init__(
        self,
        globalOptions={"isData":False, "isSignal":False},
        outputName="category_simplified",
        looseLeptons=None,
        jetsCollection=None,
        maxDeltaR=1.3,
        jetLabels=['LLP_Q','LLP_QE','LLP_QMU'],
        flags={
            'isPrompt_MU': ['isPrompt_MU'],
            'isPrompt_E': ['isPrompt_E'],
            'isPrompt_TAU': ['isPrompt_TAU'],
            'isB': ['isB', 'isBB',  'isLeptonic_B', 'isLeptonic_C'],
            'isC': ['isC', 'isCC'],
            'isUDS': ['isS', 'isUD'],
            'isG': ['isG'],
            'isPU': ['isPU'],
            'isLLP_Q': ['isLLP_RAD','isLLP_Q','isLLP_QQ', 'isLLP_TAU', 'isLLP_QTAU', 'isLLP_QQTAU'],
            'isLLP_QMU': ['isLLP_QMU', 'isLLP_QQMU', 'isLLP_MU'],
            'isLLP_QE': ['isLLP_QE','isLLP_QQE', 'isLLP_E'],
            'isUndefined': ['isUndefined'],
        },

    ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.looseLeptons = looseLeptons
        self.jetsCollection = jetsCollection
        self.maxDeltaR = maxDeltaR
        self.jetLabels = jetLabels
        self.flags = flags

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.outputName+"_index", "I")
        self.out.branch("n"+self.outputName+"_Jets", "I")
        self.out.branch("n"+self.outputName+"_lepJets", "I")
        self.out.branch("n"+self.outputName+"_resJets", "I")

        if self.globalOptions["isSignal"]:
            self.out.branch(self.outputName+"_index_truth", "I")
            self.out.branch("n"+self.outputName+"_TrueJets", "I")
            self.out.branch("n"+self.outputName+"_TrueLepJets", "I")
            self.out.branch("n"+self.outputName+"_TrueResJets", "I")

            for truth_class in self.flags.keys():
                self.out.branch(self.outputName+"_Jet_"+truth_class, "I", lenVar="n"+self.outputName+"_Jets")
                self.out.branch(self.outputName+"_lepJet_"+truth_class, "I", lenVar="n"+self.outputName+"_lepJets")
                self.out.branch(self.outputName+"_resJet_"+truth_class, "I", lenVar="n"+self.outputName+"_resJets")
        self.out.branch(self.outputName+"_Jet_category", "I", lenVar="n"+self.outputName+"_Jets")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        jets = self.jetsCollection(event)
        looseLeptons = self.looseLeptons(event)
        if self.globalOptions["isSignal"]:
            jetOrigin = Collection(event, "jetorigin")

        lepJets = []
        resJets = []

        for jet in jets:
            if self.globalOptions["isSignal"]:
                for truth_class, subclass_dict in self.flags.iteritems():
                    setattr(jet, truth_class, 0)
                    for subclass in subclass_dict:
                        if (getattr(jetOrigin[jet._index], subclass) > 0.5):
                            setattr(jet, truth_class, 1)

        if len(looseLeptons) == 0:
            resJets = jets

        elif len(jets) > 0:
            for jet in jets:
                minDeltaR = min([deltaR(jet, lepton) for lepton in looseLeptons])
                if (minDeltaR) < 0.4:
                    lepJets.append(jet)
                elif (minDeltaR < self.maxDeltaR):
                    resJets.append(jet)

        nlepJets = len(lepJets)
        nresJets = len(resJets)
        Jets = lepJets+resJets
        nJets = nlepJets+nresJets

        self.out.fillBranch("n"+self.outputName+"_Jets", nJets)
        self.out.fillBranch("n"+self.outputName+"_lepJets", nlepJets)
        self.out.fillBranch("n"+self.outputName+"_resJets", nresJets)

        setattr(event, self.outputName+"_Jets", Jets)
        setattr(event, self.outputName+"_lepJets", lepJets)
        setattr(event, self.outputName+"_resJets", resJets)

        if len(looseLeptons) > 0:
            if nJets > 0:
                if nlepJets > 0:
                    category_index = 2
                else:
                    category_index = 1
            else:
                category_index = -1
        else:
            if nJets == 1:
                category_index = 3
            elif nJets > 1:
                category_index = 4
            else:
                category_index = -1

        self.out.fillBranch(self.outputName+"_index", category_index)

        if self.globalOptions["isSignal"]:
            for truth_class in self.flags.keys():
                self.out.fillBranch(self.outputName+"_Jet_"+truth_class, [getattr(jet, truth_class) for jet in Jets])
                self.out.fillBranch(self.outputName+"_lepJet_"+truth_class, [getattr(jet, truth_class) for jet in lepJets])
                self.out.fillBranch(self.outputName+"_resJet_"+truth_class, [getattr(jet, truth_class) for jet in resJets])

            nTrueResJets = 0
            nTrueLepJets = 0
            category_index_truth = -1

            for jet in jets:
                if jet.isLLP_Q > 0:
                    nTrueResJets += 1
                if jet.isLLP_QMU > 0 or jet.isLLP_QE > 0:
                    nTrueLepJets += 1

            nTrueJets = nTrueLepJets+nTrueResJets

            if nTrueJets == 1:
                if nTrueResJets == 1:
                    category_index_truth = 1
                elif nTrueLepJets == 1:
                    category_index_truth = 2
            elif nTrueJets > 1:
                if nTrueLepJets > 0:
                    category_index_truth = 4
                else:
                    category_index_truth = 3
            self.out.fillBranch("n"+self.outputName+"_TrueJets", nTrueJets)
            self.out.fillBranch("n"+self.outputName+"_TrueLepJets", nTrueLepJets)
            self.out.fillBranch("n"+self.outputName+"_TrueResJets", nTrueResJets)
            self.out.fillBranch(self.outputName+"_index_truth", category_index_truth)

        self.out.fillBranch(self.outputName+"_Jet_category", nJets*[category_index])
        return True
