import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from utils import deltaR, getCtauLabel

class SimplifiedEventCategorization(Module):
    def __init__(
        self,
        globalOptions={"isData":False, "isSignal":False},
        outputName="category_simplified",
        looseLeptons=None,
        jetsCollection=None,
        taggerName="llpdnnx",
        maxDeltaR=1.3,
        jetLabels=['LLP_Q','LLP_QE','LLP_QMU'],
        flags={
            'isPrompt_MU': ['isPrompt_MU'],
            'isPrompt_E': ['isPrompt_E'],
            'isPrompt_TAU': ['isPrompt_TAU'],
            'isB': ['isB', 'isBB', 'isGBB', 'isLeptonic_B', 'isLeptonic_C'],
            'isC': ['isC', 'isCC', 'isGCC'],
            'isUDS': ['isS', 'isUD'],
            'isG': ['isG'],
            'isPU': ['isPU'],
            'isLLP_Q': ['isLLP_RAD','isLLP_Q','isLLP_QQ'],
            'isLLP_QMU': ['isLLP_QMU', 'isLLP_QQMU', 'isLLP_MU'],
            'isLLP_QE': ['isLLP_QE','isLLP_QQE', 'isLLP_E'],
            'isLLP_TAU': ['isLLP_TAU'],
            'isLLP_QTAU': ['isLLP_QTAU','isLLP_QQTAU'],
            'isUndefined': ['isUndefined'],
        },
    ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.looseLeptons = looseLeptons
        self.jetsCollection = jetsCollection
        self.taggerName = taggerName
        self.jetLabels = jetLabels
        self.maxDeltaR = maxDeltaR
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
        for label in self.jetLabels:
            self.out.branch("{}_lepJet_{}_{}".format(self.outputName, self.taggerName, label), "F", lenVar="n"+self.outputName+"_lepJets")
            self.out.branch("{}_resJet_{}_{}".format(self.outputName, self.taggerName, label), "F", lenVar="n"+self.outputName+"_resJets")

        if self.globalOptions["isSignal"]:
            for truth_class in self.flags.keys():
                self.out.branch(self.outputName+"_lepJet_"+truth_class, "I", lenVar="n"+self.outputName+"_lepJets")
                self.out.branch(self.outputName+"_resJet_"+truth_class, "I", lenVar="n"+self.outputName+"_resJets")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        jets = self.jetsCollection(event)
        looseLeptons = self.looseLeptons(event)

        # only for MC: truth labels
        if self.globalOptions["isSignal"]:
            jetOrigin = Collection(event, "jetorigin")
            indexFlag = []
            indexFlag.append(-10.)

        lepJets = []
        resJets = []

        for jet in jets:
            taggerScore = getattr(jet, self.taggerName)
            for label in self.jetLabels:
                setattr(jet, "{}_{}".format(self.taggerName, label), taggerScore[label])

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
                if(minDeltaR) < 0.4:
                    lepJets.append(jet)
                elif (minDeltaR < self.maxDeltaR):
                    resJets.append(jet)

        if len(lepJets) > 0:
            lepJets = sorted(lepJets, key=lambda lepJet: \
                max(getattr(lepJet, self.taggerName+"_LLP_QMU"), getattr(lepJet, self.taggerName+"_LLP_QE")),
                reverse=True)
        if len(resJets) > 0:
            resJets = sorted(resJets, key=lambda resJet: getattr(resJet, self.taggerName+"_LLP_Q"), reverse=True)

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
            if nJets == 1:
                if nresJets == 1:
                    category_index = 1
                elif nlepJets == 1:
                    category_index = 2
            elif nJets > 1:
                if nlepJets > 0:
                    category_index = 4
                else:
                    category_index = 3
            else:
                category_index = -1
        else:
            if nresJets == 1:
                category_index = 5
            elif nresJets > 1:
                category_index = 6
            else:
                category_index = -1

        self.out.fillBranch(self.outputName+"_index", category_index)

        for label in self.jetLabels:
            self.out.fillBranch("{}_lepJet_{}_{}".format(self.outputName, self.taggerName, label), [getattr(jet, "{}_{}".format(self.taggerName, label)) for jet in lepJets])
            self.out.fillBranch("{}_resJet_{}_{}".format(self.outputName, self.taggerName, label), [getattr(jet, "{}_{}".format(self.taggerName, label)) for jet in resJets])

        if self.globalOptions["isSignal"]:
            for truth_class in self.flags.keys():
                self.out.fillBranch(self.outputName+"_lepJet_"+truth_class, [getattr(jet, truth_class) for jet in lepJets])
                self.out.fillBranch(self.outputName+"_resJet_"+truth_class, [getattr(jet, truth_class) for jet in resJets])

        return True
