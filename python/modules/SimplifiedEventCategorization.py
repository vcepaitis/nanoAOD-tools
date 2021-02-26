import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from utils import deltaR

class SimplifiedEventCategorization(Module):
    def __init__(
        self,
        globalOptions={"isData":False, "isSignal":False},
        outputName="category",
        looseLeptons=None,
        jetsCollection=None,
        maxDeltaR=1.3,
        flags={
            'isLLP_Q': ['isLLP_RAD','isLLP_Q','isLLP_QQ', 'isLLP_TAU', 'isLLP_QTAU', 'isLLP_QQTAU'],
            'isLLP_QMU': ['isLLP_QMU', 'isLLP_QQMU', 'isLLP_MU'],
            'isLLP_QE': ['isLLP_QE','isLLP_QQE', 'isLLP_E'],
        },

    ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.looseLeptons = looseLeptons
        self.jetsCollection = jetsCollection
        self.maxDeltaR = maxDeltaR
        self.flags = flags

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.outputName+"_isResolved", "I")
        self.out.branch(self.outputName+"_isMerged", "I")

        if self.globalOptions["isSignal"]:
            self.out.branch(self.outputName+"_isTrueResolved", "I")
            self.out.branch(self.outputName+"_isTrueMerged", "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        jets = self.jetsCollection(event)
        looseLeptons = self.looseLeptons(event)

        isMerged = 0
        isResolved = 0

        if len(jets) > 0 and len(looseLeptons) > 0:
            minDeltaR = min([deltaR(jet, looseLeptons[0]) for jet in jets])
            if minDeltaR < 0.4:
                isMerged = 1
            elif (minDeltaR < self.maxDeltaR):
                isResolved = 1

        self.out.fillBranch(self.outputName+"_isMerged", isMerged)
        self.out.fillBranch(self.outputName+"_isResolved", isResolved)

        if self.globalOptions["isSignal"]:
            nTrueResJets = 0
            nTrueLepJets = 0
            isTrueMerged = 0
            isTrueResolved = 0

            jetOrigin = Collection(event, "jetorigin")

            for jet in jets:
                for truth_class, subclass_dict in self.flags.iteritems():
                    setattr(jet, truth_class, 0)
                    for subclass in subclass_dict:
                        if (getattr(jetOrigin[jet._index], subclass) > 0.5):
                            setattr(jet, truth_class, 1)

                if jet.isLLP_Q > 0:
                    nTrueResJets += 1
                if jet.isLLP_QMU > 0 or jet.isLLP_QE > 0:
                    nTrueLepJets += 1

            if nTrueLepJets > 0:
                isTrueMerged = 1
            elif nTrueResJets > 0:
                isTrueResolved = 1


            self.out.fillBranch(self.outputName+"_isTrueMerged", isMerged)
            self.out.fillBranch(self.outputName+"_isTrueResolved", isResolved)

        return True
