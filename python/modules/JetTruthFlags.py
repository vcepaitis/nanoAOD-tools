import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel \
    import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR


class JetTruthFlags(Module):
    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "Jet"),
        outputName="selectedJets_nominal",
        latentVariables=[],
        flags={
            'isB': ['isB', 'isBB', 'isGBB', 'isLeptonic_B', 'isLeptonic_C'],
            'isC': ['isC', 'isCC', 'isGCC'],
            'isUDS': ['isS', 'isUD'],
            'isG': ['isG'],
            'isPU': ['isPU'],
            'isLLP_QMU': ['isLLP_QMU'],
            'isLLP_QQMU': ['isLLP_QQMU'],
            'isLLP_Q': ['isLLP_Q'],
            'isLLP_QQ': ['isLLP_QQ'],
            'isLLP_MU': ['isLLP_MU'],
            'isLLP_Merged': ['isLLP_QMU', 'isLLP_QQMU'],
            'isLLP_Resolved': ['isLLP_Q', 'isLLP_QQ'],
            'isLLP_Q': ['isLLP_Q'],
            'isLLP_QQ': ['isLLP_QQ'],
            'isLLP_MU': ['isLLP_MU'],
            'isLLP_RAD': ['isLLP_RAD']

        },
        globalOptions={"isData": False}
    ):
        self.globalOptions = globalOptions
        self.flags = flags
        self.inputCollection = inputCollection
        self.latentVariables = latentVariables
        self.outputName = outputName

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if not self.globalOptions['isData']:
            for latentVariable in self.latentVariables:
                self.out.branch(self.outputName+"_"+latentVariable, "F",
                                lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        if self.globalOptions['isData']:
            return True

        jets = self.inputCollection(event)
        jetOrigin = Collection(event, "jetorigin")

        latentVariableDict = {}

        for latentVariable in self.latentVariables:
            latentVariableDict[latentVariable] = [-1.]*len(jets)

        flavors = {}

        for k in sorted(self.flags.keys()):
            flavors[k] = [-1.]*len(jets)
            self.out.branch(self.outputName+"_"+k, "F",
                            lenVar="n"+self.outputName)

        for ijet, jet in enumerate(jets):
            for k in sorted(self.flags.keys()):
                flavorFlag = 0.
                for originFlag in self.flags[k]:
                    flagValue = getattr(jetOrigin[jet._index], originFlag)
                    if (flagValue > 0.5):
                        flavorFlag = 1.
                        break
                flavors[k][ijet] = flavorFlag
                setattr(jet, k, flavorFlag > 0.5)
            for latentVariable in self.latentVariables:
                latentVariableDict[latentVariable][ijet] = getattr(
                    jetOrigin[jet._index], latentVariable)

        for k in sorted(self.flags.keys()):
            self.out.fillBranch(self.outputName+"_"+k, flavors[k])

        if not self.globalOptions['isData']:
            for latentVariable in self.latentVariables:
                self.out.fillBranch(self.outputName+"_"+latentVariable,
                                    latentVariableDict[latentVariable])

        return True
