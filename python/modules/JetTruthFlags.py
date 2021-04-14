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
        originVariables = [],
        globalVariables = [],
        flags={
            'isPrompt_MU': ['isPrompt_MU'],
            'isPrompt_E': ['isPrompt_E'],
            'isPrompt_TAU': ['isPrompt_TAU'],
            'isPrompt_PHOTON': ['isPrompt_PHOTON'],
            'isB': ['isB', 'isBB', 'isLeptonic_B', 'isLeptonic_C'],
            'isC': ['isC', 'isCC'],
            'isUDS': ['isS', 'isUD'],
            'isG': ['isG'],
            'isPU': ['isPU'],
            'isLLP_Q': ['isLLP_RAD','isLLP_Q','isLLP_QQ'],
            'isLLP_MU': ['isLLP_MU'],
            'isLLP_QMU': ['isLLP_QMU', 'isLLP_QQMU'],
            'isLLP_E': ['isLLP_E'],
            'isLLP_QE': ['isLLP_QE','isLLP_QQE'],
            'isLLP_TAU': ['isLLP_TAU'],
            'isLLP_QTAU': ['isLLP_QTAU','isLLP_QQTAU'],  
            'isUndefined': ['isUndefined'] , 
        },
        genVariables = ['displacement','displacement_xy','displacement_z'],

        globalOptions={"isData": False, "isSignal": False}
    ):
        self.globalOptions = globalOptions
        self.flags = flags
        self.inputCollection = inputCollection
        self.originVariables = originVariables
        self.globalVariables = globalVariables
        self.genVariables = genVariables
        self.outputName = outputName

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if not self.globalOptions['isData']:
            for originVariable in self.originVariables+self.globalVariables:
                self.out.branch(self.outputName+"_"+originVariable, "F",
                                lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        if self.globalOptions['isData']:
            return True

        jets = self.inputCollection(event)
        jetOrigin = Collection(event, "jetorigin")
        jetGlobal = Collection(event, "global")

        extraVariableDict = {}

        for extraVariable in self.originVariables+self.globalVariables:
            extraVariableDict[extraVariable] = [-1.]*len(jets)

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
                
            if jet._index<len(jetOrigin):
                for originVariable in self.originVariables:
                    extraVariableDict[originVariable][ijet] = getattr(
                        jetOrigin[jet._index], originVariable)
                for genVariable in self.genVariables:
                    setattr(
                        jet,
                        genVariable,
                        getattr(jetOrigin[jet._index], genVariable)
                    )
                    
            if jet._index<len(jetGlobal):
                for globalVariable in self.globalVariables:
                    extraVariableDict[globalVariable][ijet] = getattr(
                        jetGlobal[jet._index], globalVariable)

        for k in sorted(self.flags.keys()):
            self.out.fillBranch(self.outputName+"_"+k, flavors[k])

        for extraVariable in self.originVariables+self.globalVariables:
            self.out.fillBranch(self.outputName+"_"+extraVariable,
                                extraVariableDict[extraVariable])

        return True
