import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class ElectronSelection(Module):
    VETO = 0
    LOOSE = 1
    MEDIUM = 2
    TIGHT = 3

    def __init__(
        self,
        inputCollection = lambda event: Collection(event, "Electron"),
        outputName = "tightElectrons",
        triggerMatch = False,
        electronID = TIGHT,
        electronMinPt = 5.,
        electronMaxEta = 2.4,
        storeKinematics=['pt','eta'],
        storeWeights=False,
        selectLeadingOnly=False,
        globalOptions={"isData":False, "year":2016}
    ):
        
        self.globalOptions = globalOptions
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.electronMinPt = electronMinPt
        self.electronMaxEta = electronMaxEta
        self.storeKinematics = storeKinematics
        self.storeWeights = storeWeights
        self.selectLeadingOnly = selectLeadingOnly
        self.triggerMatch = triggerMatch
        self.electronID = electronID

        if triggerMatch:
            self.trigger_object = lambda event: Collection(event, "TrigObj")
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName, "I")

        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
            
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = self.inputCollection(event)
        
        selectedElectrons = []
        unselectedElectrons = []
        
        for electron in electrons:
            if electron.pt>self.electronMinPt and math.fabs(electron.eta)<self.electronMaxEta and (electron.cutBased>self.electronID):
                selectedElectrons.append(electron)
            else:
                unselectedElectrons.append(electron)

        if len(selectedElectrons) > 0:
            if self.selectLeadingOnly:
                unselectedElectrons.extend(selectedElectrons[1:])
                selectedElectrons = [selectedElectrons[0]]


  
        self.out.fillBranch("n"+self.outputName,len(selectedElectrons))
        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName+"_"+variable,map(lambda electron: getattr(electron,variable),selectedElectrons))
 

        setattr(event,self.outputName,selectedElectrons)
        setattr(event,self.outputName+"_unselected",unselectedElectrons)

        return True
        
