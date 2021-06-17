import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class PDFUncertainty(Module):
    def __init__(
        self,
        storeVariables = []
    ):
        self.storeVariables = storeVariables
        
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch("pdf_nominal","F")
        self.out.branch("pdf_up","F")
        self.out.branch("pdf_down","F")
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        pdfidx = range(1,102)
        values = np.zeros(len(pdfidx))
        for i,idx in enumerate(pdfidx):
            values[i] = getattr(event,"LHEWeights_nnpdfreplica_%i"%idx)
        down,up = np.percentile(values[1:],[15.8,84.1]) #exclude nominal weight at [0]
        
        self.out.fillBranch("pdf_nominal",1.)
        self.out.fillBranch("pdf_up",up)
        self.out.fillBranch("pdf_down",down)
        
        return True
        
