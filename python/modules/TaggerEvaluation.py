import os
import sys
import math
import json
import ROOT
import random
import numpy
import time
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import getCtauLabel

class TaggerEvaluation(Module):

    def __init__(
        self,
        modelPath,
        inputCollections = [lambda event: Collection(event, "Jet")],
        taggerName = "llpdnnx",
        predictionLabels = ["B","C","UDS","G","PU","LLP_QMU","LLP_Q"], #this is how the output array from TF is interpreted
        logctauValues = range(-2, 4),
        globalOptions = {"isData":False},
    ):
        self.globalOptions = globalOptions
        self.inputCollections = inputCollections
        self.predictionLabels = predictionLabels
        self.logctauValues = logctauValues
        self.logctau = numpy.array(logctauValues,dtype=numpy.float32)
        
        self.modelPath = modelPath
        self.taggerName = taggerName
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.setup(inputTree)
                

    def setupTFEval(self,tree,modelFile,featureDict):
        print "Building TFEval object"
        tfEval = ROOT.TFEval()
        print "Succesfully built TFEval object"

        if (not tfEval.loadGraph(modelFile)):
            sys.exit(1)
        tfEval.addOutputNodeName("prediction")
        print "--- Model: ",modelFile," ---"
        for groupName,featureCfg in featureDict.iteritems():
            if featureCfg.has_key("max"):
                print "building group ... %s, shape=[%i,%i], length=%s"%(groupName,featureCfg["max"],len(featureCfg["branches"]),featureCfg["length"])
                lengthBranch = ROOT.TFEval.createAccessor(tree.arrayReader(featureCfg["length"]))
                featureGroup = ROOT.TFEval.ArrayFeatureGroup(
                    groupName,
                    len(featureCfg["branches"]),
                    featureCfg["max"],
                    lengthBranch
                )
                for branchName in featureCfg["branches"]:
                    print branchName
                    featureGroup.addFeature(ROOT.TFEval.createAccessor(tree.arrayReader(branchName)))
                tfEval.addFeatureGroup(featureGroup)
            else:
                print "building group ... %s, shape=[%i]"%(groupName,len(featureCfg["branches"]))
                featureGroup = ROOT.TFEval.ValueFeatureGroup(
                    groupName,
                    len(featureCfg["branches"])
                )
                for branchName in featureCfg["branches"]:
                    print branchName
                    featureGroup.addFeature(ROOT.TFEval.createAccessor(tree.arrayReader(branchName)))
                tfEval.addFeatureGroup(featureGroup)
                
        return tfEval
        
    def setup(self,tree):
        #load dynamically from file
        print "loading feature dict..."
        featureDict = import_module('feature_dict').featureDict
        print "imported featurel dict..."
        self.tfEvalParametric = self.setupTFEval(tree,self.modelPath,featureDict)
        print "setup model successfully..."
        
        genFeatureGroup = ROOT.TFEval.ValueFeatureGroup("gen",1)
        self.nJets = 0
        genFeatureGroup.addFeature(ROOT.TFEval.PyAccessor(lambda: self.nJets, lambda jetIndex,batchIndex: self.logctau[batchIndex%len(self.logctau)]))
        self.tfEvalParametric.addFeatureGroup(genFeatureGroup)
        
        self._ttreereaderversion = tree._ttreereaderversion
        
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):

        jetglobal = Collection(event, "global")
        
        jetOriginIndices = set() #superset of all indices to evaluate

        for jetCollection in self.inputCollections:
            jets = jetCollection(event)
            for ijet,jet in enumerate(jets):
                if jet._index>=len(jetglobal):
                    print "Jet not filled"
                    continue
                jetOriginIndices.add(jet._index)
                
        jetOriginIndices = list(jetOriginIndices)
        
        evaluationIndices = []
        for index in jetOriginIndices:
            evaluationIndices.extend([index]*len(self.logctauValues))
                
        if event._tree._ttreereaderversion > self._ttreereaderversion:
            self.setup(event._tree)
            
        self.nJets = len(jetglobal)
        if len(jetOriginIndices)==0:
            return True
            
        evaluationIndices = numpy.array(evaluationIndices,numpy.int64)
        result = self.tfEvalParametric.evaluate(
            evaluationIndices.shape[0],
            evaluationIndices
        )

        predictionsPerIndexAndCtau = {}

        for ijet,jetIndex in enumerate(jetOriginIndices):
            predictionsPerIndexAndCtau[jetIndex] = {}

            for ictau,ctau in enumerate(self.logctauValues):
                predictionIndex = ijet*len(self.logctauValues)+ictau
                predictionsPerIndexAndCtau[jetIndex][ctau] = result.get("prediction",predictionIndex)
                
                
        for jetCollection in self.inputCollections:
            jets = jetCollection(event)

            for ijet,jet in enumerate(jets):
                taggerOutput = {}

                for ictau,ctau in enumerate(self.logctauValues):
                    taggerOutput[self.logctauValues[ictau]] = {}

                    for iclass, classLabel in enumerate(self.predictionLabels):  
                        if jet._index<len(jetglobal):
                            taggerOutput[self.logctauValues[ictau]][classLabel] = \
                                    predictionsPerIndexAndCtau[jet._index][ctau][iclass]
                        else:
                            taggerOutput[self.logctauValues[ictau]][classLabel] = -1

                setattr(jet, self.taggerName, taggerOutput)
        return True
            
