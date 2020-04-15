import os
import sys
import math
import json
import ROOT
import random
import numpy
import time
import imp

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import getCtauLabel

class TaggerEvaluation(Module):

    def __init__(
        self,
        modelPath,
        featureDictFile,
        inputCollections = [lambda event: Collection(event, "Jet")],
        taggerName = "llpdnnx",
        predictionLabels = ["B","C","UDS","G","PU","isLLP_QMU_QQMU","isLLP_Q_QQ"], #this is how the output array from TF is interpreted
        evalValues = range(-1, 4),
        globalOptions = {"isData":False},
    ):
        self.globalOptions = globalOptions
        self.inputCollections = inputCollections
        self.predictionLabels = predictionLabels
        self.evalValues = evalValues
        self.logctau = numpy.array(evalValues,dtype=numpy.float32)
        
        self.modelPath = os.path.expandvars(modelPath)
        print featureDictFile
        self.featureDict = imp.load_source(
            'feature_dict', 
            os.path.expandvars(featureDictFile)
        ).featureDict
        self.taggerName = taggerName
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.setup(inputTree)
                

    def setupTFEval(self,tree,modelFile):
        print "Building TFEval object"
        tfEval = ROOT.TFEval()
        print "Succesfully built TFEval object"

        if (not tfEval.loadGraph(modelFile)):
            sys.exit(1)
        tfEval.addOutputNodeName("prediction")
        print "--- Model: ",modelFile," ---"
        for groupName,featureCfg in self.featureDict.iteritems():
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
        self.tfEvalParametric = self.setupTFEval(tree,self.modelPath)
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
        jetglobal_indices = [global_jet.jetIdx for global_jet in jetglobal]

        jetOriginIndices = set() # superset of all indices to evaluate

        for jetCollection in self.inputCollections:
            jets = jetCollection(event)
            for ijet, jet in enumerate(jets):
                global_jet_index = jetglobal_indices.index(jet._index)
                global_jet = jetglobal[global_jet_index]

                if abs(jet.eta - global_jet.eta) < 0.01 or 
                   abs(jet.phi - global_jet_phi) < 0.01:
                       print "Warning ->> jet might be mismatched!"
                jetOriginIndices.add(global_jet_index)
                setattr(jet, "globalIdx", global_jet_index)
                
        jetOriginIndices = list(jetOriginIndices)

        
        evaluationIndices = []
        for index in jetOriginIndices:
            evaluationIndices.extend([index]*len(self.evalValues))
                
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

            for ictau,ctau in enumerate(self.evalValues):
                predictionIndex = ijet*len(self.evalValues)+ictau
                predictionsPerIndexAndCtau[jetIndex][ctau] = result.get("prediction",predictionIndex)
                
                
        for jetCollection in self.inputCollections:
            jets = jetCollection(event)

            for ijet, jet in enumerate(jets):
                taggerOutput = {}

                for ictau, ctau in enumerate(self.evalValues):
                    taggerOutput[self.evalValues[ictau]] = {}

                    for iclass, classLabel in enumerate(self.predictionLabels):  
                            taggerOutput[self.evalValues[ictau]][classLabel] = \
                                    predictionsPerIndexAndCtau[jet.globalIdx][ctau][iclass]

                setattr(jet, self.taggerName, taggerOutput)
        return True

