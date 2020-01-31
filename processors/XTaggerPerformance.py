import os
import sys
import math
import json
import ROOT
import random
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from PhysicsTools.NanoAODTools.modules import *

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--isData', dest='isData', action='store_true',default=False)
parser.add_argument('--year', dest='year', action='store',default=2016, type=int)
parser.add_argument('--input', dest='inputFiles', action='append',default=[])
parser.add_argument('output', nargs=1)

args = parser.parse_args()

print "isData:",args.isData
print "year:",args.year
print "inputs:",len(args.inputFiles)
for inputFile in args.inputFiles:
    rootFile = ROOT.TFile.Open(inputFile)
    if not rootFile:
        print "CRITICAL - file '"+inputFile+"' not found!"
        sys.exit(1)
    tree = rootFile.Get("Events")
    if not tree:
        print "CRITICAL - 'Events' tree not found in file '"+inputFile+"'!"
        sys.exit(1)
    print " - ",inputFile,", events=",tree.GetEntries()
    
print "output directory:",args.output[0]

globalOptions = {
    "isData":False,
    "year":args.year
}
analyzerChain = []

analyzerChain.append(
    JetSelection(
        outputName="selectedJets",
        storeKinematics=['pt','eta'],
        genJetMinPt = 5.,
        minRatio=0.5
    )
)

analyzerChain.append(
    JetTruthFlags(
        inputCollection = lambda event: event.selectedJets,
        outputName = "selectedJets",
        latentVariables = ["llp_mass", "llp_pt", "displacement", "displacement_xy", "displacement_z", "decay_angle", "betagamma"]
    )   
)
    
analyzerChain.append(
    EventSkim(selection=lambda event: 
        len(event.selectedJets)>0
    )
)

analyzerChain.append(
    TaggerEvaluation(
        modelPath="PhysicsTools/NanoAODTools/data/nn/weight2016_75.pb",
        logctauValues = [1.74],
        #modelPath="PhysicsTools/NanoAODTools/data/nn/da.pb",
        inputCollections=[
            lambda event: event.selectedJets,
        ],
        taggerName="llpdnnx",
    )
)


analyzerChain.append(
    JetTaggerResult(
        inputCollection = lambda event: event.selectedJets,
        taggerName = "llpdnnx",
        logctauValues = [1.74],
        predictionLabels = ["LLP_Q", "LLP_QMU"],
    )
)

 

p=PostProcessor(
    args.output[0],
    [args.inputFiles],
    modules=analyzerChain,
    maxEvents=-1,
    friend=True
)
p.run()
