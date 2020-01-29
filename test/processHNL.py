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

from modules import *


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--isData', dest='isData', action='store_true',default=False)
parser.add_argument('--year', dest='year', action='store',type=int, default=2016)
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
    "isData":args.isData,
    "year":args.year
}


muonSelection = [
    MuonSelection(
        outputName="tightMuons",
        storeKinematics=['pt','eta', 'dxy', 'dxyErr', 'dz', 'dzErr', 'phi'],
        storeWeights=True,
        muonMinPt = 25.,
        muonMaxDxy = 0.002,
        muonMaxDz = 0.01,
        muonID = MuonSelection.TIGHT,
        muonIso = MuonSelection.TIGHT,
        globalOptions=globalOptions
    ),
    MuonSelection(
        inputCollection = lambda event: event.tightMuons_unselected,
        outputName="looseMuons",
        storeKinematics=['pt','eta', 'dxy', 'dxyErr', 'dz', 'dzErr', 'phi'],
        storeWeights=False,
        muonMinPt = 5.,
        muonID = MuonSelection.LOOSE,
        muonIso = MuonSelection.NONE,
        globalOptions=globalOptions
    ),
    
    SingleMuonTriggerSelection(
        inputCollection=lambda event: event["tightMuons"],
        outputName="IsoMuTrigger",
        storeWeights=False,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.IsoMuTrigger_flag==1),
    EventSkim(selection=lambda event: event.ntightMuons>0),
    EventSkim(selection=lambda event: event.nlooseMuons>0)
]

analyzerChain = []

analyzerChain.extend(muonSelection)

analyzerChain.append(
    JetSelection(
        outputName="selectedJets",
        jetMinPt = 20.,
        jetMaxEta = 2.4,
        storeKinematics=['pt','eta', 'phi', 'mass', 'nMuons', 'nElectrons', 'muonSubtrFactor'],
    )
)


storeVariables = [
    [lambda tree: tree.branch("genweight","F"),lambda tree,event: tree.fillBranch("genweight",event.Generator_weight)],
]


analyzerChain.append(EventInfo(storeVariables=storeVariables))

analyzerChain.append(EventSkim(selection=lambda event: len(event.selectedJets)>0))


analyzerChain.append(
    InvariantSystem(
        inputCollection = lambda event: [event.looseMuons[0], event.tightMuons[0]],
        outputName = "dimuon"
    )
)

analyzerChain.append(
    LepJetFinder(
        jetCollection = lambda event: event.selectedJets,
        leptonCollection = lambda event: event.looseMuons,
    )
)
 
analyzerChain.append(
    InvariantSystem(
        inputCollection = lambda event: [event.tightMuons[0], event.lepJet[0]],
        outputName = "lepjet_muon"
    )
)


'''
analyzerChain.append(
    EventSkim(selection=lambda event: 
        event.dimuon_mass > 20 and event.dimuon_mass < 85,
    )
)
  
analyzerChain.append(
    EventSkim(selection=lambda event: 
        event.dimuon_deltaR > 1 and event.dimuon_deltaR < 5,
    )
)


if args.inputFiles[0].find("Heavy")>=0:
    storeVariables = [
        [lambda tree: tree.branch("llp","F"),lambda tree,event: tree.fillBranch("llp",Collection(event,"llpinfo").llp_mass)],
        [lambda tree: tree.branch("llp","F"),lambda tree,event: tree.fillBranch("llp",Collection(event,"llpinfo").llp_pt)],
    ]
  
    analyzerChain.append(EventInfo(storeVariables=storeVariables))
   
'''


analyzerChain.append(
    TaggerEvaluation(
        modelPath="PhysicsTools/NanoAODTools/data/nn/HNL.pb",
        inputCollections=[
            lambda event: event.selectedJets,
            lambda event: event.lepJet
        ],
        taggerName="llpdnnx_da",
    )
)


analyzerChain.append(
    JetTaggerResult(
        inputCollection = lambda event: event.lepJet,
        taggerName = "llpdnnx_da",
        logctauValues = range(-1, 3),
        outputName = "lepJet",
        predictionLabels = ["LLP"],
    )
)


analyzerChain.append(
    TaggerWorkingpoints(
        inputCollection = lambda event: event.selectedJets,
        taggerName = "llpdnnx_da",
        logctauValues = range(-1, 3),
        predictionLabels = ["LLP"],
        multiplicities = [0]
    )
)

analyzerChain.append(
    JetTruthFlags(inputCollection= lambda event: event.selectedJets,
        outputName="selectedJets"
    )
)

analyzerChain.append(
    JetTruthFlags(inputCollection= lambda event: event.lepJet,
        outputName="lepJet"
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
