import os
import sys
import math
import json
import argparse
import random
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.modules import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2       import * 
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetRecalib import *


parser = argparse.ArgumentParser()
parser.add_argument('--isData', dest='isData', action='store_true',default=False)
parser.add_argument('--year', dest='year', action='store',type=int, default=2016)
parser.add_argument('--input', dest='inputFiles', action='append',default=[])
parser.add_argument('output', nargs=1)

args = parser.parse_args()


print "isData:",args.isData
print "inputs:",len(args.inputFiles)
for inputFile in args.inputFiles:
    print "2018" in inputFile
    if "2016" in inputFile:
        year = 2016
    elif "2017" in inputFile:
        year = 2017
    elif "2018" in inputFile:
        year = 2018
    else:
        year = args.year
    rootFile = ROOT.TFile.Open(inputFile)
    if not rootFile:
        print "CRITICAL - file '"+inputFile+"' not found!"
        sys.exit(1)
    tree = rootFile.Get("Events")
    if not tree:
        print "CRITICAL - 'Events' tree not found in file '"+inputFile+"'!"
        sys.exit(1)
    print " - ",inputFile,", events=",tree.GetEntries()
    
print "year:",year
print "output directory:",args.output[0]

globalOptions = {
    "isData":args.isData,
    "year":year
}

isMC = not args.isData

minMuonPt = {2016: 25., 2017: 28., 2018: 25.}

muonSelection = [
    EventSkim(selection=lambda event: event.nTrigObj>0),
    MuonSelection(
        outputName="tightMuons",
        storeKinematics=['pt','eta', 'dxy', 'dxyErr', 'dz', 'dzErr', 'phi', 'pfRelIso04_all', 'looseId', 'tightId'],
        storeWeights=True,
        muonMinPt = minMuonPt[globalOptions["year"]],
        triggerMatch = True,
        muonID = MuonSelection.TIGHT,
        muonIso = MuonSelection.TIGHT,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.ntightMuons>0),
    MuonSelection(
        inputCollection = lambda event: [muon for muon in Collection(event, "Muon") if abs(muon.pt-event.tightMuons[0].pt)>1e-4],
        outputName="looseMuons",
        storeKinematics=['pt','eta', 'dxy', 'dxyErr', 'dz', 'dzErr', 'phi', 'pfRelIso04_all', 'looseId', 'tightId'],
        storeWeights=True,
        muonMinPt = 5.,
        muonID = MuonSelection.LOOSE,
        muonIso = MuonSelection.NONE,
        globalOptions=globalOptions
    ),
    
    SingleMuonTriggerSelection(
        inputCollection=lambda event: event["tightMuons"],
        outputName="IsoMuTrigger",
        storeWeights=True,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.IsoMuTrigger_flag==1),
    EventSkim(selection=lambda event: event.nlooseMuons>0)
]

#analyzerChain = [jetRecalibration]
analyzerChain = []

analyzerChain.extend(muonSelection)

analyzerChain.append(
    JetSelection(
        leptonCollection=lambda event: [event.tightMuons[0]],
        globalOptions=globalOptions
    )
)

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


'''
analyzerChain.append(
        JetFeatures()
)
'''

analyzerChain.append(
    TaggerEvaluation(
        modelPath="PhysicsTools/NanoAODTools/data/nn/weight2016_75.pb",
        logctauValues = range(-1, 4),
        inputCollections=[
            lambda event: event.lepJet,
        ],
        taggerName="llpdnnx",
    )
)


analyzerChain.append(
    JetTaggerResult(
        inputCollection = lambda event: event.lepJet,
        taggerName = "llpdnnx",
        outputName = "lepJet",
        logctauValues = range(-1, 4),
    )
)


analyzerChain.append(
    JetTruthFlags(inputCollection= lambda event: event.lepJet,
        globalOptions=globalOptions,
        outputName="lepJet"
    )
)


analyzerChain.append( 
     MetFilter(
     	globalOptions=globalOptions,
     	outputName="MET_filter"
     	)
)
analyzerChain.append(
     EventObservables(
       jetCollection = lambda event: event.selectedJets,
       leptonCollection = lambda event: event.tightMuons[0] 
     )
)
'''
 
analyzerChain.append(
    XGBEvaluation(
        modelPath="PhysicsTools/NanoAODTools/data/nn/bdt.model",
    )
)
'''

storeVariables = [
    [lambda tree: tree.branch("MET_pt", "F"), lambda tree,event: tree.fillBranch("MET_pt", event.MET_pt)],
    [lambda tree: tree.branch("MET_phi", "F"), lambda tree,event: tree.fillBranch("MET_phi", event.MET_phi)],
    [lambda tree: tree.branch("MET_significance", "F"), lambda tree,event: tree.fillBranch("MET_significance", event.MET_significance)],
    [lambda tree: tree.branch("PV_npvs", "I"), lambda tree,event: tree.fillBranch("PV_npvs", event.PV_npvs)],
    [lambda tree: tree.branch("PV_npvsGood", "I"), lambda tree,event: tree.fillBranch("PV_npvsGood", event.PV_npvsGood)],
    [lambda tree: tree.branch("fixedGridRhoFastjetAll", "F"), lambda tree,event: tree.fillBranch("fixedGridRhoFastjetAll", event.fixedGridRhoFastjetAll)],
    #[lambda tree: tree.branch("bdt_score", "F"), lambda tree,event: tree.fillBranch("bdt_score", event.bdt_score)]
]

if not globalOptions["isData"]:
    storeVariables.append([lambda tree: tree.branch("genweight","F"),lambda tree,event: tree.fillBranch("genweight",event.Generator_weight)])


analyzerChain.append(EventInfo(storeVariables=storeVariables))

analyzerChain.append(
    PileupWeight(
        outputName ="puweight",
        globalOptions=globalOptions
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
