import os
import sys
import math
import json
import argparse
import random
import ROOT
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor \
    import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel \
    import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.modules import *

parser = argparse.ArgumentParser()

parser.add_argument('--isData', dest='isData',
                    action='store_true', default=False)
parser.add_argument('--year', dest='year',
                    action='store', type=int, default=2016)
parser.add_argument('--input', dest='inputFiles', action='append', default=[])
parser.add_argument('output', nargs=1)

args = parser.parse_args()

print "isData:", args.isData
print "inputs:", len(args.inputFiles)

for inputFile in args.inputFiles:
    if "-2016" in inputFile:
        year = 2016
    elif "-2017" in inputFile:
        year = 2017
    elif "-2018" in inputFile:
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
    print " - ", inputFile, ", events=", tree.GetEntries()

print "year:", year
print "output directory:", args.output[0]

globalOptions = {
    "isData": args.isData,
    "year": year
}

isMC = not args.isData

minMuonPt = {2016: 25., 2017: 28., 2018: 25.}


muonSelection = [
    EventSkim(selection=lambda event: event.nTrigObj > 0),
    MuonSelection(
        outputName="tightMuon",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all', 'charge'],
        storeWeights=True,
        muonMinPt=minMuonPt[globalOptions["year"]],
        triggerMatch=True,
        muonID=MuonSelection.TIGHT,
        muonIso=MuonSelection.TIGHT,
        selectLeadingOnly=True,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.ntightMuon == 1),
    MuonSelection(
        inputCollection=lambda event: event.tightMuon_unselected,
        outputName="looseMuons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all',
                         'tightId', 'charge'],
        storeWeights=True,
        muonMinPt=5.,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.NONE,
        globalOptions=globalOptions
    ),

    SingleMuonTriggerSelection(
        inputCollection=lambda event: event.tightMuon,
        outputName="IsoMuTrigger",
        storeWeights=True,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.IsoMuTrigger_flag == 1),
    EventSkim(selection=lambda event: event.nlooseMuons > 0)
]

analyzerChain = []

analyzerChain.extend(muonSelection)


analyzerChain.append(
    InvariantSystem(
        inputCollection=lambda event: [event.looseMuons[0],
                                       event.tightMuon[0]],
        outputName="dimuon"
    )
)

analyzerChain.append(
     MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
     )
)


analyzerChain.append(
    JetSelection(
        leptonCollection=lambda event: event.tightMuon,
        outputName="selectedJets_nominal",
        storeKinematics=['pt', 'eta'],
        globalOptions=globalOptions
    )
)

analyzerChain.append(
    EventSkim(
        selection=lambda event: len(event.selectedJets_nominal) > 0
    )
)

analyzerChain.append(
    LepJetFinder(
        jetCollection=lambda event: event.selectedJets_nominal,
        leptonCollection=lambda event: event.looseMuons,
        outputName="lepJet_nominal"
    )
)


analyzerChain.append(
    EventObservables(
        jetCollection=lambda event: event.selectedJets_nominal,
        leptonCollection=lambda event: event.tightMuon[0],
        outputName="EventObservables_nominal"
    )
)
analyzerChain.append(
    EventSkim(
        selection=lambda event: len(event.lepJet_nominal) > 0
    )
)

analyzerChain.append(
    EventSkim(
        selection=lambda event: event.dimuon_mass < 80.
    )
)
'''
analyzerChain.append(
    EventSkim(
        selection=lambda event: event.MET_pt > 100.
    )
)
'''
p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    modules=analyzerChain,
    maxEvents=-1,
    friend=False
)

p.run()
