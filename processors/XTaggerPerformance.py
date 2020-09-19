import os
import sys
import math
import argparse
import random
import numpy as np
import ROOT

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor \
    import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel \
    import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.modules import *

parser = argparse.ArgumentParser()

parser.add_argument('--year', dest='year',
                    action='store', type=int, default=2016)
parser.add_argument('--isSignal', dest='isSignal',
                    action='store_true', default=False)
parser.add_argument('--input', dest='inputFiles', action='append', default=[])
parser.add_argument('output', nargs=1)

args = parser.parse_args()

print "inputs:",len(args.inputFiles)
print "isSignal:",args.isSignal

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
    "isData": False,
    "isSignal": args.isSignal,
    "year": year
}

modelPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weight2016.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weight2017.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weight2018.pb"
}

modelGunPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weightGun2016.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weightGun2017.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weightGun2018.pb"
}



isMC = True

# analyzerChain = [jetRecalibration]
analyzerChain = []

analyzerChain.append(
    JetSelection(
        jetMinPt=15.,
        jetMaxEta=2.399, #TODO: change to 2.4
        jetMinNConstituents=3,
        jetId=JetSelection.LOOSE,
        outputName="selectedJets_nominal",
        globalOptions=globalOptions
    )
)

analyzerChain.append(
    JetTruthFlags(
        inputCollection= lambda event: event.selectedJets_nominal,
        outputName="selectedJets_nominal",
        globalOptions=globalOptions
    )
)

featureDictFile = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/feature_dict.py"

analyzerChain.append(
    TaggerEvaluationProfiled(
        modelPath=modelPath[year],
        featureDictFile=featureDictFile,
        inputCollections=[
            lambda event: event.selectedJets_nominal,
        ],
        taggerName="llpdnnx",
        globalOptions=globalOptions,
        evalValues = np.linspace(-3,2,5*5+1),
    )
)

analyzerChain.append(
    TaggerEvaluationProfiled(
        modelPath=modelGunPath[year],
        featureDictFile=featureDictFile,
        inputCollections=[
            lambda event: event.selectedJets_nominal,
        ],
        taggerName="llpdnnx_gun",
        globalOptions=globalOptions,
        evalValues = np.linspace(-3,2,5*5+1),
    )
)

analyzerChain.append(
    JetTaggerResult(
        taggerName="llpdnnx",
    )
)

analyzerChain.append(
    JetTaggerResult(
        taggerName="llpdnnx_gun",
    )
)

analyzerChain.append(
    EventSkim(
        selection=lambda event: event.nselectedJets_nominal > 0
    )
)

p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    modules=analyzerChain,
    maxEvents=30000,
    cut="(nJet>0)",
    friend=True
)

p.run()

