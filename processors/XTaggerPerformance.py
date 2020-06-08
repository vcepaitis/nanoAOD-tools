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
parser.add_argument('--input', dest='inputFiles', action='append', default=[])
parser.add_argument('output', nargs=1)

args = parser.parse_args()

print "inputs:",len(args.inputFiles)

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
    "year": year
}

isMC = True

# analyzerChain = [jetRecalibration]
analyzerChain = []

featureDictFile = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200311/feature_dict.py"
modelPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200311/weight2016_attention.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200311/weight2017_attention.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200311/weight2018_attention.pb"
}

analyzerChain.append(
     MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
     )
)

analyzerChain.append(
    JetSelection(
        leptonCollection=lambda event:None,
        outputName="selectedJets_nominal",
        storeKinematics=['pt', 'eta'],
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

analyzerChain.append(
    TaggerEvaluation(
        modelPath=modelPath[year],
        evalValues = np.linspace(-1, 3, num=9),
        featureDictFile = featureDictFile,
        inputCollections=[
            lambda event: event.selectedJets_nominal,
        ],
        taggerName="llpdnnx",
        integrateDisplacementOrder=-1
    )
)

'''
analyzerChain.append(
    JetTaggerIntegral(
        taggerName="llpdnnx",
        integrateDisplacementOrder=-1,
        inputCollection=lambda event: event.selectedJets_nominal,
        outputName="selectedJets_nominal",
    )
)
'''

analyzerChain.append(
    JetTaggerResult(
        inputCollection=lambda event:event.selectedJets_nominal,
        taggerName="llpdnnx",
        outputName="selectedJets_nominal",
    )
)

analyzerChain.append(
    EventSkim(
        selection=lambda event: event.nselectedJets_nominal > 0
    )
)

storeVariables = [[lambda tree: tree.branch("genweight", "F"),
                       lambda tree,
                       event: tree.fillBranch("genweight",
                       event.Generator_weight)]
                ]

analyzerChain.append(EventInfo(storeVariables=storeVariables))



p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    modules=analyzerChain,
    maxEvents=-1,
    friend=True
)

p.run()

