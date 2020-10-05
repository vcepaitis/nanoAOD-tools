import os
import sys
import math
import argparse
import random
import ROOT
import numpy as np

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

for inputFile in args.inputFiles:
    if "-2016" in inputFile or "Run2016" in inputFile:
        year = 2016
    elif "-2017" in inputFile or "Run2017" in inputFile:
        year = 2017
    elif "-2018" in inputFile or "Run2018" in inputFile:
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

chain = []

chain.append(LeptonGenEfficiency())

storeVariables = [
]

for coupling in range(1,68):
    storeVariables.append([
        lambda tree, coupling=coupling: tree.branch('LHEWeights_coupling_%i'%coupling,'F'),
        lambda tree, event, coupling=coupling: tree.fillBranch('LHEWeights_coupling_%i'%coupling,getattr(event,"LHEWeights_coupling_%i"%coupling)),
    ])

chain.append(EventInfo(storeVariables=storeVariables))

p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    modules=chain,
    maxEvents=-1,
    friend=True
)

p.run()
