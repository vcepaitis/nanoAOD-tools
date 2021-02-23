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
parser.add_argument('--input', dest='inputFiles', action='append', default=[])
parser.add_argument('--year', dest='year', action='store', type=int, default=2016)
parser.add_argument('output', nargs=1)

args = parser.parse_args()

print "inputs:",len(args.inputFiles)
print "year:",args.year

featureDictFile = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/feature_dict.py"
modelPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2016_NominalNetwork_ref_201117.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2017_NominalNetwork_ref_201117.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2018_NominalNetwork_ref_201117.pb"
}

analyzerChain = [
    TaggerEvaluationCache(
        modelPath=modelPath[args.year],
        featureDictFile=featureDictFile,
        taggerName="llpdnnx",
        evalValues = np.linspace(-1.9,1.9,5*4),
    )
]


p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    modules=analyzerChain,
    maxEvents=-1,
    friend=True
)

p.run()
