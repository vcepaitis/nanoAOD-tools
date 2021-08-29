import os
import sys
import math
import argparse
import random
import ROOT
import json
import numpy as np

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

print "isData:",args.isData
print "inputs:",len(args.inputFiles)


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
    "isData": args.isData,
    "year": year
}

isMC = not args.isData

minMuonPt = {2016: 26., 2017: 29., 2018: 26.}
minElectronPt = {2016: 29., 2017: 34., 2018: 34.}

if isMC:
    jecTags = {2016: 'Summer16_07Aug2017_V11_MC',
               2017: 'Fall17_17Nov2017_V32_MC',
               2018: 'Autumn18_V19_MC'
               }

    jerTags = {2016: 'Summer16_25nsV1_MC',
               2017: 'Fall17_V3_MC',
               2018: 'Autumn18_V7_MC'
               }

if args.isData:
    jecTags = {2016: 'Summer16_07Aug2017All_V11_DATA',
               2017: 'Fall17_17Nov2017_V32_DATA',
               2018: 'Autumn18_V19_DATA'
               }

met_variable = {2016: lambda event: Object(event, "MET"),
                2017: lambda event: Object(event, "METFixEE2017"),
                2018: lambda event: Object(event, "MET")
                }


jesUncertaintyFile = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Autumn18_V19_MC_Uncertainty_AK4PFchs.txt"
}
jerResolutionFile = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Fall17_V3_MC_PtResolution_AK4PFchs.txt",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Autumn18_V7_MC_PtResolution_AK4PFchs.txt"
}

jerSFUncertaintyFile = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer16_25nsV1_MC_SF_AK4PFchs.txt",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Fall17_V3_MC_SF_AK4PFchs.txt",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Autumn18_V7_MC_SF_AK4PFchs.txt"
}


leptonSelection = [
    MuonSelection(
        outputName="tightMuons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all', 'charge'],
        storeWeights=True,
        muonMinPt=minMuonPt[globalOptions["year"]],
        muonMaxDxy=0.01,
        muonMaxDz=0.05,
        muonID=MuonSelection.TIGHT,
        muonIso=MuonSelection.TIGHT,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        outputName="tightElectrons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi','pfRelIso03_all', 'charge'],
        electronMinPt=minElectronPt[globalOptions["year"]],
        electronID="Iso_WP90",
        storeWeights=True,
        electronIPCuts=True,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.ntightMuons + event.ntightElectrons > 0),
    SingleMuonTriggerSelection(
        inputCollection=lambda event: event.tightMuons,
        outputName="IsoMuTrigger",
        storeWeights=True,
        globalOptions=globalOptions
    ),
    SingleElectronTriggerSelection(
        inputCollection=lambda event: event.tightElectrons,
        outputName="IsoElectronTrigger",
        storeWeights=False,
        globalOptions=globalOptions
    ),
    MuonSelection(
        inputCollection=lambda event: event.tightMuons_unselected,
        outputName="looseMuons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all',
                         'looseId', 'mediumId', 'tightId', 'charge'],
        muonMinPt=3.,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.NONE,
        globalOptions=globalOptions
    ),
    MuonSelection(
        inputCollection=lambda event: event.looseMuons,
        outputName="looseIsoMuons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all',
                         'looseId', 'mediumId', 'tightId', 'charge'],
        muonMinPt=3.,
        storeWeights=False,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.TIGHT,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        inputCollection=lambda event: event.tightElectrons_unselected,
        outputName="looseElectrons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi','pfRelIso03_all', 'charge',
                         'mvaFall17V2noIso_WP80', 'mvaFall17V2noIso_WP90', 'mvaFall17V2noIso_WPL',
                         'mvaFall17V2Iso_WP80', 'mvaFall17V2Iso_WP90', 'mvaFall17V2Iso_WPL',
                         'cutBased'],

        electronMinPt=5.,
        electronID="Custom",
        globalOptions=globalOptions
    ),
    ElectronSelection(
        inputCollection=lambda event: event.looseElectrons,
        outputName="looseIsoElectrons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi','pfRelIso03_all', 'charge',
                         'mvaFall17V2noIso_WP80', 'mvaFall17V2noIso_WP90', 'mvaFall17V2noIso_WPL',
                         'mvaFall17V2Iso_WP80', 'mvaFall17V2Iso_WP90', 'mvaFall17V2Iso_WPL',
                         'cutBased'],
        electronMinPt=5.,
        storeWeights=False,
        electronID="CustomIso",
        globalOptions=globalOptions
    ),
    LeptonCollecting(
        tightMuonsCollection=lambda event:event.tightMuons,
        tightElectronsCollection=lambda event:event.tightElectrons,
        looseMuonCollection=lambda event:event.looseMuons,
        looseElectronCollection=lambda event:event.looseElectrons,
        outputName = "Leptons"
    ),

    EventSkim(selection=lambda event: (event.IsoMuTrigger_flag + event.IsoElectronTrigger_flag) > 0),
    EventSkim(selection=lambda event: event.leadingLeptons[0].isTriggerMatched>0),
    EventSkim(selection=lambda event: event.nsubleadingLeptons==1),
]

analyzerChain = []

analyzerChain.extend(leptonSelection)


analyzerChain.append(
    InvariantSystem(
        inputCollection= lambda event:
            sorted(event.tightMuons+event.looseMuons+event.tightElectrons+event.looseElectrons,key=lambda x: -x.pt)[:2],
        outputName="dilepton"
    )
)

analyzerChain.append(EventSkim(selection=lambda event: (event.dilepton_mass > 80.) and (event.dilepton_pt > 40.)))


analyzerChain.append(
     MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
     )
)


if isMC:
    analyzerChain.append(
        JetMetUncertainties(
            metInput=met_variable[year],
            rhoInput = lambda event: event.fixedGridRhoFastjetAll,
            jetCollection = lambda event: Collection(event,"Jet"),
            lowPtJetCollection = lambda event: Collection(event,"CorrT1METJet"),
            genJetCollection = lambda event: Collection(event,"GenJet"),
            muonCollection = lambda event: Collection(event,"Muon"),
            electronCollection = lambda event: Collection(event,"Electron"),
            jesUncertaintyFile=jesUncertaintyFile[year],
            jerResolutionFileName=jerResolutionFile[year],
            jerSFUncertaintyFileName=jerSFUncertaintyFile[year],
            propagateJER = False,
            jetKeys = ['rawFactor','jetId', 'nConstituents'],
        )
    )

    for systName, jetCollection in [
        ("nominal", lambda event: event.jets_nominal),
    ]:

        analyzerChain.append(
            JetSelection(
                inputCollection=jetCollection,
                leptonCollectionDRCleaning=lambda event: event.tightMuons+event.tightElectrons,#+event.looseIsoMuons+event.looseIsoElectrons,
                leptonCollectionP4Subraction=lambda event: [], #event.looseMuons+event.looseElectrons,
                jetMinPt=15.,
                jetMaxEta=2.4,
                jetId=JetSelection.TIGHT,
                flagDA = True,
                storeKinematics=['pt', 'eta', 'phi', 'minDeltaRSubtraction', 'ptLepton', 'ptOriginal', 'ptSubtracted'],
                outputName="selectedJets_"+systName,
                globalOptions=globalOptions
            )
        )

        
    analyzerChain.append(
        EventSkim(
            selection=lambda event: event.nselectedJets_nominal > 0 and event.nselectedJets_nominal < 5
        )
    )

else:
    analyzerChain.append(
        JetSelection(
            inputCollection=lambda event: Collection(event, "Jet"),
            leptonCollectionDRCleaning=lambda event: event.tightMuons+event.tightElectrons,#+event.looseIsoMuons+event.looseIsoElectrons,
            leptonCollectionP4Subraction=lambda event: [], #event.looseMuons+event.looseElectrons,
            jetMinPt=15.,
            jetMaxEta=2.4,
            jetId=JetSelection.TIGHT,
            flagDA = True,
            storeKinematics=['pt', 'eta', 'phi', 'minDeltaRSubtraction', 'ptLepton', 'ptOriginal', 'ptSubtracted'],
            outputName="selectedJets_nominal",
            globalOptions=globalOptions
        )
    )
    

    analyzerChain.append(
        EventSkim(
            selection=lambda event: event.nselectedJets_nominal > 0 and event.nselectedJets_nominal < 5
        )
    )
    


analyzerChain.append(
    PileupWeight(
        outputName="puweight",
        globalOptions=globalOptions
    )
)

yieldFile = {
    2016: "/vols/cms/LLP/yields_201117/2016/eventyields.json",
    2017: "/vols/cms/LLP/yields_201117/2017/eventyields.json",
    2018: "/vols/cms/LLP/yields_201117/2018/eventyields.json"
}

lumi = {
    2016: 35.92*1e3,
    2017: 41.53*1e3,
    2018: 59.68*1e3,
}

analyzerChain.append(
    XsecWeight(
        yields = json.load(open(yieldFile[year])),
        lumi = lumi[year],
        globalOptions=globalOptions
    )
)

p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    cut="(nJet>0)&&((nElectron+nMuon)>0)",
    modules=analyzerChain,
    maxEvents=-1,
    friend=False
)

p.run()
