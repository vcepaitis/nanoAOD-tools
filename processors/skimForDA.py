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


yields = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/da/"+str(year)+"/eventyields.json"




leptonSelection = [
    EventSkim(selection=lambda event: event.nTrigObj > 0),
    MuonSelection(
        outputName="tightMuon",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all', 'charge'],
        storeWeights=True,
        muonMinPt=minMuonPt[globalOptions["year"]],
        muonMaxDxy=0.01,
        muonMaxDz=0.05,
        triggerMatch=True,
        muonID=MuonSelection.TIGHT,
        muonIso=MuonSelection.TIGHT,
        selectLeadingOnly=True,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        outputName="tightElectron",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi','pfRelIso03_all', 'charge'],
        electronMinPt=minElectronPt[globalOptions["year"]],
        electronID="Iso_WP90",
        storeWeights=True,
        triggerMatch=True,
        electronIPCuts=True,
        selectLeadingOnly=True,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.ntightMuon + event.ntightElectron > 0),
    SingleMuonTriggerSelection(
        inputCollection=lambda event: event.tightMuon,
        outputName="IsoMuTrigger",
        storeWeights=True,
        globalOptions=globalOptions
    ),
    SingleElectronTriggerSelection(
        inputCollection=lambda event: event.tightElectron,
        outputName="IsoElectronTrigger",
        storeWeights=False,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.IsoMuTrigger_flag + event.IsoElectronTrigger_flag > 0),
    MuonSelection(
        inputCollection=lambda event: event.tightMuon_unselected,
        outputName="looseMuons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all',
                         'looseId', 'mediumId', 'tightId', 'charge'],
        muonMinPt=5.,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.NONE,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        inputCollection=lambda event: event.tightElectron_unselected,
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


    LeptonCollecting(
        tightMuonCollection=lambda event:event.tightMuon,
        tightElectronCollection=lambda event:event.tightElectron,
        looseMuonCollection=lambda event:event.looseMuons,
        looseElectronCollection=lambda event:event.looseElectrons,
        outputName = "Leptons"
    ),
    EventSkim(selection=lambda event: event.isTriggered),
    EventSkim(selection=lambda event: event.nsubleadingLeptons<2),


]

analyzerChain = []

analyzerChain.extend(leptonSelection)


analyzerChain.append(
    InvariantSystem(
        inputCollection= lambda event:
            sorted(event.tightMuon+event.looseMuons+event.tightElectron+event.looseElectrons,key=lambda x: -x.pt)[:2],
        outputName="dilepton"
    )
)

analyzerChain.append(EventSkim(selection=lambda event: event.dilepton_mass > 91.1876-10. and event.dilepton_mass < 91.1876+10.))


analyzerChain.append(
     MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
     )
)

analyzerChain.append(
    JetSelection(
        inputCollection=lambda event: Collection(event, "Jet"),
        leptonCollectionDRCleaning=lambda event: event.leadingLeptons,
        #leptonCollectionP4Subraction=lambda event: event.subleadingLeptons,
        jetMinPt=15.,
        jetMaxEta=2.399, #TODO: change to 2.4
        jetMinNConstituents=3,
        jetId=JetSelection.TIGHT,
        outputName="Jet",
        globalOptions=globalOptions,
        flagDA=True,
    )
)

analyzerChain.append(
    EventObservables(
        metInput=met_variable[year],
        globalOptions=globalOptions,
        outputName="EventObservables_nominal"
    )
)


analyzerChain.append(EventSkim(selection=lambda event: event.EventObservables_nominal_met < 100.))


analyzerChain.append(
    PileupWeight(
        outputName="puweight",
        globalOptions=globalOptions
    )
)

analyzerChain.append(
    DataFlag(
        yields,
        globalOptions=globalOptions,
    )
)

p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    cut="(nJet<5)&&(nJet>0)&&((nElectron+nMuon)>0)",
    modules=analyzerChain,
    maxEvents=-1,
    friend=False
)

p.run()
