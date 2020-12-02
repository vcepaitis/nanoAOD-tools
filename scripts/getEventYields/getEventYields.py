import subprocess
import re
import time
import os
import ROOT
import json


file_path = "files_200311/nanoAOD_files_200311/2016"

txtFiles = sorted(os.listdir(file_path))

processDict = {}
pileupHists = {}

for txtFile in txtFiles:
    process = txtFile.split(".")[0]
    if process.find("SingleMu")>=0:
        print "skip ",process
        continue
    processDict[process] = 0.
    pileupHists[process] = ROOT.TH1F(process,"",101,0,100)
    pileupHists[process].SetDirectory(0)
    
    f = open (os.path.join(file_path, txtFile))
    print "reading ",txtFile, "...",
    for l in f:
        if len(l)>0:
            fileName = l.replace("\n","").replace("\r","")
            print "opening",fileName
            rootFile=None
            while (rootFile==None):
                rootFile = ROOT.TFile.Open(fileName)
                if not rootFile:
                    print "Cannot found file: ",rootFile, "-> retry"
                    continue
                break
            tree = rootFile.Get("Events")
            if not rootFile or not tree:
                print "Cannot found tree in file: ",rootFile
                continue
            h = ROOT.TH1F("pu","",101,0,100)
            tree.Project(h.GetName(),"Pileup_nTrueInt","genWeight")
            pileupHists[process].Add(h)
            processDict[process] += h.Integral()
            rootFile.Close()
            break
    print processDict[process]
     
with open('eventyields.json', 'w') as outfile:
    json.dump(processDict, outfile,ensure_ascii=True,indent=2,sort_keys=True)    

rootFile = ROOT.TFile("pileup.root","RECREATE")
for process in pileupHists.keys():
    pileupHists[process].SetDirectory(rootFile)
    pileupHists[process].SetName(process)
    pileupHists[process].Write()
rootFile.Close()


