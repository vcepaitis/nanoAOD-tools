from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import numpy as np
import math
import utils

class TrackAndSVSelection(Module):
    def __init__(
        self,
        jetCollection=lambda event: Collection(event, "Jet"),
        svType="adapted",
        cpfCollection=lambda event: Collection(event, "cpf"),
        outputName="hnlJet_track_weight",
        globalOptions={"isData":False, "isSignal":False, "year": 2016},
        storeWeights="False",

    ):
        self.jetCollection = jetCollection
        if svType not in ["regular", "adapted"]:
            raise ValueError("Wrong sv type")
        self.svType = svType
        self.cpfCollection = cpfCollection
        self.outputName = outputName
        self.globalOptions = globalOptions
        if self.globalOptions["year"]== 2016 :   
           self.sf = utils.getHistCanvas("PhysicsTools/NanoAODTools/data/track/track_sf_2016.root", "c1" , "ratio")

        if self.globalOptions["year"]== 2017 :   
           self.sf = utils.getHistCanvas("PhysicsTools/NanoAODTools/data/track/track_sf_2017.root", "c1" , "ratio")

        if self.globalOptions["year"]== 2018 :   
           self.sf = utils.getHistCanvas("PhysicsTools/NanoAODTools/data/track/track_sf_2018.root", "c1" , "ratio")
      
        self.storeWeights = storeWeights

    def beginJob(self):
        pass
    def endJob(self):
        pass
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch(self.outputName+"_"+self.svType+"_muon", "I")
        self.out.branch(self.outputName+"_"+self.svType+"_muon_dxy", "F")
        self.out.branch(self.outputName+"_"+self.svType+"_muon_pt", "F")
        self.out.branch(self.outputName+"_"+self.svType+"_muon_dxysig", "F")
        
        self.out.branch(self.outputName+"_"+self.svType+"_electron", "I")
        self.out.branch(self.outputName+"_"+self.svType+"_electron_pt", "F")
        self.out.branch(self.outputName+"_"+self.svType+"_electron_dxy", "F")
        self.out.branch(self.outputName+"_"+self.svType+"_electron_dxysig", "F")
        


        self.out.branch("n"+self.outputName+"_"+self.svType, "I")
        self.out.branch(self.outputName+"_"+self.svType+"_pt", "F", lenVar="nhnlJet_svMatchedTracks_"+self.svType)
        self.out.branch(self.outputName+"_"+self.svType+"_dxy", "F", lenVar="nhnlJet_svMatchedTracks_"+self.svType+"_"+self.outputName)
        self.out.branch(self.outputName+"_"+self.svType+"_dxysig", "F", lenVar="nhnlJet_svMatchedTracks_"+self.svType+"_"+self.outputName)
      
        if self.storeWeights:        
            self.out.branch("hnlJet_track_weight_"+self.svType+"_nominal", "F")
            self.out.branch("hnlJet_track_weight_"+self.svType+"_up", "F")
            self.out.branch("hnlJet_track_weight_"+self.svType+"_down", "F")


    def analyze(self, event):
          
        jets = self.jetCollection(event)
        cpfs = self.cpfCollection(event)

        # Iterate over all jets and find the matching tracks
        matchingVariable = "matchedSV_adapted" if self.svType == "adapted" else "matchedSV"
        cpfsMatchedToJetAndSV = []
        cpfMuon = 0
        cpfMuonpt = 0
        cpfMuondxy = 0
        cpfMuondxysig = 0
        cpfElectron = 0
        cpfElectronpt = 0
        cpfElectrondxy = 0
        cpfElectrondxysig = 0
        
        if len(jets) == 0:
            averageWeightUp = 1.
            ntracks = 0
            trackpt = 0
            trackdxy = 0
            trackdxysig = 0
        else:
            jet = jets[0]
            ntracks = len(cpfs)
            matchedCpfs = [cpf for cpf in cpfs if cpf.jetIdx == jet._index and getattr(cpf, matchingVariable, True)]
            if len(matchedCpfs) == 0:
                ntracks = 0
                trackpt = 0
                trackdxy = 0
                trackdxysig = 0
            else :
               trackpt = matchedCpfs[0].ptrel*jet.ptRaw
               trackdxy = matchedCpfs[0].trackSip2dVal
               trackdxysig = matchedCpfs[0].trackSip2dSig
            
            for cpf in matchedCpfs:
                cpf.pt = jet.ptRaw*cpf.ptrel
                if cpf.matchedMuon == True and cpfMuon == 0:
                  cpfMuon = 1
                  cpfMuonpt = cpf.pt
                  cpfMuondxy= cpf.trackSip2dVal
                  cpfMuondxysig = cpf.trackSip2dSig
                elif cpf.matchedElectron == True and cpfElectron == 0 :
                   cpfElectron = 1
                   cpfElectronpt = cpf.pt
                   cpfElectrondxy = cpf.trackSip2dVal
                   cpfElectrondxysig = cpf.trackSip2dSig
                  
            cpfsMatchedToJetAndSV.extend(matchedCpfs)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_muon", cpfMuon)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_muon_pt", cpfMuonpt)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_muon_dxy", cpfMuondxy)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_muon_dxysig", cpfMuondxysig)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_electron", cpfElectron)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_electron_pt", cpfElectronpt)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_electron_dxy", cpfElectrondxy)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_electron_dxysig", cpfElectrondxysig)
            

            
            self.out.fillBranch("n"+self.outputName+"_"+self.svType, len(cpfsMatchedToJetAndSV))
            self.out.fillBranch(self.outputName+"_"+self.svType+"_pt", map(lambda cpf: cpf.pt, cpfsMatchedToJetAndSV))
            self.out.fillBranch(self.outputName+"_"+self.svType+"_dxy", map(lambda cpf: cpf.trackSip2dVal, cpfsMatchedToJetAndSV))
            self.out.fillBranch(self.outputName+"_"+self.svType+"_dxysig", map(lambda cpf: cpf.trackSip2dSig, cpfsMatchedToJetAndSV))
            
            if self.storeWeights:
                matchedCpfs = [cpf for cpf in cpfs if cpf.jetIdx == jet._index]
                weightsPerTrack = []
                if len(matchedCpfs) > 0:
                    # Apply the corrections based on leading three constituents
                    matchedCpfs = matchedCpfs[:3]
                    #weightsPerTrack = np.asarray(map(lambda x: self.correctionFunction(self.correctionCoefficients, math.log10(max(1e-3, x.trackSip2dVal))), matchedCpfs)) 
                    for cpf in matchedCpfs : 
                      print "it arrives where you want and the dxy value is " , cpf.trackSip2dVal
                      binNumber = self.sf.FindBin(abs(cpf.trackSip2dVal))
 		      print "the given number of a given dxy is " , binNumber
                      weightsPerTrack.append(self.sf.GetBinContent(binNumber))

                    ptsPerTrack = np.asarray(map(lambda x: x.ptrel, matchedCpfs))
		    print "the pts per track are:  " , ptsPerTrack
                    averageWeight = np.sum(weightsPerTrack*ptsPerTrack/np.sum(ptsPerTrack))
                    print "average Weight gives,  ", averageWeight
                else:
                    averageWeight = 1.
        #### needs to be updated.
        if self.storeWeights:
            averageWeightDown = 2. - averageWeightUp
            self.out.fillBranch(self.outputName+"_"+self.svType+"_nominal",averageWeight )
            self.out.fillBranch(self.outputName+"_"+self.svType+"_up", averageWeightUp)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_down", averageWeightDown)

        return True
