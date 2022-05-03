from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import numpy as np
import math

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
        self.correctionCoefficients = {2016: [-1.54794802,  1.14272696,  0.18702605, -0.28109537],
                                       2017: [-1.84007013,  1.13222034,  0.25396792, -0.13395763], 
                                       2018: [-1.85331292,  1.14988184,  0.2155245,  -0.19477323]
                                      }[self.globalOptions["year"]]
        self.correctionFunction = lambda p, x: np.piecewise(x, [x < p[0], x >= p[0]],
                     [lambda x:p[2]*x + p[1]-p[2]*p[0], lambda x:p[3]*x + p[1]-p[3]*p[0]])
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
        else:
            jet = jets[0]
            matchedCpfs = [cpf for cpf in cpfs if cpf.jetIdx == jet._index and getattr(cpf, matchingVariable, True)]
            for cpf in matchedCpfs:
                cpf.pt = jet.ptRaw*cpf.ptrel
                if cpf.matchedMuon == 1 and cpfMuon == 0:
                  cpfMuon = 1
                  cpfMuonpt = cpf.pt
                  cpfMuondxy= cpf.trackSip2dVal
                  cpfMuondxysig = cpf.trackSip2dSig
                elif cpf.matchedElectron == True and cpfElectron == 0 :
                   cpfElectron = 1
                   cpfElectronpt = cpf.pt
                   cpfElectrondxy = cpf.trackSip2dVal
                   cpfMuondxysig = cpf.trackSip2dSig
                  
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

            if self.storeWeights:
                matchedCpfs = [cpf for cpf in cpfs if cpf.jetIdx == jet._index]
                if len(matchedCpfs) > 0:
                    # Apply the corrections based on leading three constituents
                    matchedCpfs = matchedCpfs[:3]
                    weightsPerTrack = np.asarray(map(lambda x: self.correctionFunction(self.correctionCoefficients, math.log10(max(1e-3, x.trackSip2dVal))), matchedCpfs))
                    ptsPerTrack = np.asarray(map(lambda x: x.ptrel, matchedCpfs))
                    averageWeightUp = np.sum(weightsPerTrack*ptsPerTrack/np.sum(ptsPerTrack))

                else:
                    averageWeightUp = 1.

        if self.storeWeights:
            averageWeightDown = 2. - averageWeightUp
            self.out.fillBranch(self.outputName+"_"+self.svType+"_nominal", 1.)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_up", averageWeightUp)
            self.out.fillBranch(self.outputName+"_"+self.svType+"_down", averageWeightDown)

        return True
