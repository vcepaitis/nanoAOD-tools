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
        outputName="nominal",
        globalOptions={"isData":False, "isSignal":False, "year": 2016},
        storeWeights="False",

    ):
        self.jetCollection = jetCollection
        if svType not in ["nominal", "adapted"]:
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
        self.out.branch("nhnlJets_svMatchedTracks_"+self.svType+"_"+self.outputName, "I")
        self.out.branch("hnlJets_svMatchedTracks_"+self.svType+"_"+self.outputName+"_pt", "F", lenVar="nhnlJets_svMatchedTracks_"+self.svType+"_"+self.outputName)
        self.out.branch("hnlJets_svMatchedTracks_"+self.svType+"_"+self.outputName+"_dxy", "F", lenVar="nhnlJets_svMatchedTracks_"+self.svType+"_"+self.outputName)  

        if self.storeWeights:        
            for n in range(1, 6)+[25]:   
                self.out.branch("hnlJet_track_weight_"+self.svType+"_"+self.outputName, "F")

    def analyze(self, event):
          
        jets = self.jetCollection(event)
        cpfs = self.cpfCollection(event)

        # Iterate over all jets and find the matching tracks
        cpfsMatchedToJetsAndToSVs = []
        matchingVariable = "matchedSV_adapted" if self.svType == "adapted" else "matchedSV"
        for jet in jets:
            matchedCpfs = [cpf for cpf in cpfs if cpf.jetIdx == jet._index and getattr(cpf, matchingVariable, True)]
            for cpf in matchedCpfs:
                cpf.pt = jet.ptRaw*cpf.ptrel # Almost exactly the same as using global jet
            cpfsMatchedToJetsAndToSVs.extend(matchedCpfs)
        self.out.fillBranch("nhnlJets_svMatchedTracks_"+self.svType+"_"+self.outputName, len(cpfsMatchedToJetsAndToSVs))
        self.out.fillBranch("hnlJets_svMatchedTracks_"+self.svType+"_"+self.outputName+"_pt", map(lambda cpf: cpf.pt, cpfsMatchedToJetsAndToSVs))
        self.out.fillBranch("hnlJets_svMatchedTracks_"+self.svType+"_"+self.outputName+"_dxy", map(lambda cpf: cpf.trackSip2dVal, cpfsMatchedToJetsAndToSVs))

        if self.storeWeights:
            weightsPerJet = []
            for jet in jets:
                matchedCpfs = [cpf for cpf in cpfs if cpf.jetIdx == jet._index]
                if len(matchedCpfs) > 0:
                    # Apply the corrections based on leading three constituents
                    matchedCpfs = matchedCpfs[:3]
                    weightsPerTrack = np.asarray(map(lambda x: self.correctionFunction(self.correctionCoefficients, math.log10(max(1e-3, x.trackSip2dVal))), matchedCpfs))
                    ptsPerTrack = np.asarray(map(lambda x: x.ptrel, matchedCpfs))
                    averageWeight = np.sum(weightsPerTrack*ptsPerTrack/np.sum(ptsPerTrack))
                    weightsPerJet.append(averageWeight)
                else:
                    weightsPerJet.append(1.)
            
            weight = 1.
            for w in weightsPerJet:
                weight *= w
            #print(w)

            self.out.fillBranch("hnlJet_track_weight_"+self.svType+"_"+self.outputName, weight)

        return True
