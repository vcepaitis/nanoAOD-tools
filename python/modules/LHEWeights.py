
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class LHEWeights(Module):
    def __init__(self):
        self.couplings = 68
        
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for i in range(1, self.couplings):
            self.out.branch("LHEWeights_coupling_%i"%i,"F")
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
            
    def analyze(self, event):

        for i in range(1, self.couplings):
            self.out.fillBranch("LHEWeights_coupling_%i"%i, getattr(event, "LHEWeights_coupling_%i" %i))
        
        return True




