import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel \
    import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi


class EventObservables(Module):

    def __init__(
        self,
        lepton1Object = lambda event: event.leadingLeptons[0],
        lepton2Object = None,
        jetCollection=lambda event: Collection(event, "Jet"),
        metInput=lambda event: Object(event, "MET"),
        outputName="EventObservables",
        globalOptions={"isData": False}
    ):
        self.lepton1Object = lepton1Object
        self.lepton2Object = lepton2Object
        self.jetCollection = jetCollection
        self.metInput = metInput
        self.outputName = outputName
        self.globalOptions = globalOptions

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self.outputName+"_met", "F")
        self.out.branch(self.outputName+"_met_phi", "F")
        
        self.out.branch(self.outputName+"_mtw", "F")
        self.out.branch(self.outputName+"_met_l1_dphi", "F")
        
        self.out.branch(self.outputName+"_eventShape_isotropy", "F")
        self.out.branch(self.outputName+"_eventShape_circularity", "F")
        self.out.branch(self.outputName+"_eventShape_sphericity", "F")
        self.out.branch(self.outputName+"_eventShape_aplanarity", "F")
        self.out.branch(self.outputName+"_eventShape_C", "F")
        
        self.out.branch(self.outputName+"_ht", "F")
        self.out.branch(self.outputName+"_hmass", "F")
        self.out.branch(self.outputName+"_mht", "F")
        self.out.branch(self.outputName+"_mht_met_dphi", "F")
        self.out.branch(self.outputName+"_minPhiStar", "F")
        self.out.branch(self.outputName+"_mht_l1_dphi", "F")
        self.out.branch(self.outputName+"_mht_l1_ptRatio", "F")
        

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        lepton1 = self.lepton1Object(event)
        lepton2 = None if self.lepton2Object==None else self.lepton2Object(event)
    
        jets = self.jetCollection(event)

        met = self.metInput(event)
        self.out.fillBranch(self.outputName+"_met", met.pt)
        self.out.fillBranch(self.outputName+"_met_phi", met.phi)
        
        dPhi = deltaPhi(lepton1,met)
        mtw = math.sqrt(2*lepton1.pt*met.pt*(1-math.cos(dPhi)))

        self.out.fillBranch(self.outputName+"_mtw", mtw)
        self.out.fillBranch(self.outputName+"_met_l1_dphi", dPhi)
        
        eventShapes = ROOT.EventShapes()
        eventShapes.addObject(lepton1.pt, lepton1.eta, lepton1.phi, 0.0)
        if lepton2!=None:
            eventShapes.addObject(lepton2.pt, lepton2.eta, lepton2.phi, 0.0)

        # ht and mht
        vectorSum = ROOT.TLorentzVector()
        scalarPtSum = 0.0

        for jet in jets:
            eventShapes.addObject(
                jet.p4Subtracted.Pt(), 
                jet.p4Subtracted.Eta(), 
                jet.p4Subtracted.Phi(), 
                0.0
            )
            vectorSum += jet.p4Subtracted
            scalarPtSum += jet.ptSubtracted
            
        self.out.fillBranch(self.outputName+"_eventShape_isotropy", eventShapes.isotropy())
        self.out.fillBranch(self.outputName+"_eventShape_circularity", eventShapes.circularity())
        self.out.fillBranch(self.outputName+"_eventShape_sphericity", eventShapes.sphericity())
        self.out.fillBranch(self.outputName+"_eventShape_aplanarity", eventShapes.aplanarity())
        self.out.fillBranch(self.outputName+"_eventShape_C", eventShapes.C())
               
        mht_met_dphi = math.fabs(deltaPhi(vectorSum.Phi(),met.phi))

        # minPhiStar
        minPhiStar = math.pi
        for jet in jets:
            negSum = -(vectorSum-jet.p4Subtracted)
            minPhiStar = min(minPhiStar, math.fabs(deltaPhi(negSum.Phi(), jet.phi)))

        self.out.fillBranch(self.outputName+"_ht", scalarPtSum)
        self.out.fillBranch(self.outputName+"_hmass", vectorSum.M())
        self.out.fillBranch(self.outputName+"_mht", vectorSum.Pt())
        self.out.fillBranch(self.outputName+"_mht_met_dphi", mht_met_dphi)
        self.out.fillBranch(self.outputName+"_minPhiStar", minPhiStar)
        self.out.fillBranch(self.outputName+"_mht_l1_dphi", math.fabs(deltaPhi(lepton1.phi,vectorSum.Phi())))
        self.out.fillBranch(self.outputName+"_mht_l1_ptRatio", vectorSum.Pt()/lepton1.pt)


        return True
