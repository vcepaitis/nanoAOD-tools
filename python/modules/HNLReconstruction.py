import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from utils import deltaR, deltaPhi


class HNLReconstruction(Module):
    def __init__(
        self,
        lepton1Object = lambda event: event.leadingLeptons[0],
        lepton2Object = None,
        jetCollection=lambda event: Collection(event, "Jet"),
        outputName="nominal",
        globalOptions={"isData":False, "isSignal":False},

    ):
        self.lepton1Object = lepton1Object
        self.lepton2Object = lepton2Object
        self.jetCollection = jetCollection

        self.outputName = outputName
        self.globalOptions = globalOptions

    def beginJob(self):
        pass

    def endJob(self):
        pass
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nhnlJet"+self.outputName, "I")
        self.out.branch("hnlJet_"+self.outputName+"_pt", "F" , lenVar="nhnlJet"+self.outputName)
        self.out.branch("hnlJet_"+self.outputName+"_ptsub", "F" , lenVar="nhnlJet"+self.outputName)
        self.out.branch("hnlJet_"+self.outputName+"_ptraw", "F" , lenVar="nhnlJet"+self.outputName)
        self.out.branch("hnlJet_"+self.outputName+"_ptrawsub", "F" , lenVar="nhnlJet"+self.outputName)
        self.out.branch("hnlJet_"+self.outputName+"_ptorig", "F" , lenVar="nhnlJet"+self.outputName)
        self.out.branch("hnlJet_"+self.outputName+"_ptorigsub", "F" , lenVar="nhnlJet"+self.outputName)
        
        self.out.branch("hnlJet_"+self.outputName+"_eta", "F" , lenVar="nhnlJet"+self.outputName)

        if self.lepton2Object==None:
            self.out.branch(self.outputName+"_m_l1j", "F")
            self.out.branch(self.outputName+"_dPt_l1j", "F")
            self.out.branch(self.outputName+"_dPhi_l1j", "F")
            self.out.branch(self.outputName+"_dR_l1j", "F")
            self.out.branch(self.outputName+"_dEta_l1j", "F")
        else:
            self.out.branch(self.outputName+"_m_llj", "F" , lenVar="nhnlJet"+self.outputName)
            self.out.branch(self.outputName+"_dPt_llj", "F" , lenVar="nhnlJet"+self.outputName)
            self.out.branch(self.outputName+"_dPhi_l1j", "F" , lenVar="nhnlJet"+self.outputName)
            self.out.branch(self.outputName+"_dEta_l1j", "F" , lenVar="nhnlJet"+self.outputName)
            self.out.branch(self.outputName+"_dR_l2j", "F", lenVar="nhnlJet"+self.outputName)

        if self.globalOptions["isSignal"]:
            self.out.branch("hnlJet_"+self.outputName+"_isTrueQ", "I" , lenVar="nhnlJet"+self.outputName)
            self.out.branch("hnlJet_"+self.outputName+"_isTrueQE", "I" , lenVar="nhnlJet"+self.outputName)
            self.out.branch("hnlJet_"+self.outputName+"_isTrueQMU", "I", lenVar="nhnlJet"+self.outputName)
            self.out.branch("hnlJet_"+self.outputName+"_isTrueQTAU", "I" , lenVar="nhnlJet"+self.outputName)
            self.out.branch("hnlJet_"+self.outputName+"_trueLxy", "F", lenVar="nhnlJet"+self.outputName)
            
            #number of jets in event that are LLP but are not selected
            self.out.branch(self.outputName+"_nTrueMissed", "I" , lenVar="nhnlJet"+self.outputName)


    def fillTruthInfo(self,hnlJet,otherJets):
        if hnlJet!=None:
            llpq = []
            llpqmu = []
            llpqe = []
            llpqtau = []
            llpxy = []
            for jet in hnlJet :
                llpq.append(jet.isLLP_Q)
                llpqmu.append(jet.isLLP_QMU)
                llpqe.append(jet.isLLP_QE)
                llpqtau.append(jet.isLLP_QTAU)
                llpxy.append(jet.displacement_xy)
            self.out.fillBranch("hnlJet_"+self.outputName+"_isTrueQ", llpq)
            self.out.fillBranch("hnlJet_"+self.outputName+"_isTrueQE", llpqe)
            self.out.fillBranch("hnlJet_"+self.outputName+"_isTrueQMU", llpqmu)
            self.out.fillBranch("hnlJet_"+self.outputName+"_isTrueQTAU", llpqtau)
            self.out.fillBranch("hnlJet_"+self.outputName+"_trueLxy", llpxy)

        else:
            self.out.fillBranch("hnlJet_"+self.outputName+"_isTrueQ", -1)
            self.out.fillBranch("hnlJet_"+self.outputName+"_isTrueQE", -1)
            self.out.fillBranch("hnlJet_"+self.outputName+"_isTrueQMU", -1)
            self.out.fillBranch("hnlJet_"+self.outputName+"_isTrueQTAU", -1)
            self.out.fillBranch("hnlJet_"+self.outputName+"_trueLxy", -10)
            
            
        if len(otherJets)==0:
            self.out.fillBranch(self.outputName+"_nTrueMissed", 0)
        else:
            missed = sum(map(lambda x: x.isLLP_Q+x.isLLP_QE+x.isLLP_QMU+x.isLLP_QTAU,otherJets))
            self.out.fillBranch(self.outputName+"_nTrueMissed", missed)
            
            
    def fillHNLJetInfo(self,hnlJet):
        if hnlJet!=None:
        
            pts = []
            ptssub = []
            ptraws = []
            ptorigins = []
            etas = []
            for jet in hnlJet :
              pts.append(jet.pt)
              ptssub.append(jet.p4Subtracted.Pt())
              ptraws.append(jet.ptRaw)
              ptorigins.append(jet.p4Original.Pt())
              etas.append(jet.eta)

            self.out.fillBranch("hnlJet_"+self.outputName+"_pt", pts)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptsub", ptssub)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptraw", ptraws)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptorig", ptorigins)
            self.out.fillBranch("hnlJet_"+self.outputName+"_eta", etas)
            '''
            self.out.fillBranch("hnlJet_"+self.outputName+"_pt", hnlJet.pt)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptsub", hnlJet.p4Subtracted.Pt())
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptraw", hnlJet.ptRaw)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptrawsub", hnlJet.ptRawSubtracted)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptorig", hnlJet.p4Original.Pt())
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptorigsub", hnlJet.p4OriginalSubtracted.Pt())
            self.out.fillBranch("hnlJet_"+self.outputName+"_eta", hnlJet.eta)
        
        else:
            self.out.fillBranch("hnlJet_"+self.outputName+"_pt", 0)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptsub", 0)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptraw", 0)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptrawsub", 0)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptorig", 0)
            self.out.fillBranch("hnlJet_"+self.outputName+"_ptorigsub", 0)
            self.out.fillBranch("hnlJet_"+self.outputName+"_eta", 0)
        '''

    def analyze(self, event):

        lepton1 = self.lepton1Object(event)
        lepton2 = None if self.lepton2Object==None else self.lepton2Object(event)
          
        jets = self.jetCollection(event)

        if lepton2==None:
            if len(jets) > 0:
                #take jet opposite of l1 in transverse plane
                sortedJets = sorted(jets, key=lambda jet: math.fabs(deltaPhi(jet, lepton1)), reverse=True)
                hnlJet = sortedJets[0]
                self.fillHNLJetInfo(hnlJet)
                
                self.out.fillBranch(self.outputName+"_m_l1j", (lepton1.p4()+hnlJet.p4Subtracted).M())
                self.out.fillBranch(self.outputName+"_dPt_l1j", (lepton1.p4()-hnlJet.p4Subtracted).Pt())
                self.out.fillBranch(self.outputName+"_dPhi_l1j", math.fabs(deltaPhi(lepton1, hnlJet)))
                self.out.fillBranch(self.outputName+"_dR_l1j", deltaR(lepton1, hnlJet))
                self.out.fillBranch(self.outputName+"_dEta_l1j", math.fabs(lepton1.eta-hnlJet.eta))
                setattr(event, "hnlJets_"+self.outputName, [hnlJet])
                
                if self.globalOptions["isSignal"]:
                    self.fillTruthInfo(hnlJet,sortedJets[1:])
                
            else:
                self.fillHNLJetInfo(None)
            
                self.out.fillBranch(self.outputName+"_m_l1j", 0)
                self.out.fillBranch(self.outputName+"_dPt_l1j", 0)
                self.out.fillBranch(self.outputName+"_dPhi_l1j", 0)
                self.out.fillBranch(self.outputName+"_dR_l1j", 0)
                self.out.fillBranch(self.outputName+"_dEta_l1j", 0)
                setattr(event, "hnlJets_"+self.outputName, [])
                
                if self.globalOptions["isSignal"]:
                    self.fillTruthInfo(None,[])
        else:
            if len(jets) > 0:
                #take jet closest to l2
                sortedJets = sorted(jets, key=lambda jet: deltaR(jet, lepton2), reverse=False)
                hnlJet = sortedJets[:2]
                self.fillHNLJetInfo(hnlJet)
                
                dR_l2jet = []
                dEta_l1j = []
                dPhi_l1j = []
                dPt_llj = []
                m_llj = []
                for jet in hnlJet :
                   dR_l2jet.append(deltaR(lepton2,jet))
                   dEta_l1j.append(math.fabs(lepton1.eta-jet.eta))
                   dPhi_l1j.append(math.fabs(deltaPhi(lepton1,jet)))
                   dPt_llj.append((lepton1.p4()-lepton2.p4()-jet.p4Subtracted).Pt())
                   m_llj.append((lepton1.p4()+lepton2.p4()+jet.p4Subtracted).M())
                   
                self.out.fillBranch(self.outputName+"_m_llj", m_llj)
                self.out.fillBranch(self.outputName+"_dPt_llj",dPt_llj )
                self.out.fillBranch(self.outputName+"_dPhi_l1j", dPhi_l1j)
                self.out.fillBranch(self.outputName+"_dEta_l1j", dEta_l1j)
                self.out.fillBranch(self.outputName+"_dR_l2j", dR_l2jet)
                setattr(event, "hnlJets_"+self.outputName, hnlJet)
                

                '''
                self.out.fillBranch(self.outputName+"_m_llj", (lepton1.p4()+lepton2.p4()+hnlJet.p4Subtracted).M())
                self.out.fillBranch(self.outputName+"_dPt_llj", (lepton1.p4()-lepton2.p4()-hnlJet.p4Subtracted).Pt())
                self.out.fillBranch(self.outputName+"_dPhi_l1j", math.fabs(deltaPhi(lepton1,hnlJet)))
                self.out.fillBranch(self.outputName+"_dEta_l1j", math.fabs(lepton1.eta-hnlJet.eta))
                self.out.fillBranch(self.outputName+"_dR_l2j", deltaR(lepton2,hnlJet))
                setattr(event, "hnlJets_"+self.outputName, [hnlJet])
                
                if self.globalOptions["isSignal"]:
                    self.fillTruthInfo(hnlJet,sortedJets[1:])
            
            else:
                self.fillHNLJetInfo(None)
            
                self.out.fillBranch(self.outputName+"_m_llj", 0)
                self.out.fillBranch(self.outputName+"_dPt_llj", 0)
                self.out.fillBranch(self.outputName+"_dPhi_l1j", 0)
                self.out.fillBranch(self.outputName+"_dEta_l1j", 0)
                self.out.fillBranch(self.outputName+"_dR_l2j", 0)
                setattr(event, "hnlJets_"+self.outputName, [])
                
                if self.globalOptions["isSignal"]:
                    self.fillTruthInfo(None,[])
            '''
        return True
