#!/usr/bin/env python
import os
import time
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.branchselection import BranchSelection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import InputTree
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import eventLoop
from PhysicsTools.NanoAODTools.postprocessing.framework.output import FriendOutput, FullOutput
from PhysicsTools.NanoAODTools.postprocessing.framework.preskimming import preSkim
from PhysicsTools.NanoAODTools.postprocessing.framework.jobreport import JobReport

class PostProcessor :
    def __init__(self,outputDir,inputFiles,cut=None,branchsel=None,modules=[],compression="LZMA:9",friend=False,postfix=None,
		 jsonInput=None,noOut=False,justcount=False,provenance=False,haddFileName=None,fwkJobReport=False,histFileName=None,histDirName=None, outputbranchsel=None,
		 maxEvents=-1,treeName="Events", cutFlow=False):
	self.outputDir=outputDir
	self.inputFiles=inputFiles
	self.cut=cut
	self.modules=modules
	self.compression=compression
	self.postfix=postfix
	self.json=jsonInput
	self.noOut=noOut
	self.friend=friend
	self.justcount=justcount
	self.provenance=provenance
	self.jobReport = JobReport() if fwkJobReport else None
	self.haddFileName=haddFileName
	self.histFile = None
	self.histDirName = None
	self.maxEvents = maxEvents
	self.treeName = treeName
	self.cutFlow = cutFlow
	if self.jobReport and not self.haddFileName :
		print "Because you requested a FJR we assume you want the final hadd. No name specified for the output file, will use tree.root"
		self.haddFileName="tree.root"
 	self.branchsel = BranchSelection(branchsel) if branchsel else None 
        if outputbranchsel != None:
            self.outputbranchsel = BranchSelection(outputbranchsel)
        elif outputbranchsel == None and branchsel != None:
            # Use the same branches in the output as in input
            self.outputbranchsel = BranchSelection(branchsel)
        else: 
            self.outputbranchsel = None

        self.histFileName=histFileName
        self.histDirName=histDirName
    def run(self) :
        outpostfix = self.postfix if self.postfix != None else ("_Friend" if self.friend else "_Skim")
    	if not self.noOut:
            
            if self.compression != "none":
                ROOT.gInterpreter.ProcessLine("#include <Compression.h>")
                (algo, level) = self.compression.split(":")
                compressionLevel = int(level)
                if   algo == "LZMA": compressionAlgo  = ROOT.ROOT.kLZMA
                elif algo == "ZLIB": compressionAlgo  = ROOT.ROOT.kZLIB
                else: raise RuntimeError("Unsupported compression %s" % algo)
            else:
                compressionLevel = 0 
	    print "Will write selected trees to "+self.outputDir
            if not self.justcount:
                if not os.path.exists(self.outputDir):
                    os.system("mkdir -p "+self.outputDir)
        else:
            compressionLevel = 0

	if self.noOut:
	    if len(self.modules) == 0: 
		raise RuntimeError("Running with --noout and no modules does nothing!")

        # Open histogram file, if desired 
        if (self.histFileName != None and self.histDirName == None) or (self.histFileName == None and self.histDirName != None) :
            raise RuntimeError("Must specify both histogram file and histogram directory!")
        elif self.histFileName != None and self.histDirName != None:
            self.histFile = ROOT.TFile.Open( self.histFileName, "RECREATE" )
        else :
            self.histFile = None

        for m in self.modules:
            if hasattr( m, 'writeHistFile') and m.writeHistFile :
                m.beginJob(histFile=self.histFile,histDirName=self.histDirName)
            else :
                m.beginJob()

	fullClone = (len(self.modules) == 0)
	outFileNames=[]
        t0 = time.time()
	totEntriesRead=0
	for friendList in self.inputFiles:

	    fname = friendList[0]
	    # open input file
	    inFile = ROOT.TFile.Open(fname)

	    #get input tree
	    inTree = inFile.Get(self.treeName)
	    for friend in friendList[1:]:
	        inTree.AddFriend(self.treeName,friend)
	    totEntriesRead+=inTree.GetEntries()
	    # pre-skimming
	    elist,jsonFilter = preSkim(inTree, self.json, self.cut)
	    if self.justcount:
		print 'Would select %d entries from %s'%(elist.GetN() if elist else inTree.GetEntries(), fname)
		continue
	    else:
		print 'Pre-select %d entries out of %s '%(elist.GetN() if elist else inTree.GetEntries(),inTree.GetEntries())
		
	    if fullClone:
		# no need of a reader (no event loop), but set up the elist if available
		if elist: inTree.SetEntryList(elist)
	    else:
		# initialize reader
		inTree = InputTree(inTree, elist) 

	    # prepare output file
            if not self.noOut:
                outFileName = os.path.join(self.outputDir, os.path.basename(fname).replace(".root",outpostfix+".root"))
                outFile = ROOT.TFile.Open(outFileName, "RECREATE", "", compressionLevel)
                outFileNames.append(outFileName)
                if compressionLevel: 
                    outFile.SetCompressionAlgorithm(compressionAlgo)
                # prepare output tree
                if self.friend:
                    outTree = FriendOutput(inFile, inTree, outFile)
                else:
                    outTree = FullOutput(
                        inFile,
                        inTree,
                        outFile,
                        branchSelection=self.branchsel,
                        outputbranchSelection=self.outputbranchsel,
                        fullClone=fullClone,
                        jsonFilter=jsonFilter,
                        provenance=self.provenance)
            else : 
                outFile = None
                outTree = None

	    # process events, if needed
	    if not fullClone:
		(nall, npass, timeLoop) = eventLoop(self.modules, inFile, outFile, inTree, outTree,maxEvents=self.maxEvents,cutFlow=self.cutFlow)
		print 'Processed %d preselected entries from %s (%s entries). Finally selected %d entries' % (nall, fname, inTree.GetEntries(), npass)
	    else:
                nall = inTree.GetEntries()
		print 'Selected %d entries from %s' % (outTree.tree().GetEntries(), fname)

	    # now write the output
            if not self.noOut: 
                outTree.write()
                outFile.Close()
                print "Done %s" % outFileName
	    if self.jobReport:
		self.jobReport.addInputFile(fname,nall)
		
	for m in self.modules: m.endJob()
	
	print "Total time %.1f sec. to process %i events. Rate = %.1f Hz." %((time.time()-t0), totEntriesRead, totEntriesRead/(time.time()-t0))


	if self.haddFileName :
		os.system("./haddnano.py %s %s" %(self.haddFileName," ".join(outFileNames))) #FIXME: remove "./" once haddnano.py is distributed with cms releases
	if self.jobReport :
		self.jobReport.addOutputFile(self.haddFileName)
		self.jobReport.save()
