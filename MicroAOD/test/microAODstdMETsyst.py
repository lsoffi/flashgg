import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAOD")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_asymptotic_v4')

# Fix because auto:run2_mc points to MCRUN2_74_V9::All
current_gt = process.GlobalTag.globaltag.value()
if current_gt.count("::All"):
    new_gt = current_gt.replace("::All","")
    print 'Removing "::All" from GlobalTag by hand for condDBv2: was %s, now %s' % (current_gt,new_gt)
    process.GlobalTag.globaltag = new_gt


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService")
process.RandomNumberGeneratorService.flashggRandomizedPhotons = cms.PSet(
          initialSeed = cms.untracked.uint32(16253245)
        )
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/data/Run2015D/DoubleEG/MINIAOD/05Oct2015-v1/50000/DEDE4FB0-556F\-E511-9F9F-0025905B85E8.root")) 

process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/RunIISpring15MiniAODv2/GluGluHToGG_M-125_13TeV_powheg_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/048F7B0F-0A6E-E511-B710-00259073E37A.root"))
process.MessageLogger.cerr.threshold = 'ERROR' # can't get suppressWarning to work: disable all warnings for now

#===================================Standard Flashgg Sequence + Noise filters===============================================#
process.load("flashgg/MicroAOD/flashggMicroAODSequence_cff")

#================================ Get the most recent JEC ==================================================================#                     
# Setup the private SQLite -- Ripped from PhysicsTools/PatAlgos/test/corMETFromMiniAOD.py                                                         
usePrivateSQlite=True                                                                                                                             
isMC=True                                                                                                                                         
if usePrivateSQlite:                                                                                                                              
    from CondCore.DBCommon.CondDBSetup_cfi import *                                                                                               
    import os
    era = "Summer15_25nsV7"
    if isMC :
        era += "_MC"
    else :
        era += "_DATA"
    dBFile = os.path.expandvars(era+".db")
    process.jec = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        connect = cms.string("sqlite_file:"+dBFile),
        toGet =  cms.VPSet(
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                label= cms.untracked.string("AK4PF")
            ),
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                label= cms.untracked.string("AK4PFchs")
            ),
        )                                                                                                                                         
    )                                                                                                                                             
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')                                                                                  
#===========================================================================================================================#                     
#============================================Apply MET correction and syst.=================================================#
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process, metType="PF",
                           jetCollUnskimmed="slimmedJets",
                           jetColl="slimmedJets",
                           electronColl="slimmedElectrons",
                           photonColl="slimmedPhotons",
                           muonColl="slimmedMuons",
                           tauColl="slimmedTaus",
                           reclusterJets = False,
                           pfCandColl = "packedPFCandidates",
                           postfix="",
                           isData=not isMC,
                           )
#===========================================================================================================================#
#=======================================HCAL Noise Filters==================================================================#
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

#===========================================================================================================================#

#===================================Use the 0 vtx===========================================================================#
from flashgg.MicroAOD.flashggDiPhotons_cfi import flashggDiPhotons
process.flashggDiPhotons0vtx = flashggDiPhotons.clone()
process.flashggDiPhotons0vtx.VertexSelectorName = "FlashggZerothVertexSelector"
#===========================================================================================================================#


from flashgg.MicroAOD.flashggMicroAODOutputCommands_cff import microAODDefaultOutputCommand
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myMicroAODOutputFile.root'),
                               outputCommands = microAODDefaultOutputCommand
                               )
process.p = cms.Path( process.flashggMicroAODSequence+
                      process.HBHENoiseFilterResultProducer+ #produces HBHE bools baseline
                      process.ApplyBaselineHBHENoiseFilter+  #reject events based           
                      process.ApplyBaselineHBHEIsoNoiseFilter
                                #+process.METReader
                      )
process.e = cms.EndPath(process.out)


from flashgg.MicroAOD.MicroAODCustomize import customize
customize(process)

if "DY" in customize.datasetName or "SingleElectron" in customize.datasetName or "DoubleEG" in customize.datasetName:
    customize.customizeHLT(process)
#===========================================================================================================================#
