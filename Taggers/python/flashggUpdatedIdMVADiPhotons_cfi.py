import FWCore.ParameterSet.Config as cms

flashggUpdatedIdMVADiPhotons = cms.EDProducer("FlashggDiPhotonWithUpdatedPhoIdMVAProducer",
                                              src                      = cms.InputTag("flashggDiPhotons"),
                                              rhoFixedGridCollection   = cms.InputTag('fixedGridRhoAll'),
                                              #photonIdMVAweightfile_EB = cms.FileInPath("flashgg/MicroAOD/data/MVAweights_76X_25ns_r9shift_barrel.xml"),
                                              photonIdMVAweightfile_EB = cms.FileInPath("flashgg/MicroAOD/data/MVAweights_76X_25ns_barrel.xml"),
                                              photonIdMVAweightfile_EE = cms.FileInPath("flashgg/MicroAOD/data/MVAweights_76X_25ns_endcap.xml"),
                                              correctionFile           = cms.FileInPath("flashgg/MicroAOD/data/transformation_76X_v2.root"),
                                              Debug                    = cms.bool(False)
                                              )
