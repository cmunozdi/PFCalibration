# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --runUnscheduled --conditions auto:run1_mc -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT,VALIDATION:@standardValidationNoHLT+@miniAODValidation,DQM:@standardDQMFakeHLT+@miniAODDQM --eventcontent RECOSIM,MINIAODSIM,DQM -n 100 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('ana',eras.Run3)

# import of standard configurations
# process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
# process.load('Configuration.EventContent.EventContent_cff')
# process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
# process.load('Configuration.StandardSequences.MagneticField_cff')
# process.load('Configuration.StandardSequences.RawToDigi_cff')
# process.load('Configuration.StandardSequences.L1Reco_cff')
# process.load('Configuration.StandardSequences.Reconstruction_cff')
# process.load('Configuration.StandardSequences.RecoSim_cff')
# process.load('CommonTools.ParticleFlow.EITopPAG_cff')
# process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
# process.load('Configuration.StandardSequences.PATMC_cff')
# process.load('Configuration.StandardSequences.Validation_cff')
# process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 0.2 to 200 GeV withPU samples (x6) 24000 events
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/038ca65d-9cac-440a-b76b-349f3632b238.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/053de016-e900-492a-8912-c6aade92cbcb.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/06a9d0c9-54d6-4a2c-89c3-95a78e3eb442.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/07311d06-f9cc-46a5-9bad-38123d541a80.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0783d628-02b6-47e6-b6b0-7c24bb790276.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/07b54731-53bd-4b01-bea8-427f25202783.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/08412684-d18d-48a4-bda1-19b69c3d56d5.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/08561436-66c2-450d-bb6d-27d97c63b5b2.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/08736dc4-d644-4c06-883d-fc0c740e678e.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/08a7ce0c-c817-4897-8ef3-43a25add9f56.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/08fb30c5-2218-46df-8d8c-d1401304a002.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/09281253-627e-40fb-a4c6-108a07202b7a.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/093aaec7-9dc0-429a-94a6-1a07cab97571.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/09d82275-eaac-4266-bd69-65739da1f857.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0a0538d6-7f3c-416b-bda7-1f6d13906713.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0a0680b1-93ce-4121-89fa-4dc31b2a8750.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0a1b4ce1-0064-4b01-b10f-1ef2c67258ca.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0a1df484-491f-45ea-bc04-ed7964dad973.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0a238f79-81a1-4c05-8401-80845e269861.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0aaa81a1-07e1-4b5a-b0fe-e417a0845e3d.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0aaddc8e-d3e6-484d-a5ac-197873280220.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0aba87bd-49a3-4ef8-950f-f3000026f60c.root",
        "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/142X_mcRun3_2025_realistic_v7-v1/2820000/0ac5dab2-b67a-45bd-b7af-5c9d92d668da.root",

        # # 0.2 to 200 GeV NoPUv2 samples (x6)
        # "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/NoPU_142X_mcRun3_2025_realistic_v7-v1/2820000/03fcaf73-8422-4e6d-90eb-cc9561abe797.root",
        # "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/NoPU_142X_mcRun3_2025_realistic_v7-v1/2820000/04979205-cc24-48ac-a906-63a5b2249475.root",

        # # 0.2 to 200 GeV NoPUv1 samples (x6)
        # "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/NoPU_142X_mcRun3_2025_realistic_v4-v2/2540000/07587c2d-9ba8-4e46-8264-3e488bc9e578.root",
        # "/store/mc/Run3Winter25Reco/Pi_Par-E-0p2to200_PGun/GEN-SIM-RECO/NoPU_142X_mcRun3_2025_realistic_v4-v2/2540000/0789ce43-17cc-4199-bab1-229fd4c51d4f.root",
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    # numberOfThreads = cms.untracked.uint32(16),  # Number of threads
    # numberOfStreams = cms.untracked.uint32(0)   # Number of concurrent event streams
)

# Production Info
# process.configurationMetadata = cms.untracked.PSet(
#     annotation = cms.untracked.string('step3 nevts:100'),
#     name = cms.untracked.string('Applications'),
#     version = cms.untracked.string('$Revision: 1.19 $')
# )

# Output definition

# process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('GEN-SIM-RECO'),
#         filterName = cms.untracked.string('')
#     ),
#     fileName = cms.untracked.string('file:step3.root'),
#     outputCommands = process.RECOSIMEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

# process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
#     compressionAlgorithm = cms.untracked.string('LZMA'),
#     compressionLevel = cms.untracked.int32(4),
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('MINIAODSIM'),
#         filterName = cms.untracked.string('')
#     ),
#     dropMetaData = cms.untracked.string('ALL'),
#     eventAutoFlushCompressedSize = cms.untracked.int32(-900),
#     fastCloning = cms.untracked.bool(False),
#     fileName = cms.untracked.string('file:step3_inMINIAODSIM.root'),
#     outputCommands = process.MINIAODSIMEventContent.outputCommands,
#     overrideBranchesSplitLevel = cms.untracked.VPSet(
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patJets_slimmedJets__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
#             splitLevel = cms.untracked.int32(99)
#         )
#     ),
#     overrideInputFileSplitLevels = cms.untracked.bool(True),
#     splitLevel = cms.untracked.int32(0)
# )

# process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('DQMIO'),
#         filterName = cms.untracked.string('')
#     ),
#     fileName = cms.untracked.string('file:step3_inDQM.root'),
#     outputCommands = process.DQMEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

# Additional output definition

# # Other statements
# process.mix.playback = True
# process.mix.digitizers = cms.PSet()
# for a in process.aliases: delattr(process, a)
# process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2018_realistic_v10', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '105X_mcRun2_asymptotic_v2', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '142X_mcRun3_2025_realistic_v7', '')#'126X_mcRun3_2023_forPU64_v1', '')''133X_mcRun3_2024_realistic_v9
#142X_mcRun3_2025_realistic_v7

process.pfChargedHadronAnalyzer = cms.EDAnalyzer(
    "PFChargedHadronAnalyzer",
    genParticles = cms.InputTag("genParticles"),
    PFCandidates = cms.InputTag("particleFlow"),
    PFSimParticles = cms.InputTag("particleFlowSimParticle"),
    EcalPFClusters = cms.InputTag("particleFlowClusterECAL"),
    HcalPFClusters = cms.InputTag("particleFlowClusterHCAL"),
    ptMin = cms.double(.1),                     # Minimum pt cmunozdi changed this from 1. to 0.1 to include 0p2 to 10 GeV samples                   
    pMin = cms.double(.1),                      # Minimum p  cmunozdi changed this from 1. to 0.1 to include 0p2 to 10 GeV samples                                                                     
    nPixMin = cms.int32(2),                     # Nb of pixel hits                                                                   
    nHitMin = cms.vint32(14,17,20,17,10),       # Nb of track hits                                                                   
    nEtaMin = cms.vdouble(1.4,1.6,2.0,2.4,2.6), # in these eta ranges                                                                
    hcalMin = cms.double(0.),                   # Minimum hcal energy                                                                cmunozdi changed this from 0.5 to 0. to include 0p2 to 10 GeV samples
    ecalMax = cms.double(1E9),                  # Maximum ecal energy                                                                
    verbose = cms.untracked.bool(True),         # not used.                                                                          
    #rootOutputFile = cms.string("PGun__2_200GeV__81X_upgrade2017_realistic_v22.root"),# the root tree                               
    rootOutputFile = cms.string("step3.root"),# the root tree                                                       
#    IsMinBias = cms.untracked.bool(False)                                                                                           
)




process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cff")
#process.load("RecoParticleFlow.Configuration.HepMCCopy_cfi")                                                                        

process.particleFlowSimParticle.ParticleFilter = cms.PSet(
        # Allow *ALL* protons with energy > protonEMin                                                                               
        protonEMin = cms.double(5000.0),
        # Particles must have abs(eta) < etaMax (if close enough to 0,0,0)                                                           
        etaMax = cms.double(5.3),
        # Charged particles with pT < pTMin (GeV/c) are not simulated                                                                
        chargedPtMin = cms.double(0.0),
        # Particles must have energy greater than EMin [GeV]                                                                         
        EMin = cms.double(0.0),
        # half-length of the ECAL endcap inner surface                                                                                                                                                      
        rMax = cms.double(129.),
        zMax = cms.double(317.),
        # List of invisible particles (abs of pdgid)                                                                                                                                                        
        invisibleParticles = cms.vint32()
)
process.genReReco = cms.Sequence(#process.generator+                                                                                 
                                 #process.genParticles+                                                                              
                                 #process.genJetParticles+                                                                           
                                 #process.recoGenJets+                                                                               
                                 #process.genMETParticles+                                                                           
                                 #process.recoGenMET+                                                                                
                                 process.particleFlowSimParticle)



# Path and EndPath definitions


process.EDA = cms.EndPath(process.pfChargedHadronAnalyzer)
process.gRR = cms.EndPath(process.genReReco)



# process.raw2digi_step = cms.Path(process.RawToDigi)
# process.L1Reco_step = cms.Path(process.L1Reco)
# process.reconstruction_step = cms.Path(process.reconstruction)
# process.recosim_step = cms.Path(process.recosim)
# process.eventinterpretaion_step = cms.Path(process.EIsequence)
# process.EDA = cms.EndPath(process.pfChargedHadronAnalyzer)
# process.gRR = cms.EndPath(process.genReReco)




# process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
# process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
# process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
# process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
# process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
# process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
# process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
# process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
# process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
# process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
# process.Flag_METFilters = cms.Path(process.metFilters)
# process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
# process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
# process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
# process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
# process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
# process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
# process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
# process.Flag_ecalBadCalibFilter = cms.Path()
# process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
# process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
# process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
# process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
# process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
# process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
# process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
# process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
# process.prevalidation_step = cms.Path(process.prevalidationNoHLT)
# process.prevalidation_step1 = cms.Path(process.prevalidationMiniAOD)
# process.validation_step = cms.EndPath(process.validationNoHLT)
# process.validation_step1 = cms.EndPath(process.validationMiniAOD)
# process.dqmoffline_step = cms.EndPath(process.DQMOfflineFakeHLT)
# process.dqmoffline_1_step = cms.EndPath(process.DQMOfflineMiniAOD)
# process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
# process.dqmofflineOnPAT_1_step = cms.EndPath(process.PostDQMOfflineMiniAOD)
# process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
# process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)
# process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# # Schedule definition
# #process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.prevalidation_step,process.prevalidation_step1,process.validation_step,process.validation_step1,process.dqmoffline_step,process.dqmoffline_1_step,process.dqmofflineOnPAT_step,process.dqmofflineOnPAT_1_step,process.RECOSIMoutput_step,process.MINIAODSIMoutput_step,process.DQMoutput_step)

# process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.eventinterpretaion_step,process.gRR,process.EDA)


# process.schedule = cms.Schedule(process.gRR,process.EDA)



# process.schedule.associate(process.patTask)



# from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
# associatePatAlgosToolsTask(process)

# # customisation of the process.

# # Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
# from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

# #call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
# process = setCrossingFrameOn(process)

# # End of customisation functions
# #do not add changes to your config after this point (unless you know what you are doing)
# from FWCore.ParameterSet.Utilities import convertToUnscheduled
# process=convertToUnscheduled(process)

# # customisation of the process.

# # Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
# from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

# #call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
# process = miniAOD_customizeAllMC(process)

# # End of customisation functions

# # Customisation from command line

# #Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
# from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
# process = customiseLogErrorHarvesterUsingOutputCommands(process)

# # Add early deletion of temporary data products to reduce peak memory need
# from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
# process = customiseEarlyDelete(process)
# # End adding early deletion
