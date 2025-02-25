from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#No PU
# config.General.requestName = 'Offline_NTuples_2025_0p2to10_v7_1'
# config.General.requestName = 'Offline_NTuples_2025_0p2to200_v7_1'
# config.General.requestName = 'Offline_NTuples_2025_200to500_v7_1'
# config.General.requestName = 'Offline_NTuples_2025_500to5000_v7_1'
#With PU
# config.General.requestName = 'Offline_NTuples_2025_withPU_0p2to200_v7_1'
config.General.requestName = 'Offline_NTuples_2025_withPU_200to500_v7_1'
config.General.workArea = 'crab_projects'

#optional
#config.General.transferOutputs
#config.General.transferLogs
#config.General.failureLimit =

#Expert use
#config.General.instance
#config.General.activity

config.section_("JobType")
config.JobType.psetName = 'myEDAna.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'step3_RAW2DIGI_L1Reco_RECO_EI_PAT_VALIDATION_DQM.py'
config.JobType.outputFiles = ['step3.root']
#config.JobType.eventsPerLumi = 2000
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 2000

config.section_("Data")
config.Data.allowNonValidInputDataset = True
# No PU
# config.Data.inputDataset = '/Pi_Par-E-0p2to10_PGun/Run3Winter25Reco-NoPU_142X_mcRun3_2025_realistic_v7-v2/GEN-SIM-RECO'
# config.Data.inputDataset = '/Pi_Par-E-0p2to200_PGun/Run3Winter25Reco-NoPU_142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RECO'
# config.Data.inputDataset = '/Pi_Par-E-200to500_PGun/Run3Winter25Reco-NoPU_142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RECO'
# config.Data.inputDataset = '/Pi_Par-E-500to5000_PGun/Run3Winter25Reco-NoPU_142X_mcRun3_2025_realistic_v7-v2/GEN-SIM-RECO'
# With PU
# config.Data.inputDataset = '/Pi_Par-E-0p2to200_PGun/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RECO'
config.Data.inputDataset = '/Pi_Par-E-200to500_PGun/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RECO'


#MIKKO's new request for 0.2 to 10 GeV
#config.Data.inputDataset = '/SinglePionGun_E0p2to10/Run3Winter24Reco-NoPU_133X_mcRun3_2024_realistic_v9-v2/GEN-SIM-RECO'
#config.Data.inputDataset = '/SinglePionGun_E500to5000/Run3Winter24Reco-NoPU_133X_mcRun3_2024_realistic_v9-v2/GEN-SIM-RECO'


#config.Data.primaryDataset = ''
#config.Data.splitting = 'EventBased'
# config.Data.userInputFiles = open('/afs/cern.ch/user/b/bkansal/work/PFcalibration_2019/CMSSW_10_6_0_pre2/src/PFCalibration/PFChargedHadronAnalyzer/test/step2_2_200_10million_1060_2018.txt').readlines()
# config.Data.ignoreLocality = True
config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 180
config.Data.totalUnits = -1
#NJOBS = 5000
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
#config.Data.publishDBS = '' default for the moment
#config.Data.outLFN = '/home/spandey/t3store/PF_PGun'
#No PU
# config.Data.outLFNDirBase = '/store/user/cmunozdi/OFFLINE_NTUPLES/2025/NoPU_v7/RawNTuples/2025_E0p2to10GeV_v7_1/'
# config.Data.outLFNDirBase = '/store/user/cmunozdi/OFFLINE_NTUPLES/2025/NoPU_v7/RawNTuples/2025_E0p2to200GeV_v7_1/'
# config.Data.outLFNDirBase = '/store/user/cmunozdi/OFFLINE_NTUPLES/2025/NoPU_v7/RawNTuples/2025_E200to500GeV_v7_1/'
# config.Data.outLFNDirBase = '/store/user/cmunozdi/OFFLINE_NTUPLES/2025/NoPU_v7/RawNTuples/2025_E500to5000GeV_v7_1/'
#With PU
# config.Data.outLFNDirBase = '/store/user/cmunozdi/OFFLINE_NTUPLES/2025/WithPU_v7/RawNTuples_withPU/2025_E0p2to200GeV_v7/'
config.Data.outLFNDirBase = '/store/user/cmunozdi/OFFLINE_NTUPLES/2025/WithPU_v7/RawNTuples_withPU/2025_E200to500GeV_v7/'


config.section_("Site")
config.Site.storageSite = 'T3_CH_CERNBOX'
# config.Site.blacklist = ['T2_CH_*', 'T2_US_*']
# # config.Site.whitelist = ['T2_UK_London_IC','T2_US_*']
# config.Site.whitelist = ['T2_UA_KIPT','T2_IN_TIFR']

# config.Site.whitelist = ['T2_AT_Vienna','T2_UK_London_IC']


#config.section_("User")
#config.section_("Debug")
