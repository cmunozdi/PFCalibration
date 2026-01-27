# PFCalibration  
For Single Pion sample Generation 2026:

(For bash)
```
cmsrel CMSSW_15_0_18
cd CMSSW_15_0_18/src
cmsenv
git cms-init
git cms-addpkg RecoParticleFlow/PFSimProducer
```
Now change this line: https://github.com/cms-sw/cmssw/blob/CMSSW_14_2_1/RecoParticleFlow/PFSimProducer/plugins/PFSimParticleProducer.cc#L125 to
```cpp
  edm::ParameterSetDescription psd0;
  psd0.setUnknown();
  desc.add<edm::ParameterSetDescription>("ParticleFilter", psd0);
```
and finally
```
git clone -b PFHC26 https://gitlab.cern.ch/cms-pf/reconstruction/pfhccode.git
```
or
```
git clone -b PFHC26 https://github.com/cmunozdi/PFCalibration.git
```
and continue by compiling this
```
scram b -j 40
cd PFCalibration/PFChargedHadronAnalyzer/test/
cmsenv
```

Using centerally generated Single pion reco sample for Ultralegacy 2016:
```  
cmsRun myEDAna.py 

(for Crab job submission)
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crab_step_analyser.py
```
For PFcalibration use PFCalibration/PFChargedHadronAnalyzer/test/calibration_code/calibChris.C  
Note: Please run calibChris.C code on ROOT version 5.xx (it might crash in ROOT version 6.xx)

