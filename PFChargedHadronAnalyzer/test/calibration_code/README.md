# Instructions to run and compile Main_calib.cc

This code is required to dervive calibrations for PF charged hadrons. It is also used to get summary response and resolution plots.
Note: Please run code in ROOT version > 6.xx

You need to add the input root file path in line [995](https://github.com/bkansal/PFCalibration/blob/Run3with_126XGT/PFChargedHadronAnalyzer/test/calibration_code/Main_calib.cc#L995) of Main_calib.cc.

## To get summary response and resolution plots as a function of E(true)
Response as a function of E(true) are categorised into three regions: barrel ($\eta<1.55$), endcap within tracker ($1.55<\eta<2.5$) and endcap outside tracker ($2.5<\eta<2.75$) for each EH and H hadrons. Similarly, for resolution as a function of E(true) will have total six plots.

1. You need to select only one _region_ from the lines [38-40](https://github.com/bkansal/PFCalibration/blob/Run3with_126XGT/PFChargedHadronAnalyzer/test/calibration_code/Main_calib.cc#L38-L41) in the Main_calib.cc code.

2. Then search "summary" in the same code and you will get list of commented functions. (you can find from lines [2522-2555](https://github.com/bkansal/PFCalibration/blob/Run3with_126XGT/PFChargedHadronAnalyzer/test/calibration_code/Main_calib.cc#L2522-L2555))

    i)  There we have in total 12 `drawGausFit` functions which are used to fit 1D raw, energy corrected and eta corrected energy response distributions for different regions in the various E(true) bins. `drawGausFit` function is defined in lines [166-388](https://github.com/bkansal/PFCalibration/blob/Run3with_126XGT/PFChargedHadronAnalyzer/test/calibration_code/Main_calib.cc#L167-L388) and these 1D fitted distributions for various E(true) bins will be saved in projections_*.root file. Final summary plot will be saved in resp_reso_*.root file.
    
    ii) There are 14 calibration coefficients plots for H barrel, H endcap, EH barrel & EH endcap in lines [2578-2601](https://github.com/bkansal/PFCalibration/blob/Run3with_126XGT/PFChargedHadronAnalyzer/test/calibration_code/Main_calib.cc#L2578-L2601).
    
   You need to uncomment lines to get the corresponding plot depending on the region.
 
3. To run the code : 
```
	make
	./PFCalib
```

For example :

If you want to look into the calibration coefficients for H barrel then you need to :
1. mention the _region_. (char* _region_ = (char*)"barrel")
2. uncomment the H barrel calibration coefficients functions only.
Note: comment out all the other plots.    
3. Now, a run and complile the code

If you want to look into the final corrected response wrt true energy for H barrel hadrons then you need to :
1. mention the _region_ . (char* _region_ = (char*)"barrel")
2. uncomment the drawGausFit(corrEtaBarrelHcal,responseCor,resolutionCor) function.
Note: comment out all the other plots.
3. Now, run and complile the code.


## To get summary response plots as a function of $|\eta|$
These response plots as function of $\eta$ covers full eta coverage for both EH and H hadrons.

1. You need to select `Full` region from the line [41](https://github.com/bkansal/PFCalibration/blob/Run3with_126XGT/PFChargedHadronAnalyzer/test/calibration_code/Main_calib.cc#L41) in the Main_calib.cc code.

2. Then search "summary" in the same code and you will get list of commented functions. (you can find from lines [2562-2568](https://github.com/bkansal/PFCalibration/blob/Run3with_126XGT/PFChargedHadronAnalyzer/test/calibration_code/Main_calib.cc#L2562-L2568))
    
    i)  There are 6 `drawEtaDependence` functions which are used to fit 1D raw, energy corrected and eta corrected energy response distributions for EH and H hadrons in the eta fine bin spectrum. These 1D fitted distributions for various E(true) bins will be saved in projectionseta_*.root file and final summary plot will be saved in resp_reso_*_wrtEta.root file,

    ii) You need to uncomment the required lines to get the corresponding plot depending on the correction.

3. To run the code :
```
   	make
	./PFCalib
```

For example :
If you want to look into the final corrected response wrt eta for H hadrons then you need to :
1. mention the _region_ . (char* _region_ = (char*)"Full")
2. uncomment the drawEtaDependence(corrEtaDependenceH, responseEtaEtaH) function.
Note: comment out all the other plots.
3. Now, run and complile the code.
