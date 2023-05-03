# Instructions to run and compile Main_calib.cc

This code is required to dervive calibrations for PF charged hadrons. It is also used to get summary response and resolution plots.
Note: Please run code in ROOT version > 6.xx

You need to add the input root file path in line [995](https://github.com/bkansal/PFCalibration/blob/Run3with_126XGT/PFChargedHadronAnalyzer/test/calibration_code/Main_calib.cc#L995) of Main_calib.cc.

## To get summary response and resolution plots as a function of E(true)
Response and resolution as a function of E(true) are divided into three regions: barrel ($\eta<1.55$), endcap within tracker ($1.55<\eta<2.5$) and endcap outside tracker ($2.5<\eta<2.75$)for each EH and H hadrons.

1. You need to select only one _region_ from the lines [38-40](https://github.com/bkansal/PFCalibration/blob/Run3with_126XGT/PFChargedHadronAnalyzer/test/calibration_code/Main_calib.cc#L38-L41) in the Main_calib.cc code.

2. Then search "summary" in the code you will get list of commented functions. (you can find from lines 2301-2376)

    i)   There we have in total 12 drawGausFit functions for response wrt true energy. 
    
    ii)  There are 6 drawEtaDependence functions for response wrt abs(Eta).
    
    iii) There are 14 calibration coefficients plots for H barrel, H endcap, EH barrel & EH endcap.
    
   You need to uncomment lines to get the corresponding plot.
 
3. To run the code : 
```
	make
	./PFCalib
```

For example :

If you want to look into the calibration coefficients for H barrel then you need to :
1. mention the _region_. ( char* _region_ = (char*)"barrel")
2. uncomment the H barrel calibration coefficients functions only.
Note: comment out all the other plots.    
3. Now, run and complile the code

If you want to look into the final corrected response wrt true energy for H barrel hadrons then you need to :
1. mention the _region_ . (char* _region_ = (char*)"barrel")
2. uncomment the drawGausFit(corrEtaBarrelHcal,responseCor,resolutionCor) function.
Note: comment out all the other plots.    
3. Now, run and complile the code.

If you want to look into the final corrected response wrt eta for H hadrons then you need to :
1. mention the _region_ . (char* _region_ = (char*)"Full")
2. uncomment the drawEtaDependence(corrEtaDependenceH, responseEtaEtaH) function.
Note: comment out all the other plots.    
3. Now, run and complile the code.
