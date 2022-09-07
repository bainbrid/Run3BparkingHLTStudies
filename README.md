# Run3BparkingHLTStudies

The instructions below have been verified with python3 (CMSSW_12_0_X onwards).

# 1. log-in to lxplus 
```
cmsrel CMSSW_12_1_0
cd CMSSW_12_1_0/src
cmsenv 
```

# 2. Setup package

```
git clone git@github.com:gitytakahas/Run3BparkingHLTStudies.git
cd Run3BparkingHLTStudies
```

# 3. (Can be skipped) derive efficiecny for single-mu after folding in analysis efficiency. 

```
cd $CMSSW_BASE/src/Run3BparkingHLTStudies/single-mu/
python3 draw_roc.py --pr --weight 
```

This will print out relevant numbers to be later used to make the ROC curve. 
If you just need L1 times HLT trigger eff. just remove ```--weight``` option.
These numbers have been already embedded into the ```draw_roc.py``` script, which is the one used for drawing the ROC curve.


# 4. (Can be skipped) create efficiency map for di-e after folding in analysis efficiency.

```
cd $CMSSW_BASE/src/Run3BparkingHLTStudies/
python3 compare_eff.py --weight
```

This will create the efficiency map to be later used to make the ROC curve. 
If you don't want to fold in analysis efficiency, you can just remove ```--weight``` option. 
Obtained maps are already stored in eos so this step can be skipped unless you want to create the new one.
If you want to make comparison plots for several basic distributions at the HLT level, just add ```--plot``` option.


# 5. (Can be skipped) create rate map for di-electron trigger

```
cd $CMSSW_BASE/src/Run3BparkingHLTStudies/
python3 compare_rate.py
```

This will create the HLT rate map (as a function of number of pileup) to be later used to make the ROC curve. 
Obtained maps are already stored in eos so this step can be skipped unless you want to create the new one.


# 6. Create ROC curve 

```
cd $CMSSW_BASE/src/Run3BparkingHLTStudies/
python3 draw_roc.py --weight
```

The script will loop over relevant npu: 
   * pu = 17 (0.6E34)
   * pu = 25 (0.9E34)
   * pu = 30 (1.1E34)
   * pu = 36 (1.3E34)
   * pu = 42 (1.5E34)
   * pu = 48 (1.7E34)
   * pu = 56 (2.0E34)

and create the ROOT files containing each ROC curve, to be later used in the next step (estimate.py). If you want to have different npu, you can just edit the last few lines of ```draw_roc.py```.


# 5. Estimate number of Kee candidates

Provides the menus and various tables of useful metrics:
```
python3 metrics.py 
```

Provides estimates per LHC fill (based on a pseudo LHC luminosity profiles):
```
python3 estimate2.py 
```
