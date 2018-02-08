# Low-energy excess analyzer
This module is a [LArSoft](http://www.larsoft.org) analyzer that builds a ROOT `TTree` with information from &nu;<sub>e</sub> candidates, reconstructed with the [Pandora framework](https://github.com/PandoraPFA).
It looks for neutrino PFParticles with at least one daughter shower and at least one daughter track.

There are two FCL files, one for data (`run_PandoraOnly_data.fcl`) and one for Monte Carlo (`run_PandoraOnly.fcl`).

## Requirements

- `uboonecode v06_26_01_10` (MCC8.6)
- [XGBoost](http://xgboost.readthedocs.io/en/latest/). Follow the instructions in [DocDB 10685](https://microboone-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=10685).

### Instructions to set up xgboost

- install XGboost in the $MRB_INSTALL/xgboost directory (usually it corresponds to the localProducts_*/xgboost directory) following the instructions here http://xgboost.readthedocs.io/en/latest/build.html, meaning:
```
cd $MRB_INSTALL
git clone --recursive https://github.com/dmlc/xgboost
cd xgboost
make -j4
```
- copy the compiled library libxgboost.so to the localProducts folder
```
cp localProducts_*/xgboost/libxgboost.so localProducts_*/uboonecode/v06_26_01_10/slf6.x86_64.e10.prof/lib/
```
- copy the BDT model:
```
cp /uboone/app/users/srsoleti/v06_26_00_2/localProducts_larsoft_v06_30_00_e10_prof/multiclass_pandora_mcc8.model localProducts_*/
```

Now you can compile using at usual:
- got to the `srcs` folder
- `mrb z` to remove all the previous compiled libraries
- `mrbsetenv` to reset the compiler
- `mrb uc` to add the libraries to the cmake file
- `mrb i -j 32` to compile
