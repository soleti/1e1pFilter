# Low-energy excess analyzer
This module is a [LArSoft](http://www.larsoft.org) analyzer that builds a ROOT `TTree` with information from &nu;<sub>e</sub> candidates, reconstructed with the [Pandora framework](https://github.com/PandoraPFA).
It looks for neutrino PFParticles with at least one daughter shower and at least one daughter track.

There are two FCL files, one for data (`run_PandoraOnly_data.fcl`) and one for Monte Carlo (`run_PandoraOnly.fcl`).

## Requirements

- `uboonecode v06_26_01_10` (MCC8.6)
- [XGBoost](http://xgboost.readthedocs.io/en/latest/). Follow the instructions in [DocDB 10685](https://microboone-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=10685).

### How to setup uboonecode

- go to your app folder: `/uboone/app/users/USER`
- create Binaries and MyLarsoft: `mkdir Binaries/MyLarsoft/`
- `source /grid/fermiapp/products/uboone/setup_uboone.sh` (must be redone every session)
- look at the available versions of uboonecode: `ups list -aK+ uboonecode`
- suppose you want `v06_26_01_13`
- setup uboonecode: `setup uboonecode v06_26_01_13 -q e10:prof` (must be redone every session)
- create a new dev: `mrb newDev` (in case you need references: https://cdcvs.fnal.gov/redmine/projects/mrb/wiki/MrbRefereceGuide)
- `cd scrs`
- clone the repository of uboonecode, checking out the correct version, which is one of the commit (tag) of the repository `mrb g -t v06_26_01_13 uboonecode`
- `mrbsetenv` (must be redone everytime)
- `mrb i -j 32`

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
cp localProducts_*/xgboost/lib/libxgboost.so localProducts_*/uboonecode/v06_26_01_10/slf6.x86_64.e10.prof/lib/
```
- copy the BDT model:
```
cp /uboone/app/users/wvdp/Binaries/MyLarsoft/localProducts_larsoft_v06_26_01_09_e10_prof/multiclass_pandoraNu_mcc86.model localProducts_*/
```
- add on your .bash_profile the correct exports
```
XGBOOSTSYS="${MRB_INSTALL}/xgboost" 
export XGBOOST_LIB=${XGBOOSTSYS}/lib
export XGBOOST_INC=${XGBOOSTSYS}/include
export RABIT_INC=${XGBOOSTSYS}/rabit/include
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${XGBOOSTSYS}/lib
```
- source again your .bash_profile
```source ~/.bash_profile```

Now you can compile using at usual:
- got to the `srcs` folder
- `mrb z` to remove all the previous compiled libraries
- `mrbsetenv` to reset the compiler
- `mrb uc` to add the libraries to the cmake file
- `mrb i -j 32` to compile
