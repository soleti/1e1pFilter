# Low-energy excess analyzer
This module is a [LArSoft](http://www.larsoft.org) analyzer that builds a ROOT `TTree` with information from &nu;<sub>e</sub> candidates, reconstructed with the [Pandora framework](https://github.com/PandoraPFA).
It looks for neutrino PFParticles with at least one daughter shower and at least one daughter track, or at least two daughter showers in the full Pandora PFParticle hierarchy.

There are two FCL files, one for data (`run_PandoraOnly_data_bnb.fcl`, `run_PandoraOnly_data_extbnb.fcl`) and one for Monte Carlo (`run_PandoraOnly.fcl`).

## Requirements
- `uboonecode v06_26_01_17` (MCC8.7)
- Checkout feature `feature/alister1_EventWeightTreeUtility`
- Checkout feature `feature/kwoodruf_DecisionTreeID`
- Checkout feature `wvdp_lightcharge_v6_26`
- Copy XGBoost model from `/uboone/data/users/kwoodruf/treefiles/multiclass_pandoraNu_mcc86.model` into the `localProducts` folder
- Copy `/uboone/app/users/srsoleti/v06_26_01_11/localProducts_larsoft_v06_26_01_10_e10_prof/uboonecode/v06_26_01_12/slf6.x86_64.e10.prof/lib/libxgboost.so` into the equivalent directory in your local installation
- [XGBoost](http://xgboost.readthedocs.io/en/latest/). Follow the instructions in [DocDB 10685](https://microboone-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=10685).
