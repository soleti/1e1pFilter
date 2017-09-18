////////////////////////////////////////////////////////////////////////
// Class:       PandoraLEEAnalyzer
// Module Type: analyzer
// File:        PandoraLEEAnalyzer.h
//
////////////////////////////////////////////////////////////////////////

#include <math.h>

#include <fstream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "larsim/MCCheater/BackTracker.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// uncomment the lines below as you use these objects

#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TTree.h"
#include "TVector3.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "uboone/EventWeight/MCEventWeight.h"

#include "TEfficiency.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "ElectronEventSelectionAlg.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "EnergyHelper.h"
#include "GeometryHelper.h"

#include "SpaceChargeMicroBooNE.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace lee {

class PandoraLEEAnalyzer : public art::EDAnalyzer {
public:
  explicit PandoraLEEAnalyzer(fhicl::ParameterSet const &pset);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.
  virtual ~PandoraLEEAnalyzer();

  // Plugins should not be copied or assigned.
  PandoraLEEAnalyzer(PandoraLEEAnalyzer const &) = delete;
  PandoraLEEAnalyzer(PandoraLEEAnalyzer &&) = delete;
  PandoraLEEAnalyzer &operator=(PandoraLEEAnalyzer const &) = delete;
  PandoraLEEAnalyzer &operator=(PandoraLEEAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;
  void endSubRun(const art::SubRun &sr);
  void reconfigure(fhicl::ParameterSet const &pset) override;

private:
  std::string _hitfinderLabel = "pandoraCosmicHitRemoval";
  std::string _geantModuleLabel = "largeant";
  std::string _pfp_producer = "pandoraNu";
  std::string _spacepointLabel = "pandoraNu";
  std::string _mctruthLabel = "generator";

  std::vector<double> _predict_p;
  std::vector<double> _predict_mu;
  std::vector<double> _predict_pi;
  std::vector<double> _predict_em;
  std::vector<double> _predict_cos;

  const anab::CosmicTagID_t TAGID_P = anab::CosmicTagID_t::kGeometry_YY;
  const anab::CosmicTagID_t TAGID_MU = anab::CosmicTagID_t::kGeometry_YZ;
  const anab::CosmicTagID_t TAGID_PI = anab::CosmicTagID_t::kGeometry_ZZ;
  const anab::CosmicTagID_t TAGID_EM = anab::CosmicTagID_t::kGeometry_XX;
  const anab::CosmicTagID_t TAGID_CS = anab::CosmicTagID_t::kGeometry_XY;

  lee::ElectronEventSelectionAlg fElectronEventSelectionAlg;

  EnergyHelper energyHelper;
  GeometryHelper geoHelper;
  PandoraInterfaceHelper pandoraHelper;

  TFile *myTFile;
  TTree *myTTree;
  TTree *myPOTTTree;

  int _interaction_type;

  bool m_isData;
  bool m_printDebug;

  const int k_cosmic = 1;
  const int k_nu_e = 2;
  const int k_nu_mu = 3;
  const int k_nc = 4;
  const int k_dirt = 5;
  const int k_data = 6;
  const int k_other = 0;
  const int k_mixed = 7;
  double _energy;
  int _true_nu_is_fiducial;
  double _nu_energy;

  int _n_tracks;
  int _n_showers;
  double _gain;
  double _vx;
  double _vy;
  double _vz;

  double _true_vx;
  double _true_vy;
  double _true_vz;

  double _true_vx_sce;
  double _true_vy_sce;
  double _true_vz_sce;

  int _nu_matched_tracks;
  int _nu_matched_showers;

  int _nu_pdg;

  int _category;
  int _run;
  int _subrun;
  int _event;
  int _n_candidates;
  int _n_true_nu;
  int _run_sr;
  int _subrun_sr;
  int _n_matched;
  double _pot;
  int _event_passed;
  double _distance;

  std::vector<int> _flash_passed;
  std::vector<int> _track_passed;
  std::vector<int> _shower_passed;
  std::vector<int> _primary_indexes;

  int _n_primaries;
  int _chosen_candidate;

  double _bnbweight;

  std::vector<std::vector<double>> _shower_dQdx;
  std::vector<std::vector<double>> _shower_dEdx;

  std::vector<double> _shower_open_angle;
  std::vector<double> _shower_dir_x;
  std::vector<double> _shower_dir_y;
  std::vector<double> _shower_dir_z;

  std::vector<double> _shower_start_x;
  std::vector<double> _shower_start_y;
  std::vector<double> _shower_start_z;

  std::vector<double> _shower_theta;
  std::vector<double> _shower_phi;

  std::vector<double> _shower_energy;

  std::vector<double> _track_dir_x;
  std::vector<double> _track_dir_y;
  std::vector<double> _track_dir_z;
  std::vector<int> _track_is_fiducial;
  std::vector<int> _shower_is_fiducial;

  std::vector<double> _track_start_x;
  std::vector<double> _track_start_y;
  std::vector<double> _track_start_z;

  std::vector<double> _track_end_x;
  std::vector<double> _track_end_y;
  std::vector<double> _track_end_z;

  std::vector<double> _track_theta;
  std::vector<double> _track_phi;

  std::vector<double> _track_length;
  std::vector<double> _track_id;

  std::vector<double> _track_energy;

  std::vector<int> _nu_daughters_pdg;
  std::vector<double> _nu_daughters_E;

  std::vector<double> _nu_daughters_px;
  std::vector<double> _nu_daughters_py;
  std::vector<double> _nu_daughters_pz;

  std::vector<double> _nu_daughters_vx;
  std::vector<double> _nu_daughters_vy;
  std::vector<double> _nu_daughters_vz;

  std::vector<double> _nu_daughters_endx;
  std::vector<double> _nu_daughters_endy;
  std::vector<double> _nu_daughters_endz;

  std::vector<int> _matched_tracks;
  std::vector<int> _matched_showers;

  double m_dQdxRectangleWidth;
  double m_dQdxRectangleLength;

  size_t choose_candidate(std::vector<size_t> &candidates,
                          const art::Event &evt);
  void clear();
  art::Ptr<recob::Track>
  get_longest_track(std::vector<art::Ptr<recob::Track>> &tracks);
  art::Ptr<recob::Shower>
  get_most_energetic_shower(std::vector<art::Ptr<recob::Shower>> &showers);
  /**
  * @brief Determines if a PFParticle is matched with a MCParticle coming from
  * a neutrino interaction or a cosmic ray
  *
  * @param evt current art Event
  * @param neutrino_pdg array of PDG codes for neutrino-matched PFParticles
  * @param neutrino_pf array of neutrino-matched PFParticles
  * @param cosmic_pdg array of PDG codes for cosmic-matched PFParticles
  * @param cosmic_pf array of cosmic-matched PFParticles
  */
  void categorizePFParticles(
    art::Event const &evt,
    std::vector<int> &neutrino_pdg,
    std::vector<art::Ptr<recob::PFParticle>> &neutrino_pf,
    std::vector<int> &cosmic_pdg,
    std::vector<art::Ptr<recob::PFParticle>> &cosmic_pf);

};

}
