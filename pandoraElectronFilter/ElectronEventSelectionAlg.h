////////////////////////////////////////////////////////////////////////
// Class:       ElectronEventSelectionAlg
// Module Type: filter
// File:        ElectronEventSelectionAlg.h
//
////////////////////////////////////////////////////////////////////////

#ifndef ELECTRON_EVENT_SELECTION_ALG_H
#define ELECTRON_EVENT_SELECTION_ALG_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// #include "art/Framework/Principal/Run.h"
// #include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
// #include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"

//#include "uboone/SpaceChargeServices/SpaceChargeServiceMicroBooNE.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightCharge.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"

#include "GeometryHelper.h"
#include "PandoraInterfaceHelper.h"

#include <numeric>


namespace lee
{

class ElectronEventSelectionAlg
{
public:
  // ElectronEventSelectionAlg(){}
  // ~ElectronEventSelectionAlg(){}

  /**
    * @brief Main Event Selection Function
    * @details Decides whether or not this event is an electron neutrino candidate
    *
    * @param evt art::Event containing the information for this event.
    * @return True or False.  True == event passed cuts, false == event failed cuts
    */
  bool eventSelected(const art::Event &evt);

  /**
    * @brief Configure all of the parameters of this class
    *
    * @param p fcl parameter set
    */
  void reconfigure(fhicl::ParameterSet const &p);

  void GetFlashLocation(std::vector<double> pePerOpChannel,
                        double &Ycenter,
                        double &Zcenter,
                        double &Ywidth,
                        double &Zwidth);

  /**
    * @brief Checks if there is a flash within the 3.2-4.8 ms window and compatible with the center of charge
    *
    * @param std::vector<size_t> pfplist : list of primary neutrino candidates that need to be tested
    * @param evt art Event
    * @return -1 if not passed, otherwise index of the flash
    */
  const std::map<size_t, int> opticalfilter(const art::Event &evt,
                                            const std::vector<size_t> &pfplist,
                                            const art::ValidHandle<std::vector<recob::PFParticle>> pfparticle_handle);

  /**
    * @brief Checks if candidates are compatible with the center of charge, the best one is selected using flashmatching
    *
    * @param std::vector<size_t> pfplist : list of primary neutrino candidates that need to be tested
    * @param evt art Event
    * @return negative code if not passed, otherwise index of the flash is returned
    */
  const std::map<size_t, int> flashBasedSelection(const art::Event &evt,
                                                  const std::vector<size_t> &pfplist,
                                                  const art::ValidHandle<std::vector<recob::PFParticle>> pfparticle_handle);

  
  /**
    * @brief Checks if there is a flash within the flash_window_start - flash_window_end window with enough PE, returns the index of -1.
    */
  void flashFinder(const art::Event &evt);

  /**
	* @brief Creates a photon cluster for a neutrino pfp hierarchy
	* PFParticle
	*
	* @param evt art Event
	* @param pfplist list of pfp indices
	* @return flashana::QCluster_t object containing the photons
	*/
  const flashana::QCluster_t collect3DHits(
      const art::Event &evt,
      const std::vector<size_t> &pfplist);

  /**
    * @brief Return the true coordinates corrected by the space-charge effect
    *
    * @param xyz TVector3 of the true position
    * @return TVector3 of the space-charge corrected position
    */
  TVector3 spaceChargeTrueToReco(const TVector3 &xyz);

  /**
    * @brief Reset internal variables
    */
  void clear();

public:
  // Access functions for the saved data:

  /**
    * @brief Returns the number of neutrino candidates from pandora, regardless of whether the passed
    * @return Number of candidates
    */

  /**
    * @brief Return a list of the selected pfparticle top level neutrino candidate indexes
    */
  const std::vector<size_t> &get_primary_indexes() const { return _primary_indexes; }

  /**
    * @brief Return the number of neutrino candidates
    */
  const int &get_n_neutrino_candidates() const { return _n_neutrino_candidates; }

  /**
    * @brief Informs whether a particular candidate passed or failed the algorithm
    * @return Vector of bool, one-to-one with get_primary_indexes
    */
  const std::map<size_t, bool> &get_neutrino_candidate_passed() const { return _neutrino_candidate_passed; }

  /**
    * @brief Informs whether the selected candidate passed or failed the fiducial cut
    * @return bool
    */
  bool get_neutrino_candidate_fiducial() const { return _neutrino_candidate_fiducial; }

  /**
    * @brief Return the index of the flash matched with the pfparticle or failure code, -4 (no intime flash), -3 (PE cut), -2 (precuts), -1 (matching).
    * @details [long description]
    * @return [description]
    */
  const std::map<size_t, int> &get_op_flash_indexes() const { return _op_flash_indexes; }

  /**
    * @brief Return the pandora calculated vertex indexed by pfparticle id number
    * @details [long description]
    * @return [description]
    */
  const std::map<size_t, TVector3> &get_neutrino_vertex() const { return _neutrino_vertex; }

  /**
    * @brief Return number of showers for this pfparticle
    * @details [long description]
    * @return [description]
    */
  const std::map<size_t, int> &get_n_showers() const { return _n_showers; }
  const std::map<size_t, int> &get_n_primary_showers() const { return _n_primary_showers; }

  /**
    * @brief Return number of tracks for pfparticle index
    * @details [long description]
    * @return [description]
    */
  const std::map<size_t, int> &get_n_tracks() const { return _n_tracks; }
  const std::map<size_t, int> &get_n_primary_tracks() const { return _n_primary_tracks; }

  /**
    * @brief Return the list of pfparticle indexes that are showers that are associated with primary pfparticle indexes
    * @details [long description]
    * @return [description]
    */
  const std::map<size_t, std::vector<size_t>> &
  get_pfp_id_showers_from_primary() const { return _pfp_id_showers_from_primary; }

  /**
    * @brief Return the list of pfparticle indexes that are tracks that are associated with primary pfparticle indexes
    * @details [long description]
    * @return [description]
    */
  const std::map<size_t, std::vector<size_t>> &
  get_pfp_id_tracks_from_primary() const { return _pfp_id_tracks_from_primary; }

  /**
    * @brief Return the list of total PE,time of the flashes and the brightest inbeam flash
    * @details [long description]
    * @return [description]
    */
  const std::vector<double> &get_flash_PE() const { return _flash_PE; }
  const std::vector<double> &get_flash_time() const { return _flash_time; }
  const double &get_flash_PE_max() const { return _flash_PE_max; }

  /**
    * @brief Return the positionaL variables of the brightest in time flash
    */
  const double &get_flash_y() const { return _flash_y; }
  const double &get_flash_z() const { return _flash_z; }
  const double &get_flash_sigma_y() const { return _flash_sy; }
  const double &get_flash_sigma_z() const { return _flash_sz; }

  const std::vector<int> &get_flash_matchid() const { return _flash_matchid; }
  const std::vector<double> &get_flash_x() const { return _flash_x; }
  const std::vector<double> &get_TPC_x() const { return _TPC_x; }
  const std::vector<double> &get_flash_score() const { return _flash_score; }
  const std::vector<double> &get_flash_hypo_PE() const { return _flash_hypo_PE; }

  /**
    * @brief Return the  info about the charge distributions of the candidates
    */
  const std::map<size_t, double> &get_candidadates_charge_x() const { return _charges_x; }
  const std::map<size_t, double> &get_candidadates_charge_y() const { return _charges_y; }
  const std::map<size_t, double> &get_candidadates_charge_z() const { return _charges_z; }
  const std::map<size_t, double> &get_candidadates_charge_total() const { return _charges_total; }

protected:
  // Variables that are used to determine the selection and might be worth passing
  // to an analyzer module:

  int _n_neutrino_candidates;
  std::vector<size_t> _primary_indexes;
  std::map<size_t, bool> _neutrino_candidate_passed;
  bool _neutrino_candidate_fiducial;
  std::map<size_t, int> _op_flash_indexes;
  std::map<size_t, TVector3> _neutrino_vertex;

  std::map<size_t, int> _n_showers;
  std::map<size_t, int> _n_primary_showers;
  std::map<size_t, std::vector<size_t>> _pfp_id_showers_from_primary;
  std::map<size_t, int> _n_tracks;
  std::map<size_t, int> _n_primary_tracks;
  std::map<size_t, std::vector<size_t>> _pfp_id_tracks_from_primary;

  std::vector<double> _flash_PE;   // PE of all flashes
  std::vector<double> _flash_time; // time of all flashes
  double _flash_PE_max; // I there was a flash in beamtime, return the one with the highest PE. 

  double _flash_y; // info of the brightest intime flash
  double _flash_sy;
  double _flash_z;
  double _flash_sz;

  std::map<size_t, double> _charges_x;     // closest x spacepoint of the candidates
  std::map<size_t, double> _charges_y;     // center of charge y
  std::map<size_t, double> _charges_z;     // center of charge z
  std::map<size_t, double> _charges_total; // total uncorrected collectionplane charge

  std::vector<int> _flash_matchid;  // The corresponding pfp ID of the matched candidate.
  std::vector<double> _TPC_x;       // candidate lowest x for the matched candidates, in order of score
  std::vector<double> _flash_x;     // flashmatched x for the matched candidates, in order of score
  std::vector<double> _flash_score; // flashmatch score for the matched candidates, in order of score
  std::vector<double> _flash_hypo_PE; // the total PE of the flash hypothesis
protected:
  // Configurable variables from the fcl file:
  int m_nTracks;
  double m_fidvolXstart;
  double m_fidvolXend;

  double m_fidvolYstart;
  double m_fidvolYend;

  double m_fidvolZstart;
  double m_fidvolZend;

  double m_fractionsigmaflashwidth;
  double m_absoluteflashdist;

  double m_startbeamtime;
  double m_endbeamtime;
  double m_PE_threshold;

  // Prematching cuts
  double m_cut_zwidth;
  double m_cut_sigzwidth;
  double m_cut_ywidth;
  double m_cut_sigywidth;
  double m_charge_light_ratio_lower;
  double m_charge_light_ratio_upper;

  // Fixing the PMT wrong ids
  bool _do_opdet_swap;              ///< If true swaps reconstructed OpDets according to _opdet_swap_map
  std::vector<int> _opdet_swap_map; ///< The OpDet swap map for reco flashes

  int m_flash_index;

  bool m_flashmatching; ///< Use flashmatching or the old set of cuts
  bool m_FM_all;        ///< Use flashmatching for all candidates or only the ones passing the topological cut
  bool m_topological;   ///< Use a topological cut of 1 shower + N tracks
  bool m_fiducial;      ///< Use a fiducial cut

  // std::map<unsigned short, double> m_ly_map;

  std::string m_pfp_producer;

  std::string fOpticalFlashFinderLabel;

  // Helper class for geometry functions:
  GeometryHelper geoHelper;

  // Helper class for dealing with pandora heirarchy:
  PandoraInterfaceHelper pandoraHelper;

  flashana::FlashMatchManager m_mgr;
  art::ServiceHandle<geo::Geometry> m_geo;

  double m_ly_proton = 32000;
  double m_ly_electron = 33333;
  double m_ly_muon = 40000;
  double m_ly_gamma = 33333;

  double m_chE_proton = 1.0 / 50;
  double m_chE_electron = 1.0 / 65;
  double m_chE_muon = 1.0 / 110;
  double m_chE_gamma = 1.0 / 65;

  std::map<unsigned short, double> m_ly_map;
};

} // lee

#endif // ELECTRON_EVENT_SELECTION_ALG_H
