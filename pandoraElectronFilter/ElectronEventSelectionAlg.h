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

#include "TVector3.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
// #include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "larcore/Geometry/Geometry.h"
#include "uboone/SpaceChargeServices/SpaceChargeServiceMicroBooNE.h"


namespace lee {

  typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
  typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
  typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
  typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;

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
  bool eventSelected(const art::Event & evt);


  /**
   * @brief Configure all of the parameters of this class
   *
   * @param p fcl parameter set
   */
  void reconfigure(fhicl::ParameterSet const & p) ;


  // General Worker functions

  /**
   * @brief Determine if the specified point is in the fiducial volume
   *
   * @param x vector of length 3
   * @return True or false
   */
  bool is_fiducial(const std::vector<double> & x) const;

  /**
   * @brief Determine if the specified point is in the fiducial volume
   *
   * @param x TVector3 of 3D location
   * @return True or false
   */
  bool is_fiducial(const TVector3 & x) const;


  /**
  * @brief Determine if the specified point is in the fiducial volume
  *        Not recommended, no array size checking is done.
  *
  * @param x array of 3D location
  * @return True or false
  */
  bool is_fiducial(double x[3]) const;

  /**
   * @brief Compute the 3D distance between two points
   *
   * @param a First Point
   * @param b Second Point
   * @return Returns SQRT( (a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2 )
   */
  double distance(const std::vector<double> & a , const std::vector<double> & b) const ;

  /**
   * @brief Compute the 3D distance between two points
   *
   * @param a First Point
   * @param b Second Point
   * @return Returns SQRT( (a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2 )
   */
  double distance(const TVector3 & a , const TVector3 & b) const ;


  TVector3 calculateChargeCenter(size_t ipf,
                                 const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
                                 const art::Event & evt);

  bool opticalfilter(  size_t ipf,
                       const std::vector<recob::PFParticle> & pfparticles,
                       TVector3 _this_center_of_charge,
                       size_t _selected_flash,
                       const art::Event & evt);

  void traversePFParticleTree(const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
                              size_t top_index,
                              std::vector<size_t> & unordered_daugthers );

  static void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles, lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits);
  /**
  *  @brief Perform matching between true and reconstructed particles
  *
  *  @param recoParticlesToHits the mapping from reconstructed particles to hits
  *  @param trueHitsToParticles the mapping from hits to true particles
  *  @param matchedParticles the output matches between reconstructed and true particles
  *  @param matchedHits the output matches between reconstructed particles and hits
  *  @param recoVeto the veto list for reconstructed particles
  *  @param trueVeto the veto list for true particles
  */
  static void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles, lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits, PFParticleSet &recoVeto, MCParticleSet &trueVeto, bool _recursiveMatching);

  void GetRecoToTrueMatches(art::Event const & e, std::string _pfp_producer, std::string _spacepointLabel, std::string _hitfinderLabel, std::string _geantModuleLabel, lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits);


  /**
   * @brief Reset internal variables
   */
  void clear();

  TVector3 spaceChargeTrueToReco(const TVector3 & xyz);


public:

  // Access functions for the saved data:

  /**
   * @brief Returns the number of neutrino candidates from pandora, regardless of whether the passed
   * @return Number of candidates
   */

  /**
   * @brief Return a list of the selected pfparticle top level neutrino candidate indexes
   */
  const std::vector<size_t> & get_primary_indexes() const {return _primary_indexes;}

  /**
   * @brief Return the number of neutrino candidates
   */
  const size_t & get_n_neutrino_candidates() const {return _n_neutrino_candidates;}

  /**
   * @brief Informs whether a particular candidate passed or failed the algorithm
   * @return Vector of bool, one-to-one with get_primary_indexes
   */
  const std::map<size_t, bool> & get_neutrino_candidate_passed() const {return _neutrino_candidate_passed;}

  /**
   * @brief Return the calculated center of charge, as indexed by the pfparticle id number (it's a map)
   * @details [long description]
   * @return [description]
   */
  const std::map<size_t, TVector3> & get_center_of_charge() const {return _center_of_charge;}

  /**
   * @brief Return the index of the flash matched with the pfparticle
   * @details [long description]
   * @return [description]
   */
  const std::map<size_t, int > & get_op_flash_indexes() const {return _op_flash_indexes;}

  /**
   * @brief Return the pandora calculated vertex indexed by pfparticle id number
   * @details [long description]
   * @return [description]
   */
  const std::map<size_t, TVector3> & get_neutrino_vertex() const {return _neutrino_vertex;}

  /**
   * @brief Return number of showers for this pfparticle
   * @details [long description]
   * @return [description]
   */
  const std::map<size_t, int> & get_n_showers() const {return _n_showers;}

  /**
   * @brief Return number of tracks for pfparticle index
   * @details [long description]
   * @return [description]
   */
  const std::map<size_t, int> & get_n_tracks() const {return _n_tracks;}

  /**
   * @brief Return the list of pfparticle indexes that are showers that are associated with primary pfparticle indexes
   * @details [long description]
   * @return [description]
   */
  const std::map<size_t,  std::vector<size_t> > &
  get_pfp_id_showers_from_primary() const {return _pfp_id_showers_from_primary;}


   /**
    * @brief Return the list of pfparticle indexes that are tracks that are associated with primary pfparticle indexes
    * @details [long description]
    * @return [description]
    */
  const std::map<size_t,  std::vector<size_t> > &
  get_pfp_id_tracks_from_primary() const {return _pfp_id_tracks_from_primary;}


protected:

  // Variables that are used to determine the selection and might be worth passing
  // to an analyzer module:


  size_t _n_neutrino_candidates;
  std::vector<size_t> _primary_indexes;
  std::map<size_t, bool> _neutrino_candidate_passed;
  std::map<size_t, TVector3> _center_of_charge;
  std::map<size_t, int > _op_flash_indexes;
  std::map<size_t, TVector3> _neutrino_vertex;
  std::map<size_t, int> _n_showers;
  std::map<size_t,  std::vector < size_t > > _pfp_id_showers_from_primary;
  std::map<size_t, int> _n_tracks;
  std::map<size_t, std::vector < size_t > > _pfp_id_tracks_from_primary;


protected:

  // Configurable variables from the fcl file:
  int m_nTracks;
  double m_fidvolXstart;
  double m_fidvolXend;

  double m_fidvolYstart;
  double m_fidvolYend;

  double m_fidvolZstart;
  double m_fidvolZend;

  double m_trackLength;

  double m_fractionsigmaflashwidth;
  double m_absoluteflashdist;

  std::string fOpticalFlashFinderLabel;


};

} // lee

#endif // ELECTRON_EVENT_SELECTION_ALG_H
