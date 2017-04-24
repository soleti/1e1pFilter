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

#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"

#include "TVector3.h"


#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
// #include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "larcore/Geometry/Geometry.h"




namespace lee {


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
                              std::vector<size_t> unordered_daugthers );

  /**
   * @brief Reset internal variables
   */
  void clear();


public:

  // Access functions for the saved data:

  const size_t & get_n_neutrino_candidates() const {return _n_neutrino_candidates;}
  const std::vector<bool> & get_neutrino_candidate_passed() const {return _neutrino_candidate_passed;}
  const std::vector<TVector3> & get_center_of_charge() const {return _center_of_charge;}
  const std::vector<int > & get_op_flash_indexes() const {return _op_flash_indexes;}
  const std::vector<TVector3> & get_neutrino_vertex() const {return _neutrino_vertex;}
  const std::vector<int> & get_n_showers() const {return _n_showers;}
  const std::vector<int> & get_n_tracks() const {return _n_tracks;}

protected:

  // Variables that are used to determine the selection and might be worth passing
  // to an analyzer module:


  size_t _n_neutrino_candidates;
  std::vector<bool> _neutrino_candidate_passed;
  std::vector<TVector3> _center_of_charge;
  std::vector<int > _op_flash_indexes;
  std::vector<TVector3> _neutrino_vertex;
  std::vector<int> _n_showers;
  std::vector<int> _n_tracks;


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
