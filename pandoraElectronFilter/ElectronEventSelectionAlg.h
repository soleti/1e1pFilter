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

#include "uboone/SpaceChargeServices/SpaceChargeServiceMicroBooNE.h"

#include "GeometryHelper.h"
#include "PandoraInterfaceHelper.h"


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




    /**
    * @brief Checks if there is a flash within the 3.2-4.8 ms window and compatible with the center of charge
    *
    * @param ipf Index of the PFParticle
    * @param pfparticles PFParticles handle
    * @param _this_center_of_charge Position of the center of charge
    * @param _selected_flash index of the selected flash
    * @param evt art Event
    * @return True or false
    */
    bool opticalfilter(size_t ipf, const std::vector<recob::PFParticle> & pfparticles, TVector3 _this_center_of_charge, int & _selected_flash, const art::Event & evt);

    
    /**
    * @brief Return the true coordinates corrected by the space-charge effect
    *
    * @param xyz TVector3 of the true position
    * @return TVector3 of the space-charge corrected position
    */
    TVector3 spaceChargeTrueToReco(const TVector3 & xyz);

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

    std::string _pfp_producer = "pandoraNu";

    std::string fOpticalFlashFinderLabel;

    // Helper class for geometry functions:
    GeometryHelper geoHelper;

    // Helper class for dealing with pandora heirarchy:
    PandoraInterfaceHelper pandoraHelper;

  };

} // lee

#endif // ELECTRON_EVENT_SELECTION_ALG_H
