////////////////////////////////////////////////////////////////////////
// Class:       PandoraInterfaceHelper
// Module Type: filter
// File:        PandoraInterfaceHelper.h
//
////////////////////////////////////////////////////////////////////////

#ifndef PANDORAINTERFACEHELPER_H
#define PANDORAINTERFACEHELPER_H

#include "HelperBase.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larsim/MCCheater/BackTracker.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace lar_pandora {
  typedef std::map<art::Ptr<recob::PFParticle>, art::Ptr<simb::MCParticle>> PFParticlesToMCParticles;
}

namespace lee {

typedef std::map<art::Ptr<recob::PFParticle>, unsigned int>
    RecoParticleToNMatchedHits;
typedef std::map<art::Ptr<simb::MCParticle>, RecoParticleToNMatchedHits>
    ParticleMatchingMap;
typedef std::set<art::Ptr<recob::PFParticle>> PFParticleSet;
typedef std::set<art::Ptr<simb::MCParticle>> MCParticleSet;

class PandoraInterfaceHelper : public HelperBase {
public:
  PandoraInterfaceHelper();
  ~PandoraInterfaceHelper(){}


  /**
  * @brief Travers the tree of the daughters of a PFParticle
  *
  * @param pfparticles PFParticles handle
  * @param top_index Index of the parent
  * @param unordered_daugthers Vector of PFParticles daughters
  */
  void traversePFParticleTree(
      const art::ValidHandle<std::vector<recob::PFParticle>> pfparticles,
      size_t top_index, std::vector<size_t> &unordered_daugthers,
      std::string _pfp_producer = "pandoraNu");

  /**
   * @brief Measures the three-dimensional center of the deposited charge for a
   * PFParticle
   *
   * @param ipf Index of the PFParticle
   * @param pfparticles PFParticles handle
   * @param evt art Event
   * @return TVector3 of the charge center
   */
  TVector3 calculateChargeCenter(
      size_t ipf,
      const art::ValidHandle<std::vector<recob::PFParticle>> pfparticles,
      const art::Event &evt,
      std::string _pfp_producer = "pandoraNu");

  /**
  *  @brief Perform matching between true and reconstructed particles
  *
  *  @param recoParticlesToHits the mapping from reconstructed particles to hits
  *  @param trueHitsToParticles the mapping from hits to true particles
  *  @param matchedParticles the output matches between reconstructed and true
  * particles
  *  @param matchedHits the output matches between reconstructed particles and
  * hits
  *  @param recoVeto the veto list for reconstructed particles
  *  @param trueVeto the veto list for true particles
  */
  void GetRecoToTrueMatches(
      const lar_pandora::PFParticlesToHits &recoParticlesToHits,
      const lar_pandora::HitsToMCParticles &trueHitsToParticles,
      lar_pandora::MCParticlesToPFParticles &matchedParticles,
      lar_pandora::MCParticlesToHits &matchedHits,
      std::string _pfp_producer = "pandoraNu");

  void GetRecoToTrueMatches(
      const lar_pandora::PFParticlesToHits &recoParticlesToHits,
      const lar_pandora::HitsToMCParticles &trueHitsToParticles,
      lar_pandora::MCParticlesToPFParticles &matchedParticles,
      lar_pandora::MCParticlesToHits &matchedHits, PFParticleSet &recoVeto,
      MCParticleSet &trueVeto, bool _recursiveMatching,
      std::string _pfp_producer = "pandoraNu");

  void
  GetRecoToTrueMatches(art::Event const &e, std::string _pfp_producer,
                       std::string _spacepointLabel,
                       std::string _hitfinderLabel,
                       std::string _geantModuleLabel,
                       lar_pandora::MCParticlesToPFParticles &matchedParticles,
                       lar_pandora::MCParticlesToHits &matchedHits);

  void get_daughter_tracks(std::vector<size_t> pf_ids, const art::Event &evt,
                           std::vector<art::Ptr<recob::Track>> &tracks,
                           std::string _pfp_producer = "pandoraNu");

  void get_daughter_showers(std::vector<size_t> pf_ids, const art::Event &evt,
                            std::vector<art::Ptr<recob::Shower>> &showers,
                            std::string _pfp_producer = "pandoraNu");


    /// Configure function parameters
     /**
     *  @brief Configure function parameters (call this function first)
     *
     *  @param e the art::Event
     *  @param _pfp_producer the PFParticle producer label
     *  @param _spacepoint_producer the SpacePoint producer label
     *  @param _hitfinder_producer the Hit producer label
     *  @param _geant_producer The Geant4 producer label
     */
    void Configure(art::Event const & e, std::string _pfp_producer, std::string _spacepoint_producer, std::string _hitfinder_producer, std::string _geant_producer);

    /**
    *  @brief Returns matching between true and reconstructed particles
    *
    *  @param matchedParticles the output matches between reconstructed and true particles
    */
   void GetRecoToTrueMatches(lar_pandora::PFParticlesToMCParticles & matchedParticles);
 protected:

   lar_pandora::HitsToMCParticles _hit_to_mcps_map; ///< A map from recon hits to MCParticles
   lar_pandora::PFParticlesToHits _pfp_to_hits_map; ///< A map from PFParticles to recon hits
 private:

   bool _configured = false;

   bool _debug      = false;
   bool _verbose    = false;

};
}

#endif
