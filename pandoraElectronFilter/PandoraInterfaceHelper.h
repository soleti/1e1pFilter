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

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace lee {

typedef std::map<art::Ptr<recob::PFParticle>, unsigned int>
    RecoParticleToNMatchedHits;
typedef std::map<art::Ptr<simb::MCParticle>, RecoParticleToNMatchedHits>
    ParticleMatchingMap;
typedef std::set<art::Ptr<recob::PFParticle>> PFParticleSet;
typedef std::set<art::Ptr<simb::MCParticle>> MCParticleSet;

class PandoraInterfaceHelper : public HelperBase {
public:
  PandoraInterfaceHelper() = default;
  ~PandoraInterfaceHelper() = default;


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

private:
};
}

#endif
