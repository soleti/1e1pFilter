#ifndef PANDORAINTERFACEHELPER_CXX
#define PANDORAINTERFACEHELPER_CXX

#include "PandoraInterfaceHelper.h"

namespace lee {
  PandoraInterfaceHelper::PandoraInterfaceHelper(){

    _configured = false;

  }
void PandoraInterfaceHelper::get_daughter_tracks(
    std::vector<size_t> pf_ids, const art::Event &evt,
    std::vector<art::Ptr<recob::Track>> &tracks, std::string _pfp_producer) {
      try {

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);

  art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt,
                                               _pfp_producer);

  for (auto const &pf_id : pf_ids) {
      auto const &track_obj = track_per_pfpart.at(pf_id);
      tracks.push_back(track_obj);

  }
} catch (...) {
  std::cout << "[PandoraLEE] "
            << "Error getting daughter tracks" << std::endl;
}
}






void PandoraInterfaceHelper::get_daughter_showers(
    std::vector<size_t> pf_ids, const art::Event &evt,
    std::vector<art::Ptr<recob::Shower>> &showers, std::string _pfp_producer) {

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);

  art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt,
                                                 _pfp_producer);

  for (auto const &pf_id : pf_ids) {
    try {

      auto const &shower_obj = shower_per_pfpart.at(pf_id);
      showers.push_back(shower_obj);
    } catch (...) {
      std::cout << "[PandoraLEE] "
                << "Error getting the shower" << std::endl;
    }
  }
}



  void PandoraInterfaceHelper::Configure(art::Event const & e,
                             std::string _pfp_producer,
                             std::string _spacepoint_producer,
                             std::string _hitfinder_producer,
                             std::string _geant_producer) {

    // Collect hits
    lar_pandora::HitVector hitVector;
    lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinder_producer, hitVector);

    // Collect PFParticles and match Reco Particles to Hits
    lar_pandora::PFParticleVector  recoParticleVector;
    lar_pandora::PFParticleVector  recoNeutrinoVector;
    lar_pandora::PFParticlesToHits pfp_to_hits_map;
    lar_pandora::HitsToPFParticles recoHitsToParticles;

    lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e,
                                                          _pfp_producer,
                                                          _spacepoint_producer,
                                                          pfp_to_hits_map,
                                                          recoHitsToParticles,
                                                          lar_pandora::LArPandoraHelper::kUseDaughters, // Consider daughters as independent pfps
                                                          true); // Use clusters to go from pfp to hits

    if (_verbose) {
      std::cout << "[PandoraInterfaceHelper] RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
      std::cout << "[PandoraInterfaceHelper] RecoParticles: " << recoParticleVector.size() << std::endl;
    }

    // Collect MCParticles and match True Particles to Hits
    lar_pandora::MCParticleVector     trueParticleVector;
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;
    lar_pandora::MCParticlesToHits    trueParticlesToHits;
    lar_pandora::HitsToMCParticles    hit_to_mcps_map;

    if (!e.isRealData()) {
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, trueParticleVector);
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, truthToParticles, particlesToTruth);
      lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e,
                                                            _geant_producer,
                                                            hitVector,
                                                            trueParticlesToHits,
                                                            hit_to_mcps_map,
                                                            lar_pandora::LArPandoraHelper::kUseDaughters); // Consider daughters as independent mcps
    }

    if (_verbose) {
      std::cout << "[PandoraInterfaceHelper] TrueParticles: " << particlesToTruth.size() << std::endl;
      std::cout << "[PandoraInterfaceHelper] TrueEvents: " << truthToParticles.size() << std::endl;
    }

    // Now set the things we need for the future
    _hit_to_mcps_map = hit_to_mcps_map;
    _pfp_to_hits_map = pfp_to_hits_map;

    if (_debug) {
      std::cout << "[PandoraInterfaceHelper] This is event " << e.id().event() << std::endl;
      art::ServiceHandle<cheat::BackTracker> bt;
      std::cout << "[PandoraInterfaceHelper] Number of MCParticles matched to hits: " << trueParticlesToHits.size() << std::endl;
      for (const auto & iter : trueParticlesToHits) {
        const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth((iter.first)->TrackId());
        std::cout << "[PandoraInterfaceHelper] MCParticle with pdg " << (iter.first)->PdgCode()
                  << " and origin " << (mc_truth->Origin() == 1 ? "neutrino" : "cosmic")
                  << " has " << (iter.second).size() << " hits ass." << std::endl;
        if (mc_truth->Origin() == 1) {
          lar_pandora::HitVector hits = (iter.second);
          /*
          for (const auto & hit : hits){
            std::cout << "[PandoraInterfaceHelper]   > Hit on plane " << hit->View()
                      << " on wire " << hit->WireID()
                      << " with time " << hit->PeakTime() << std::endl;
          }
      */
        }
      }
    }

    _configured = true;
  }

  //___________________________________________________________________________________________________
  void PandoraInterfaceHelper::GetRecoToTrueMatches(lar_pandora::PFParticlesToMCParticles &matchedParticles)
  {
    bool _debug = true;

    if (!_configured) {
      std::cout << "Call to " << __PRETTY_FUNCTION__ << " whitout having done configuration. Abort." << std::endl;
      throw std::exception();
    }

    // Loop over the reco particles
    for (auto iter1 : _pfp_to_hits_map) {

      // The PFParticle
      const art::Ptr<recob::PFParticle> recoParticle = iter1.first;

      if (_debug) std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] Looking at PFP with ID " << recoParticle->Self() << std::endl;

      // The PFParticle's hits
      const lar_pandora::HitVector &hitVector = iter1.second;

      if (_debug) std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] \t This PFP has " << hitVector.size() << " hits." << std::endl;

      lar_pandora::MCParticlesToHits truthContributionMap;

      // Loop over all the hits associated to this reco particle
      for (auto hit : hitVector) {

        // Find the MCParticle that share this same hit (if any)
        auto iter3 = _hit_to_mcps_map.find(hit);
        if (_hit_to_mcps_map.end() == iter3)
          continue;

        // If exists, get the MCParticle
        const art::Ptr<simb::MCParticle> trueParticle = iter3->second;

        if (_debug) std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] \t Found a hit shared with MCParticle with PDG " << trueParticle->PdgCode() << std::endl;

        // This map will contain all the true particles that match some or all of the hits of the reco particle
        truthContributionMap[trueParticle].push_back(hit);
      }

      // Now we want to find the true particle that has more hits in common with this reco particle than the others
      lar_pandora::MCParticlesToHits::const_iterator mIter = truthContributionMap.end();

      for (lar_pandora::MCParticlesToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
           iter4 != iterEnd4; ++iter4)
      {
        if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
        {
            mIter = iter4;
        }
      }

      if (truthContributionMap.end() != mIter)
      {
        const art::Ptr<simb::MCParticle> trueParticle = mIter->first;

        if (_debug) std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] \t >>> Match found with MCParticle with PDG " << trueParticle->PdgCode() << std::endl;

        // Emplace into the output map
        matchedParticles[recoParticle] = trueParticle;
      }

    } // _pfp_to_hits_map loop ends

  }

void PandoraInterfaceHelper::GetRecoToTrueMatches(
    art::Event const &e, std::string _pfp_producer,
    std::string _spacepointLabel, std::string _geantModuleLabel,
    std::string _hitfinderLabel,
    lar_pandora::MCParticlesToPFParticles &matchedParticles,
    lar_pandora::MCParticlesToHits &matchedParticleHits) {

  bool _debug = true;

  // --- Collect hits
  lar_pandora::HitVector hitVector;
  lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinderLabel, hitVector);

  // --- Collect PFParticles and match Reco Particles to Hits
  lar_pandora::PFParticleVector recoParticleVector;
  lar_pandora::PFParticleVector recoNeutrinoVector;
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer,
                                                    recoParticleVector);
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector,
                                                           recoNeutrinoVector);
  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(
      e, _pfp_producer, _spacepointLabel, recoParticlesToHits,
      recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

  if (_debug) {
    std::cout << "[PandoraInterfaceHelper] "
              << "  RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
    std::cout << "[PandoraInterfaceHelper] "
              << "  RecoParticles: " << recoParticleVector.size() << std::endl;
  }

  // --- Collect MCParticles and match True Particles to Hits
  lar_pandora::MCParticleVector trueParticleVector;
  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;
  lar_pandora::MCParticlesToHits trueParticlesToHits;
  lar_pandora::HitsToMCParticles trueHitsToParticles;

  if (!e.isRealData()) {
    lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel,
                                                      trueParticleVector);
    lar_pandora::LArPandoraHelper::CollectMCParticles(
        e, _geantModuleLabel, truthToParticles, particlesToTruth);
    lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(
        e, _geantModuleLabel, hitVector, trueParticlesToHits,
        trueHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);
  }

  if (_debug) {
    std::cout << "[PandoraInterfaceHelper] "
              << "  TrueParticles: " << particlesToTruth.size() << std::endl;
    std::cout << "[PandoraInterfaceHelper] "
              << "  TrueEvents: " << truthToParticles.size() << std::endl;
  }

  GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles,
                       matchedParticles, matchedParticleHits);
}

//___________________________________________________________________________________________________
void PandoraInterfaceHelper::GetRecoToTrueMatches(
    const lar_pandora::PFParticlesToHits &recoParticlesToHits,
    const lar_pandora::HitsToMCParticles &trueHitsToParticles,
    lar_pandora::MCParticlesToPFParticles &matchedParticles,
    lar_pandora::MCParticlesToHits &matchedHits, std::string _pfp_producer) {
  PFParticleSet recoVeto;
  MCParticleSet trueVeto;
  bool _recursiveMatching = true;

  GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles,
                       matchedParticles, matchedHits, recoVeto, trueVeto,
                       _recursiveMatching);
}

//___________________________________________________________________________________________________
void PandoraInterfaceHelper::GetRecoToTrueMatches(
    const lar_pandora::PFParticlesToHits &recoParticlesToHits,
    const lar_pandora::HitsToMCParticles &trueHitsToParticles,
    lar_pandora::MCParticlesToPFParticles &matchedParticles,
    lar_pandora::MCParticlesToHits &matchedHits, PFParticleSet &vetoReco,
    MCParticleSet &vetoTrue, bool _recursiveMatching,
    std::string _pfp_producer) {
  bool foundMatches(false);

  // Loop over the reco particles
  for (lar_pandora::PFParticlesToHits::const_iterator
           iter1 = recoParticlesToHits.begin(),
           iterEnd1 = recoParticlesToHits.end();
       iter1 != iterEnd1; ++iter1) {
    const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
    if (vetoReco.count(recoParticle) > 0)
      continue;

    const lar_pandora::HitVector &hitVector = iter1->second;

    lar_pandora::MCParticlesToHits truthContributionMap;

    // Loop over all the hits associated to this reco particle
    for (lar_pandora::HitVector::const_iterator iter2 = hitVector.begin(),
                                                iterEnd2 = hitVector.end();
         iter2 != iterEnd2; ++iter2) {
      const art::Ptr<recob::Hit> hit = *iter2;

      lar_pandora::HitsToMCParticles::const_iterator iter3 =
          trueHitsToParticles.find(hit);
      if (trueHitsToParticles.end() == iter3)
        continue;

      const art::Ptr<simb::MCParticle> trueParticle = iter3->second;
      if (vetoTrue.count(trueParticle) > 0)
        continue;

      // This map will contain all the true particles that match some or all of
      // the hits of the reco particle
      truthContributionMap[trueParticle].push_back(hit);
    }

    // Now we want to find the true particle that has more hits in common with
    // this reco particle than the others
    lar_pandora::MCParticlesToHits::const_iterator mIter =
        truthContributionMap.end();

    for (lar_pandora::MCParticlesToHits::const_iterator
             iter4 = truthContributionMap.begin(),
             iterEnd4 = truthContributionMap.end();
         iter4 != iterEnd4; ++iter4) {
      if ((truthContributionMap.end() == mIter) ||
          (iter4->second.size() > mIter->second.size())) {
        mIter = iter4;
      }
    }

    if (truthContributionMap.end() != mIter) {
      const art::Ptr<simb::MCParticle> trueParticle = mIter->first;

      lar_pandora::MCParticlesToHits::const_iterator iter5 =
          matchedHits.find(trueParticle);

      if ((matchedHits.end() == iter5) ||
          (mIter->second.size() > iter5->second.size())) {
        matchedParticles[trueParticle] = recoParticle;
        matchedHits[trueParticle] = mIter->second;
        foundMatches = true;
      }
    }
  } // recoParticlesToHits loop ends

  if (!foundMatches)
    return;

  for (lar_pandora::MCParticlesToPFParticles::const_iterator
           pIter = matchedParticles.begin(),
           pIterEnd = matchedParticles.end();
       pIter != pIterEnd; ++pIter) {
    vetoTrue.insert(pIter->first);
    vetoReco.insert(pIter->second);
  }

  if (_recursiveMatching)
    GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles,
                         matchedParticles, matchedHits, vetoReco, vetoTrue,
                         _recursiveMatching);
}

void PandoraInterfaceHelper::traversePFParticleTree(
    const art::ValidHandle<std::vector<recob::PFParticle>> pfparticles,
    size_t top_index, std::vector<size_t> &unordered_daugthers,
    std::string _pfp_producer) {

  // This is a tree-traversal algorithm.  It returns the index of the top
  // particle, plus the index
  // of all daughter particles.

  // This is a recursive algorithm, so it needs a break clause:
  unordered_daugthers.push_back(top_index);

  if (pfparticles->at(top_index).Daughters().size() == 0) {
    return;
  }

  // Else, go through the tree:
  for (size_t i = 0; i < pfparticles->at(top_index).Daughters().size(); i++) {
    traversePFParticleTree(pfparticles,
                           pfparticles->at(top_index).Daughters().at(i),
                           unordered_daugthers);
  }

  return;
}

// Method to calculate the total the center for a parent particle (index of
// neutrino pfp)
TVector3 PandoraInterfaceHelper::calculateChargeCenter(
    size_t top_particle_index,
    const art::ValidHandle<std::vector<recob::PFParticle>> pfparticles,
    const art::Event &evt,
    std::string _pfp_producer) {

  // First, get the indexes of pfparticles that are in the hierarchy of this
  // particle:
  std::vector<size_t> daughters;
  daughters.reserve(50);
  traversePFParticleTree(pfparticles, top_particle_index, daughters);

  // Get the associations from pfparticle to spacepoint
  auto const &spacepoint_handle =
      evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);

  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticles, evt,
                                                       _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt,
                                              _pfp_producer);

  // Variables for the total weight and center of charge
  double totalweight = 0;
  std::vector<double> chargecenter;
  chargecenter.resize(3);

  // Loop over the pfparticles, get their space points, and compute the weighted
  // average:

  for (auto &_i_pfp : daughters) {

    // Get the associated spacepoints:
    std::vector<art::Ptr<recob::SpacePoint>> spcpnts =
        spcpnts_per_pfpart.at(_i_pfp);

    // Loop over the spacepoints and get the associated hits:
    for (auto &_sps : spcpnts) {
      auto xyz = _sps->XYZ();
      std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
      // Add the hits to the weighted average, if they are collection hits:
      for (auto &hit : hits) {
        if (hit->View() == geo::kZ) {
          // Collection hits only
          double weight = hit->Integral();
          // std::cout << "[GeometryHelper] " << "Hit Integral: " <<
          // hit->Integral() << std::endl;
          // std::cout << "[GeometryHelper] " << "Hit PeakAmplitude:
          // " << hit->PeakAmplitude() << std::endl;
          // std::cout << "[GeometryHelper] " << "Hit SummedADC: " <<
          // hit->SummedADC() << std::endl;
          chargecenter[0] += (xyz[0]) * weight;
          chargecenter[1] += (xyz[1]) * weight;
          chargecenter[2] += (xyz[2]) * weight;
          totalweight += weight;
          // break; // Exit the loop over hits
        } // if collection

      } // hits

    } // spacepoints

  } // pfparticles

  // Normalize;
  chargecenter[0] /= totalweight;
  chargecenter[1] /= totalweight;
  chargecenter[2] /= totalweight;

  // Store the data:
  TVector3 _center_of_charge;
  _center_of_charge.SetX(chargecenter[0]);
  _center_of_charge.SetY(chargecenter[1]);
  _center_of_charge.SetZ(chargecenter[2]);

  return _center_of_charge;
}

} // namespace lee

#endif
