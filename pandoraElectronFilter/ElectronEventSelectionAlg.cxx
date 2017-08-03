#ifndef ELECTRON_EVENT_SELECTION_ALG_CXX
#define ELECTRON_EVENT_SELECTION_ALG_CXX

#include "ElectronEventSelectionAlg.h"

namespace lee {

  void ElectronEventSelectionAlg::GetRecoToTrueMatches(art::Event const & e,
                                             std::string _pfp_producer,
                                             std::string _spacepointLabel,
                                             std::string _geantModuleLabel,
                                             std::string _hitfinderLabel,
                                             lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                             lar_pandora::MCParticlesToHits &matchedParticleHits)
  {

     bool _debug = true;

    // --- Collect hits
    lar_pandora::HitVector hitVector;
    lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinderLabel, hitVector);

     // --- Collect PFParticles and match Reco Particles to Hits
    lar_pandora::PFParticleVector  recoParticleVector;
    lar_pandora::PFParticleVector  recoNeutrinoVector;
    lar_pandora::PFParticlesToHits recoParticlesToHits;
    lar_pandora::HitsToPFParticles recoHitsToParticles;

    lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

    if (_debug) {
      std::cout << "[ElectronEventSelectionAlg] " << "  RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
      std::cout << "[ElectronEventSelectionAlg] " << "  RecoParticles: " << recoParticleVector.size() << std::endl;
    }

    // --- Collect MCParticles and match True Particles to Hits
    lar_pandora::MCParticleVector     trueParticleVector;
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;
    lar_pandora::MCParticlesToHits    trueParticlesToHits;
    lar_pandora::HitsToMCParticles    trueHitsToParticles;

    if (!e.isRealData()) {
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, trueParticleVector);
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, truthToParticles, particlesToTruth);
      lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e, _geantModuleLabel, hitVector, trueParticlesToHits, trueHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);
    }

    if (_debug) {
      std::cout << "[ElectronEventSelectionAlg] " << "  TrueParticles: " << particlesToTruth.size() << std::endl;
      std::cout << "[ElectronEventSelectionAlg] " << "  TrueEvents: " << truthToParticles.size() << std::endl;
    }


    GetRecoToTrueMatches(recoParticlesToHits,
                         trueHitsToParticles,
                         matchedParticles,
                         matchedParticleHits);
  }


  //___________________________________________________________________________________________________
  void ElectronEventSelectionAlg::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits,
                                             const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                             lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                             lar_pandora::MCParticlesToHits &matchedHits)
  {
      PFParticleSet recoVeto; MCParticleSet trueVeto;
      bool _recursiveMatching = true;

      GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, recoVeto, trueVeto, _recursiveMatching);
  }

  //___________________________________________________________________________________________________
  void ElectronEventSelectionAlg::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits,
                                             const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                             lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                             lar_pandora::MCParticlesToHits &matchedHits,
                                             PFParticleSet &vetoReco,
                                             MCParticleSet &vetoTrue,
                                             bool _recursiveMatching)
  {
      bool foundMatches(false);

      // Loop over the reco particles
      for (lar_pandora::PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
          iter1 != iterEnd1; ++iter1)
      {
          const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
          if (vetoReco.count(recoParticle) > 0)
              continue;

          const lar_pandora::HitVector &hitVector = iter1->second;

          lar_pandora::MCParticlesToHits truthContributionMap;

          // Loop over all the hits associated to this reco particle
          for (lar_pandora::HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
          {
              const art::Ptr<recob::Hit> hit = *iter2;

              lar_pandora::HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
              if (trueHitsToParticles.end() == iter3)
                  continue;

              const art::Ptr<simb::MCParticle> trueParticle = iter3->second;
              if (vetoTrue.count(trueParticle) > 0)
                  continue;

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

              lar_pandora::MCParticlesToHits::const_iterator iter5 = matchedHits.find(trueParticle);

              if ((matchedHits.end() == iter5) || (mIter->second.size() > iter5->second.size()))
              {
                  matchedParticles[trueParticle] = recoParticle;
                  matchedHits[trueParticle] = mIter->second;
                  foundMatches = true;
              }
          }
      } // recoParticlesToHits loop ends

      if (!foundMatches)
          return;

      for (lar_pandora::MCParticlesToPFParticles::const_iterator pIter = matchedParticles.begin(), pIterEnd = matchedParticles.end();
          pIter != pIterEnd; ++pIter)
      {
          vetoTrue.insert(pIter->first);
          vetoReco.insert(pIter->second);
      }

      if (_recursiveMatching)
          GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, vetoReco, vetoTrue, _recursiveMatching);
  }


  //_____________________

void ElectronEventSelectionAlg::clear() {
  _n_neutrino_candidates = 0.0;
  _primary_indexes.clear();
  _neutrino_candidate_passed.clear();
  _center_of_charge.clear();
  _op_flash_indexes.clear();
  _neutrino_vertex.clear();
  _n_showers.clear();
  _n_tracks.clear();
  _pfp_id_showers_from_primary.clear();
  _pfp_id_tracks_from_primary.clear();
}

bool ElectronEventSelectionAlg::is_fiducial(const std::vector<double> & x) const {
  if (x.size() != 3) {
    return false;
  }

  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {0., 2.*geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(), 0., geo->DetLength()};

  bool is_x = x[0] > (bnd[0] + m_fidvolXstart) && x[0] < (bnd[1] - m_fidvolXend);
  bool is_y = x[1] > (bnd[2] + m_fidvolYstart) && x[1] < (bnd[3] - m_fidvolYend);
  bool is_z = x[2] > (bnd[4] + m_fidvolZstart) && x[2] < (bnd[5] - m_fidvolZend);
  return is_x && is_y && is_z;
}

bool ElectronEventSelectionAlg::is_fiducial(const TVector3 & x) const {

  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {0., 2.*geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(), 0., geo->DetLength()};

  bool is_x = x[0] > (bnd[0] + m_fidvolXstart) && x[0] < (bnd[1] - m_fidvolXend);
  bool is_y = x[1] > (bnd[2] + m_fidvolYstart) && x[1] < (bnd[3] - m_fidvolYend);
  bool is_z = x[2] > (bnd[4] + m_fidvolZstart) && x[2] < (bnd[5] - m_fidvolZend);
  return is_x && is_y && is_z;
}

bool ElectronEventSelectionAlg::is_fiducial(double x[3]) const {

  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {0., 2.*geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(), 0., geo->DetLength()};

  bool is_x = x[0] > (bnd[0] + m_fidvolXstart) && x[0] < (bnd[1] - m_fidvolXend);
  bool is_y = x[1] > (bnd[2] + m_fidvolYstart) && x[1] < (bnd[3] - m_fidvolYend);
  bool is_z = x[2] > (bnd[4] + m_fidvolZstart) && x[2] < (bnd[5] - m_fidvolZend);
  return is_x && is_y && is_z;
}


double ElectronEventSelectionAlg::distance(const std::vector<double> & a , const std::vector<double> & b) const {

  if (a.size() != 3 || b.size() != 3) {
    return -1;
  }

  double d = 0;

  for (int i = 0; i < 3; i++)
  {
    d += pow((a[i] - b[i]), 2);
  }

  return sqrt(d);
}

double ElectronEventSelectionAlg::distance(const TVector3 & a , const TVector3 & b) const {

  return (a - b).Mag();
}


TVector3 ElectronEventSelectionAlg::spaceChargeTrueToReco(const TVector3 & xyz) {
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  geo::Point_t point(xyz);
  // auto correction = sce->GetPosOffsets(point);
  TVector3 correctedPoint(xyz);
  correctedPoint.SetX(xyz.X() - sce->GetPosOffsets(point).X() + 0.7);
  correctedPoint.SetX(xyz.Y() + sce->GetPosOffsets(point).Y() );
  correctedPoint.SetX(xyz.Z() + sce->GetPosOffsets(point).Z() );
  return correctedPoint;
}



void ElectronEventSelectionAlg::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  m_nTracks = p.get<int>("nTracks", 1);

  m_trackLength = p.get<int>("trackLength", 100);

  m_fidvolXstart = p.get<double>("fidvolXstart", 10);
  m_fidvolXend = p.get<double>("fidvolXstart", 10);

  m_fidvolYstart = p.get<double>("fidvolYstart", 20);
  m_fidvolYend = p.get<double>("fidvolYend", 20);

  m_fidvolZstart = p.get<double>("fidvolZstart", 10);
  m_fidvolZend = p.get<double>("fidvolZend", 50);

  m_fractionsigmaflashwidth = p.get<double>("fractionsigmaflashwidth", 2.0);
  m_absoluteflashdist = p.get<double>("absoluteflashdist", 50.0);

  fOpticalFlashFinderLabel = p.get<std::string>("OpticalFlashFinderLabel", "simpleFlashBeam");
}


void ElectronEventSelectionAlg::traversePFParticleTree(
  const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
  size_t top_index,
  std::vector<size_t> & unordered_daugthers )
{

  // This is a tree-traversal algorithm.  It returns the index of the top particle, plus the index
  // of all daughter particles.

  // This is a recursive algorithm, so it needs a break clause:
  unordered_daugthers.push_back(top_index);

  if (pfparticles->at(top_index).Daughters().size() == 0) {
    return;
  }

  // Else, go through the tree:
  for (size_t i = 0; i < pfparticles->at(top_index).Daughters().size(); i ++) {
    traversePFParticleTree(pfparticles, pfparticles->at(top_index).Daughters().at(i), unordered_daugthers );
  }

  return;

}


// Method to calculate the total the center for a parent particle (index of neutrino pfp)
TVector3 ElectronEventSelectionAlg::calculateChargeCenter(
  size_t top_particle_index,
  const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
  const art::Event & evt)
{

  // First, get the indexes of pfparticles that are in the hierarchy of this particle:
  std::vector<size_t> daughters;
  daughters.reserve(50);
  traversePFParticleTree( pfparticles, top_particle_index, daughters);


  // Get the associations from pfparticle to spacepoint
  art::InputTag pandoraNu_tag { "pandoraNu" };
  auto const& spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(pandoraNu_tag);

  art::FindManyP<recob::SpacePoint > spcpnts_per_pfpart ( pfparticles, evt, pandoraNu_tag );
  art::FindManyP<recob::Hit > hits_per_spcpnts ( spacepoint_handle, evt, pandoraNu_tag );


  // Variables for the total weight and center of charge
  double totalweight = 0;
  std::vector<double> chargecenter;
  chargecenter.resize(3);

  // Loop over the pfparticles, get their space points, and compute the weighted average:

  for (auto & _i_pfp : daughters) {

    // Get the associated spacepoints:
    std::vector<art::Ptr < recob::SpacePoint > > spcpnts = spcpnts_per_pfpart.at(_i_pfp);

    // Loop over the spacepoints and get the associated hits:
    for (auto & _sps : spcpnts) {
      auto xyz = _sps->XYZ();
      std::vector<art::Ptr<recob::Hit> > hits = hits_per_spcpnts.at(_sps.key());
      // Add the hits to the weighted average, if they are collection hits:
      for (auto & hit : hits) {
        if (hit->View() == geo::kZ) {
          // Collection hits only
          double weight = hit->Integral();
          // std::cout << "[ElectronEventSelectionAlg] " << "Hit Integral: " << hit->Integral() << std::endl;
          // std::cout << "[ElectronEventSelectionAlg] " << "Hit PeakAmplitude: " << hit->PeakAmplitude() << std::endl;
          // std::cout << "[ElectronEventSelectionAlg] " << "Hit SummedADC: " << hit->SummedADC() << std::endl;
          chargecenter[0] += (xyz[0]) * weight;
          chargecenter[1] += (xyz[1]) * weight;
          chargecenter[2] += (xyz[2]) * weight;
          totalweight += weight;
          //break; // Exit the loop over hits
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


bool ElectronEventSelectionAlg::opticalfilter(
  size_t ipf,
  const std::vector<recob::PFParticle> & pfparticles,
  TVector3 _this_center_of_charge,
  int &_selected_flash,
  const art::Event & evt)
{

  bool pass = false;

  art::InputTag optical_tag   {fOpticalFlashFinderLabel};
  auto const& optical_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(optical_tag);

  double par1 = m_fractionsigmaflashwidth;     // fraction of flash sigma that is requested
  double par2 = m_absoluteflashdist;    // or also ok if it is closer than this number in cm

  // Run the calculation of charge center:
  int pass_index = -1;
  for (unsigned int ifl = 0; ifl < optical_handle->size(); ++ifl)
  {
    recob::OpFlash const& flash = optical_handle->at(ifl);
    if ((flash.Time() > 4.8 || flash.Time() < 3.2)) continue;
    bool sigma    = (flash.ZCenter() + flash.ZWidth() / par1) > _this_center_of_charge.Z() &&
                    (flash.ZCenter() - flash.ZWidth() / par1) < _this_center_of_charge.Z();
    bool absolute = std::abs(flash.ZCenter() - _this_center_of_charge.Z()) < par2;
    // std::cout << "[ElectronEventSelectionAlg] " << "The flash time is " << flash.Time()
    //           << ", Zcentre: " << flash.ZCenter()
    //           << " and the Zwidth: " << flash.ZWidth()
    //           << std::endl;
    // std::cout << "[ElectronEventSelectionAlg] " << "Z Center of charge is " << _this_center_of_charge.Z() << std::endl;
    if (sigma || absolute)
    {
      pass = true;
      pass_index = ifl;
    }
  }


  if (pass_index != -1) {
    _selected_flash = pass_index;

  }

  return pass;


}


bool ElectronEventSelectionAlg::eventSelected(const art::Event & evt)
{

  art::InputTag pandoraNu_tag { "pandoraNu" };

  auto const& shower_handle = evt.getValidHandle<std::vector<recob::Shower>>(pandoraNu_tag);

  auto const& showers_handle(*shower_handle);
  for (auto & shower: showers_handle ) {
    for (size_t i = 0; i < shower.Energy().size();i++) {
      std::cout << "[ElectronEventSelectionAlg] " << "E " << shower.Energy()[i] << std::endl;
    }
  }

  clear();

  // Get the list of pfparticles:
  auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( pandoraNu_tag );

  // Are there any pfparticles?
  if (pfparticle_handle -> size() == 0) {
    std::cout << "[ElectronEventSelectionAlg] " << "NO RECO DATA PRODUCTS" << std::endl;
    return false;
  }

  // Get the list of primary pfparticles that are also neutrinos (numu or nue)
  for (size_t _i_pfp = 0; _i_pfp < pfparticle_handle -> size(); _i_pfp ++) {

    if ( ( abs( pfparticle_handle->at(_i_pfp).PdgCode() ) == 12 ||
           abs( pfparticle_handle->at(_i_pfp).PdgCode() ) == 14 )
         && pfparticle_handle->at(_i_pfp).IsPrimary() ) {
      _primary_indexes.push_back(_i_pfp);
    }

  }

  // If there are no particles flagged as primary, return false
  if (_primary_indexes.size() == 0) {
    return false;
  }

  _n_neutrino_candidates = _primary_indexes.size();
  std::cout << "[ElectronEventSelectionAlg] " << "Primary PFParticles " << _n_neutrino_candidates << std::endl;
  // For each of the primary particles, determine if it and it's daughters pass the cuts:

  // Need associations from pfparticle to vertex
  art::FindOneP< recob::Vertex > vertex_per_pfpart(pfparticle_handle, evt, pandoraNu_tag);


  for (auto & _i_primary : _primary_indexes ) {
    std::cout << "[ElectronEventSelectionAlg] " << "Primary PDG " << pfparticle_handle->at(_i_primary).PdgCode() << std::endl;
    std::cout << "[ElectronEventSelectionAlg] " << "N. of Daughters " << pfparticle_handle->at(_i_primary).NumDaughters() << std::endl;

    _neutrino_candidate_passed[_i_primary] = false;
    _center_of_charge[_i_primary] = TVector3(0,0,0);
    _op_flash_indexes[_i_primary] = 0;
    _neutrino_vertex[_i_primary] = TVector3(0,0,0);
    _n_showers[_i_primary] = 0;
    _pfp_id_showers_from_primary[_i_primary] = std::vector<size_t>();
    _n_tracks[_i_primary] = 0;
    _pfp_id_tracks_from_primary[_i_primary] = std::vector<size_t>();

    // First, does this event pass the optical filter?
    // Calculate the center of charge, and then pass it to the optical module:

    _center_of_charge[_i_primary]
      = calculateChargeCenter(_i_primary, pfparticle_handle, evt);

    int _selected_flash = -1;
    bool _flash_passed = opticalfilter(_i_primary,
                                       * pfparticle_handle,
                                       _center_of_charge[_i_primary],
                                       _selected_flash,
                                       evt);

    std::cout << "[ElectronEventSelectionAlg] " << "Flash passed? " << _flash_passed << std::endl;
    _op_flash_indexes[_i_primary] = _selected_flash;

    if (! _flash_passed) {
      _neutrino_candidate_passed[_i_primary] = false;
      continue;
    }


    // Get the neutrino vertex and check if it's fiducial:
    std::vector<double> neutrino_vertex;
    neutrino_vertex.resize(3);
    try {
      auto const& neutrino_vertex_obj = vertex_per_pfpart.at(_i_primary);
      neutrino_vertex_obj->XYZ(&neutrino_vertex[0]); // PFParticle neutrino vertex coordinates

      // Save it as a TVector3:
      _neutrino_vertex.at(_i_primary).SetX(neutrino_vertex[0]);
      _neutrino_vertex.at(_i_primary).SetY(neutrino_vertex[1]);
      _neutrino_vertex.at(_i_primary).SetZ(neutrino_vertex[2]);

      if (! is_fiducial(_neutrino_vertex.at(_i_primary))) {
        _neutrino_candidate_passed[_i_primary] = false;
        std::cout << "[ElectronEventSelectionAlg] " << "Neutrino vertex not within fiducial volume" << std::endl;
        continue;
      }
    } catch (...) {
      std::cout << "[ElectronEventSelectionAlg] " << "NO VERTEX AVAILABLE " << std::endl;
      return false;
    }



    // Loop over the neutrino daughters and check if there is a shower and a track
    int showers = 0;
    int tracks = 0;
    int longer_tracks = 0;

    for (auto const& pfdaughter : pfparticle_handle->at(_i_primary).Daughters())
    {
      std::cout << "[ElectronEventSelectionAlg] " << "Daughter PDG " << pfparticle_handle->at(pfdaughter).PdgCode() << std::endl;

      if (pfparticle_handle->at(pfdaughter).PdgCode() == 11)
      {
        try {
          art::FindOneP< recob::Shower > shower_per_pfpart(pfparticle_handle, evt, pandoraNu_tag);
          auto const& shower_obj = shower_per_pfpart.at(pfdaughter);

          bool contained_shower = false;
          std::vector<double> start_point;
          std::vector<double> end_point;
          start_point.resize(3);
          end_point.resize(3);

          double shower_length = shower_obj->Length();
          for (int ix = 0; ix < 3; ix++)
          {
            start_point[ix] = shower_obj->ShowerStart()[ix];
            end_point[ix] = shower_obj->ShowerStart()[ix] + shower_length * shower_obj->Direction()[ix];
          }

          contained_shower = is_fiducial(start_point) && is_fiducial(end_point);
          if (contained_shower || !contained_shower) {
            _pfp_id_showers_from_primary[_i_primary].push_back(pfdaughter);
            showers++;
          }
          // std::cout << "[ElectronEventSelectionAlg] " << "Shower energy array length " << shower_obj->Energy().size() << std::endl;
          // std::cout << "[ElectronEventSelectionAlg] " << "Shower best plane " << shower_obj->best_plane() << std::endl;
          // for (size_t i = 0; i < shower_obj->Energy().size();i++) {
          //   std::cout << "[ElectronEventSelectionAlg] " << "E " << shower_obj->Energy()[i] << std::endl;
          // }
        } catch (...) {
          std::cout << "[ElectronEventSelectionAlg] " << "NO SHOWERS AVAILABLE" << std::endl;
          return false;
        }
      }



      if (pfparticle_handle->at(pfdaughter).PdgCode() == 13)
      {
        try {

          art::FindOneP< recob::Track > track_per_pfpart(pfparticle_handle, evt, pandoraNu_tag);
          auto const& track_obj = track_per_pfpart.at(pfdaughter);

          tracks++;
          std::cout << "Tracks " << tracks << std::endl;
          std::cout << "track obj length " << track_obj->Length() << std::endl;
          std::cout << "pfdaughter " << pfdaughter << std::endl;
          std::cout << "i primary " << _i_primary << std::endl;

          _pfp_id_tracks_from_primary[_i_primary].push_back(pfdaughter);

          if (track_obj->Length() > m_trackLength) {

            longer_tracks++;
          }


        } catch (...) {
          std::cout << "[ElectronEventSelectionAlg] " << "NO TRACKS AVAILABLE" << std::endl;
          return false;
        }
        // h_track_length->Fill(track_obj->Length());
      }

      _n_tracks[_i_primary] = tracks;
      _n_showers[_i_primary] = showers;

      std::cout << "[ElectronEventSelectionAlg] " << "Showers tracks " << showers << " " << tracks << std::endl;

      if (showers >= 1 && tracks >= m_nTracks)
      {
        //closest_distance = std::min(distance(neutrino_vertex,true_neutrino_vertex),closest_distance);
        _neutrino_candidate_passed[_i_primary] = true;
      }

    } // end for pfparticle daughters


  }


  // Last, determine if any primary particles passed:
  for (auto  val : _neutrino_candidate_passed) {
    if (val.second) {
      std::cout << "[ElectronEventSelectionAlg] " << "EVENT SELECTED" << std::endl;
      return true;
    }
  }


  return false;
}


} // lee


#endif // ELECTRON_EVENT_SELECTION_ALG_CXX
