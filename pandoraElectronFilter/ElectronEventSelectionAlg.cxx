#ifndef ELECTRON_EVENT_SELECTION_ALG_CXX
#define ELECTRON_EVENT_SELECTION_ALG_CXX

#include "ElectronEventSelectionAlg.h"

namespace lee {

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

TVector3 ElectronEventSelectionAlg::spaceChargeTrueToReco(const TVector3 &xyz) {
  //auto const *sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  //geo::Point_t point(xyz);
  // auto correction = sce->GetPosOffsets(point);
  TVector3 correctedPoint(xyz);
  //correctedPoint.SetX(xyz.X() - sce->GetPosOffsets(point).X() + 0.7);
  //correctedPoint.SetX(xyz.Y() + sce->GetPosOffsets(point).Y());
 // correctedPoint.SetX(xyz.Z() + sce->GetPosOffsets(point).Z());
  return correctedPoint;
}

void ElectronEventSelectionAlg::reconfigure(fhicl::ParameterSet const &p) {
  // Implementation of optional member function here.
  m_nTracks = p.get<int>("nTracks", 1);

  m_trackLength = p.get<int>("trackLength", 100);

  m_fidvolXstart = p.get<double>("fidvolXstart", 10);
  m_fidvolXend = p.get<double>("fidvolXend", 10);

  m_fidvolYstart = p.get<double>("fidvolYstart", 20);
  m_fidvolYend = p.get<double>("fidvolYend", 20);

  m_fidvolZstart = p.get<double>("fidvolZstart", 10);
  m_fidvolZend = p.get<double>("fidvolZend", 50);

  m_isCosmicInTime = p.get<bool>("isCosmicInTime", false);

  geoHelper.setFiducialVolumeCuts(m_fidvolXstart, m_fidvolXend, m_fidvolYstart,
                                  m_fidvolYend, m_fidvolZstart, m_fidvolZend);

  m_fractionsigmaflashwidth = p.get<double>("fractionsigmaflashwidth", 2.0);
  m_absoluteflashdist = p.get<double>("absoluteflashdist", 50.0);

  fOpticalFlashFinderLabel =
      p.get<std::string>("OpticalFlashFinderLabel", "simpleFlashBeam");
}



bool ElectronEventSelectionAlg::opticalfilter(
    size_t ipf, const std::vector<recob::PFParticle> &pfparticles,
    TVector3 _this_center_of_charge, int &_selected_flash,
    const art::Event &evt) {

  bool pass = false;

  art::InputTag optical_tag{fOpticalFlashFinderLabel};
  auto const &optical_handle =
      evt.getValidHandle<std::vector<recob::OpFlash>>(optical_tag);

  double par1 =
      m_fractionsigmaflashwidth; // fraction of flash sigma that is requested
  double par2 =
      m_absoluteflashdist; // or also ok if it is closer than this number in cm

  // Run the calculation of charge center:
  int pass_index = -1;
  std::cout << "[ElectronEventSelectionAlg] " << "The optical handle size " <<
  optical_handle->size() << std::endl;

  for (unsigned int ifl = 0; ifl < optical_handle->size(); ++ifl) {
    recob::OpFlash const &flash = optical_handle->at(ifl);
    std::cout << "[ElectronEventSelectionAlg] " << flash.Time() << std::endl;
    if ((flash.Time() > 4.8 || flash.Time() < 3.2))
      continue;
    bool sigma =
        (flash.ZCenter() + flash.ZWidth() / par1) >
            _this_center_of_charge.Z() &&
        (flash.ZCenter() - flash.ZWidth() / par1) < _this_center_of_charge.Z();

    bool absolute =
        std::abs(flash.ZCenter() - _this_center_of_charge.Z()) < par2;
    std::cout << "[ElectronEventSelectionAlg] " << "The flash time is " <<
    flash.Time()
              << ", Zcentre: " << flash.ZCenter()
              << " and the Zwidth: " << flash.ZWidth()
              << std::endl;
    std::cout << "[ElectronEventSelectionAlg] " << "Z Center of charge is "
    << _this_center_of_charge.Z() << std::endl;
    if (sigma || absolute) {
      pass = true;
      pass_index = ifl;
    }
  }

  if (pass_index != -1) {
    _selected_flash = pass_index;
  }

  return pass;
}

bool ElectronEventSelectionAlg::eventSelected(const art::Event &evt) {

  clear();

  // Get the list of pfparticles:
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);

  // Are there any pfparticles?
  if (pfparticle_handle->size() == 0) {
    std::cout << "[ElectronEventSelectionAlg] "
              << "NO RECO DATA PRODUCTS" << std::endl;
    return false;
  }

  // Get the list of primary pfparticles that are also neutrinos (numu or nue)
  for (size_t _i_pfp = 0; _i_pfp < pfparticle_handle->size(); _i_pfp++) {
    if (lar_pandora::LArPandoraHelper::IsNeutrino(art::Ptr<recob::PFParticle>(pfparticle_handle, _i_pfp))) {
      _primary_indexes.push_back(_i_pfp);
    }
  }

  // If there are no particles flagged as primary, return false
  if (_primary_indexes.size() == 0) {
    return false;
  }

  _n_neutrino_candidates = _primary_indexes.size();
  std::cout << "[ElectronEventSelectionAlg] "
            << "Primary PFParticles " << _n_neutrino_candidates << std::endl;
  // For each of the primary particles, determine if it and it's daughters pass
  // the cuts:

  // Need associations from pfparticle to vertex
  art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt,
                                                 _pfp_producer);

  for (auto &_i_primary : _primary_indexes) {
    std::cout << "[ElectronEventSelectionAlg] "
              << "Primary PDG " << pfparticle_handle->at(_i_primary).PdgCode()
              << std::endl;
    std::cout << "[ElectronEventSelectionAlg] "
              << "N. of Daughters "
              << pfparticle_handle->at(_i_primary).NumDaughters() << std::endl;

    _neutrino_candidate_passed[_i_primary] = false;
    _center_of_charge[_i_primary] = TVector3(0, 0, 0);
    _op_flash_indexes[_i_primary] = 0;
    _neutrino_vertex[_i_primary] = TVector3(0, 0, 0);
    _n_showers[_i_primary] = 0;
    _pfp_id_showers_from_primary[_i_primary] = std::vector<size_t>();
    _n_tracks[_i_primary] = 0;
    _pfp_id_tracks_from_primary[_i_primary] = std::vector<size_t>();

    // First, does this event pass the optical filter?
    // Calculate the center of charge, and then pass it to the optical module:

    _center_of_charge[_i_primary] =
        pandoraHelper.calculateChargeCenter(_i_primary, pfparticle_handle, evt);

    int _selected_flash = -1;
    bool _flash_passed =
        opticalfilter(_i_primary, *pfparticle_handle,
                      _center_of_charge[_i_primary], _selected_flash, evt);

    std::cout << "[ElectronEventSelectionAlg] "
              << "Flash passed? " << _flash_passed << " " << _selected_flash << std::endl;

    _op_flash_indexes[_i_primary] = _selected_flash;

    if (!_flash_passed) {
      _neutrino_candidate_passed[_i_primary] = false;
      continue;
    }

    // Get the neutrino vertex and check if it's fiducial:
    std::vector<double> neutrino_vertex;
    neutrino_vertex.resize(3);
    try {
      auto const &neutrino_vertex_obj = vertex_per_pfpart.at(_i_primary);
      neutrino_vertex_obj->XYZ(
          &neutrino_vertex[0]); // PFParticle neutrino vertex coordinates

      // Save it as a TVector3:
      _neutrino_vertex.at(_i_primary).SetX(neutrino_vertex[0]);
      _neutrino_vertex.at(_i_primary).SetY(neutrino_vertex[1]);
      _neutrino_vertex.at(_i_primary).SetZ(neutrino_vertex[2]);

      if (!geoHelper.isFiducial(_neutrino_vertex.at(_i_primary))) {
        _neutrino_candidate_passed[_i_primary] = false;
        std::cout << "[ElectronEventSelectionAlg] "
                  << "Neutrino vertex not within fiducial volume" << std::endl;
        continue;
      }
    } catch (...) {
      std::cout << "[ElectronEventSelectionAlg] "
                << "NO VERTEX AVAILABLE " << std::endl;
      _neutrino_candidate_passed[_i_primary] = false;
      continue;
    }

    // Loop over the neutrino daughters and check if there is a shower and a
    // track
    int showers = 0;
    int tracks = 0;

    for (auto const &pfdaughter :
         pfparticle_handle->at(_i_primary).Daughters()) {
      std::cout << "[ElectronEventSelectionAlg] "
                << "Daughter PDG "
                << pfparticle_handle->at(pfdaughter).PdgCode() << std::endl;

      if (pfparticle_handle->at(pfdaughter).PdgCode() == 11) {
        try {
          art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt,
                                                         _pfp_producer);
          auto const &shower_obj = shower_per_pfpart.at(pfdaughter);

          bool contained_shower = false;
          std::vector<double> start_point;
          std::vector<double> end_point;
          start_point.resize(3);
          end_point.resize(3);

          double shower_length = shower_obj->Length();
          for (int ix = 0; ix < 3; ix++) {
            start_point[ix] = shower_obj->ShowerStart()[ix];
            end_point[ix] = shower_obj->ShowerStart()[ix] +
                            shower_length * shower_obj->Direction()[ix];
          }

          contained_shower = geoHelper.isFiducial(start_point) && geoHelper.isFiducial(end_point);
          if (contained_shower || !contained_shower) {
            _pfp_id_showers_from_primary[_i_primary].push_back(pfdaughter);
            showers++;
          }
          // std::cout << "[ElectronEventSelectionAlg] " << "Shower energy array
          // length " << shower_obj->Energy().size() << std::endl;
          // std::cout << "[ElectronEventSelectionAlg] " << "Shower best plane "
          // << shower_obj->best_plane() << std::endl;
          // for (size_t i = 0; i < shower_obj->Energy().size();i++) {
          //   std::cout << "[ElectronEventSelectionAlg] " << "E " <<
          //   shower_obj->Energy()[i] << std::endl;
          // }
        } catch (...) {
          std::cout << "[ElectronEventSelectionAlg] "
                    << "NO SHOWERS AVAILABLE" << std::endl;
        }
      }

      if (pfparticle_handle->at(pfdaughter).PdgCode() == 13) {
        try {

          art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt,
                                                       _pfp_producer);
          auto const &track_obj = track_per_pfpart.at(pfdaughter);

          std::cout << "track obj length " << track_obj->Length() << std::endl;
          std::cout << "pfdaughter " << pfdaughter << std::endl;
          std::cout << "i primary " << _i_primary << std::endl;

          _pfp_id_tracks_from_primary[_i_primary].push_back(pfdaughter);


          tracks++;
          std::cout << "Tracks " << tracks << std::endl;


        } catch (...) {
          std::cout << "[ElectronEventSelectionAlg] "
                    << "NO TRACKS AVAILABLE" << std::endl;
        }
        // h_track_length->Fill(track_obj->Length());
      }

      _n_tracks[_i_primary] = tracks;
      _n_showers[_i_primary] = showers;

      std::cout << "[ElectronEventSelectionAlg] "
                << "Showers tracks " << showers << " " << tracks << std::endl;

      if (showers >= 1 && tracks >= m_nTracks) {
        // closest_distance =
        // std::min(distance(neutrino_vertex,true_neutrino_vertex),closest_distance);
        _neutrino_candidate_passed[_i_primary] = true;
      }

    } // end for pfparticle daughters
  }

  // Last, determine if any primary particles passed:
  for (auto val : _neutrino_candidate_passed) {
    if (val.second) {
      std::cout << "[ElectronEventSelectionAlg] "
                << "EVENT SELECTED" << std::endl;
      return true;
    }
  }



  return false;
}

} // lee

#endif // ELECTRON_EVENT_SELECTION_ALG_CXX
