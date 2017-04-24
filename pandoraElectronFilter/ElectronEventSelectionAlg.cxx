#ifndef ELECTRON_EVENT_SELECTION_ALG_CXX
#define ELECTRON_EVENT_SELECTION_ALG_CXX

#include "ElectronEventSelectionAlg.h"

namespace lee {

void ElectronEventSelectionAlg::clear() {
  _n_neutrino_candidates = 0.0;
  _neutrino_candidate_passed.clear();
  _center_of_charge.clear();
  _op_flash_indexes.clear();
  _neutrino_vertex.clear();
  _n_showers.clear();
  _n_tracks.clear();
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
  std::vector<size_t> unordered_daugthers )
{

  // This is a tree-traversal algorithm.  It returns the index of the top particle, plus the index
  // of all daughter particles.

  // This is a recursive algorithm, so it needs a break clause:

  if (pfparticles->at(top_index).Daughters().size() == 0) {
    unordered_daugthers.push_back(top_index);
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
      std::vector<art::Ptr<recob::Hit> > hits = hits_per_spcpnts.at(_sps.key());
      // Add the hits to the weighted average, if they are collection hits:
      for (auto & hit : hits) {
        if (hit->View() == geo::kZ) {
          // Collection hits only
          double weight = hit->Integral();
          auto xyz = _sps->XYZ();
          chargecenter[0] += (xyz[0]) * weight;
          chargecenter[1] += (xyz[1]) * weight;
          chargecenter[2] += (xyz[2]) * weight;
          totalweight += weight;
          break; // Exit the loop over hits
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
  size_t _selected_flash,
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
    bool sigma    = flash.ZCenter() + flash.ZWidth() / par1 > _this_center_of_charge.Z() && flash.ZCenter() + flash.ZWidth() / par1 < _this_center_of_charge.Z();
    bool absolute = std::abs(flash.ZCenter() - _this_center_of_charge.Z()) < par2;
    //std::cout << "The flash time is " << flash.Time() << ", Zcentre: " << flash.ZCenter() << " and the Zwidth: " << flash.ZWidth() << std::endl;
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
  bool pass = false;
  return pass;

  art::InputTag pandoraNu_tag { "pandoraNu" };

  // Get the list of pfparticles:
  auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( pandoraNu_tag );

  // Are there any pfparticles?
  if (pfparticle_handle -> size() == 0) {
    std::cout << "NO RECO DATA PRODUCTS" << std::endl;
    return false;
  }

  // Get the list of primary pfparticles that are also neutrinos (numu or nue)
  std::vector<size_t> _primary_indexes;
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
  _neutrino_candidate_passed.resize(_n_neutrino_candidates);
  _center_of_charge.resize(_n_neutrino_candidates);
  _op_flash_indexes.resize(_n_neutrino_candidates);
  _neutrino_vertex.resize(_n_neutrino_candidates);
  _n_showers.resize(_n_neutrino_candidates);
  _n_tracks.resize(_n_neutrino_candidates);

  // For each of the primary particles, determine if it and it's daughters pass the cuts:

  // Need associations from pfparticle to vertex
  art::FindOneP< recob::Vertex > vertex_per_pfpart(pfparticle_handle, evt, pandoraNu_tag);


  for (size_t _i_primary = 0; _i_primary < _primary_indexes.size(); _i_primary ++ ) {

    // First, does this event pass the optical filter?
    // Calculate the center of charge, and then pass it to the optical module:

    _center_of_charge[_i_primary]
      = calculateChargeCenter(_primary_indexes.at(_i_primary), pfparticle_handle, evt);

    int _selected_flash = -1;
    bool _flash_passed = opticalfilter(_primary_indexes.at(_i_primary),
                                       * pfparticle_handle,
                                       _center_of_charge[_i_primary],
                                       _selected_flash,
                                       evt);

    if (! _flash_passed) {
      _neutrino_candidate_passed[_i_primary] = false;
      continue;
    }

    // Get the neutrino vertex and check if it's fiducial:
    std::vector<double> neutrino_vertex;
    neutrino_vertex.resize(3);
    auto const& neutrino_vertex_obj = vertex_per_pfpart.at(_primary_indexes[_i_primary]);
    neutrino_vertex_obj->XYZ(&neutrino_vertex[0]); // PFParticle neutrino vertex coordinates

    // Save it as a TVector3:
    _neutrino_vertex.at(_i_primary).SetX(neutrino_vertex[0]);
    _neutrino_vertex.at(_i_primary).SetY(neutrino_vertex[1]);
    _neutrino_vertex.at(_i_primary).SetZ(neutrino_vertex[2]);

    if (! is_fiducial(_neutrino_vertex.at(_i_primary))) {
      _neutrino_candidate_passed[_i_primary] = false;
      continue;
    }



    // Loop over the neutrino daughters and check if there is a shower and a track
    int showers = 0;
    int tracks = 0;
    for (auto const& pfdaughter : pfparticle_handle->at(_primary_indexes[_i_primary]).Daughters())
    {


      if (pfparticle_handle->at(pfdaughter).PdgCode() == 11)
      {
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
        // TODO flash position check
        if (contained_shower) showers++;

      }

      if (pfparticle_handle->at(pfdaughter).PdgCode() == 13)
      {
        art::FindOneP< recob::Track > track_per_pfpart(pfparticle_handle, evt, pandoraNu_tag);
        auto const& track_obj = track_per_pfpart.at(pfdaughter);

        if (track_obj->Length() < m_trackLength) tracks++;
        // h_track_length->Fill(track_obj->Length());
      }

      _n_tracks[_i_primary] = tracks;
      _n_showers[_i_primary] = showers;

      if (showers >= 1 && tracks >= m_nTracks)
      {
        //closest_distance = std::min(distance(neutrino_vertex,true_neutrino_vertex),closest_distance);
        _neutrino_candidate_passed[_i_primary] = true;
      }

    } // end for pfparticle daughters


  }


  // Last, determine if any primary particles passed:
  for (auto  val : _neutrino_candidate_passed) {
    if (val) {
      std::cout << "EVENT SELECTED" << std::endl;
      return true;
    }
  }
  return false;
}


} // lee


#endif // ELECTRON_EVENT_SELECTION_ALG_CXX
