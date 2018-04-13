#ifndef FEATUREHELPER_CXX
#define FEATUREHELPER_CXX

#include "FeatureHelper.h"

namespace lee
{

FeatureHelper::FeatureHelper()
{
  // UGLY AND DANGEROUS!
  geoHelper.setFiducialVolumeCuts(10, 10, 20, 20, 10, 50);

  //detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

/////////////////// HELP FUNCTIONS ///////////////////

// Calculate the dot product between two arrays in 3d
double dot_product(double dir1[], double dir2[])
{
  std::cout << "dot product called" << std::endl;
  double nominator = (dir1[0] * dir2[0]) + (dir1[1] * dir2[1]) + (dir1[2] * dir2[2]);
  double denominator = sqrt((dir1[0] * dir1[0]) + (dir1[1] * dir1[1]) + (dir1[2] * dir1[2])) * sqrt((dir2[0] * dir2[0]) + (dir2[1] * dir2[1]) + (dir2[2] * dir2[2]));
  return nominator / denominator;
}

float dot_product2(float dir1[], float dir2[])
{
  //std::cout << "dot product 2 called" << std::endl;
  float nominator = (dir1[0] * dir2[0]) + (dir1[1] * dir2[1]) + (dir1[2] * dir2[2]);
  float denominator = sqrt((dir1[0] * dir1[0]) + (dir1[1] * dir1[1]) + (dir1[2] * dir1[2]));
  return nominator / denominator;
}

// Method to return the pdgcode the parent is matched to given pf_id of the daughter
int daughterMatchHelper(const art::Event &evt, std::string _pfp_producer,
                        std::vector<size_t> &_nu_shower_ids, std::vector<size_t> &_nu_track_ids,
                        std::vector<int> &_matched_showers, std::vector<int> &_matched_tracks,
                        int pf_id_daughter)
{
  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  size_t id_parent = pfparticle_handle->at(pf_id_daughter).Parent();
  recob::PFParticle const &pfparent = pfparticle_handle->at(id_parent);
  size_t pdg_parent = pfparent.PdgCode();
  if (pdg_parent == 11)
  {
    std::cout << "[MatchDebug] Particle is daughter of shower with id " << id_parent << std::endl;
    for (size_t i = 0; i < _nu_shower_ids.size(); ++i)
    {
      if (_nu_shower_ids[i] == id_parent)
      {
        return _matched_showers[i];
      }
    }
  }
  else if (pdg_parent == 13)
  {
    std::cout << "[MatchDebug] Particle is daughter of track with id " << id_parent << std::endl;
    for (size_t i = 0; i < _nu_track_ids.size(); ++i)
    {
      if (_nu_track_ids[i] == id_parent)
      {
        return _matched_tracks[i];
      }
    }
  }
  else
  {
    std::cout << "[MatchDebug] Particle is PRIMARY, HELP!" << std::endl;
    return 0;
  }
  std::cout << "[MatchDebug] Particle is PRIMARY, Strange Help!" << std::endl;
  return 0;
}

bool isContained(double x, double y, double z, double border)
{
  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {
      0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
      0., geo->DetLength()};

  bool is_x = x > (bnd[0] + border) && x < (bnd[1] - border);
  bool is_y = y > (bnd[2] + border) && y < (bnd[3] - border);
  bool is_z = z > (bnd[4] + border) && z < (bnd[5] - border);
  return is_x && is_y && is_z;
}

//////////////////////////////////////////////////////

int FeatureHelper::true_thresholds_1eX(int &_true_nu_is_fiducial,
                                       std::vector<int> &_nu_daughters_pdg, std::vector<double> &_nu_daughters_E)
{
  if (_true_nu_is_fiducial)
  {
    for (uint i = 0; i < _nu_daughters_pdg.size(); ++i)
    {
      if (_nu_daughters_pdg[i] == 11)
      {
        if (_nu_daughters_E[i] > min_e_kinE)
        {
          return 1;
        }
      }
    }
  }
  return 0;
}

void FeatureHelper::true_closest_electron_matched(std::vector<int> &_matched_showers, std::vector<int> &_matched_tracks,
                                                  double &_true_vx_sce, double &_true_vy_sce, double &_true_vz_sce,
                                                  std::vector<double> &_shower_start_x, std::vector<double> &_shower_start_y, std::vector<double> &_shower_start_z,
                                                  std::vector<double> &_track_start_x, std::vector<double> &_track_start_y, std::vector<double> &_track_start_z,
                                                  std::vector<int> &_shower_cle, std::vector<int> &_track_cle)
{
  _shower_cle.resize(_matched_showers.size(), 0);
  _track_cle.resize(_matched_tracks.size(), 0);
  // If there is no electron matched object closer than 10cm from the true vertex, than do not consisder any object.
  float dist = 10;
  int sh_i = -1;
  int tr_i = -1;
  if (std::find(_matched_tracks.begin(), _matched_tracks.end(), 11) != _matched_tracks.end())
  {
    std::vector<double> vtxdist = reco_vtxdistance(_true_vx_sce, _true_vy_sce, _true_vz_sce, _track_start_x, _track_start_y, _track_start_z);
    /* _matched_tracks contains 11 */
    for (uint i = 0; i < _matched_tracks.size(); ++i)
    {
      if (_matched_tracks[i] == 11 && dist > vtxdist[i])
      {
        dist = vtxdist[i];
        tr_i = i;
      }
    }
  }
  if (std::find(_matched_showers.begin(), _matched_showers.end(), 11) != _matched_showers.end())
  {
    std::vector<double> vtxdist = reco_vtxdistance(_true_vx_sce, _true_vy_sce, _true_vz_sce, _shower_start_x, _shower_start_y, _shower_start_z);
    /* _matched_showers contains 11 */
    for (uint i = 0; i < _matched_showers.size(); ++i)
    {
      if (_matched_showers[i] == 11 && dist > vtxdist[i])
      {
        dist = vtxdist[i];
        sh_i = i;
        tr_i = -1;
      }
    }
  }
  if (sh_i != -1)
  {
    _shower_cle[sh_i] = 1;
  }
  else if (tr_i != -1)
  {
    _track_cle[tr_i] = 1;
  }
}

int FeatureHelper::reco_bdt_track_precut(std::vector<double> &_predict_mu, std::vector<double> &_predict_cos, int &_n_tracks)
{
  for (int i = 0; i < _n_tracks; ++i)
  {
    if (_predict_mu[i] > max_muon_score || _predict_cos[i] > max_cosmic_score)
    {
      return 0;
    }
  }
  return 1;
}

// Method to calculate the centre of charge position for vector of pf_ids
float FeatureHelper::reco_calculateChargeCenter(size_t &object_id, const art::Event &evt, std::string _pfp_producer,
                                                float chargecenter[])
{
  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);
  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);

  float totalweight = 0;

  // Get the associated spacepoints:
  std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(object_id);

  // Loop over the spacepoints and get the associated hits:
  for (auto &_sps : spcpnts)
  {
    auto xyz = _sps->XYZ();
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
    // Add the hits to the weighted average, if they are collection hits:
    for (auto &hit : hits)
    {
      if (hit->View() == geo::kZ)
      {
        // Collection hits only
        float weight = hit->Integral();
        chargecenter[0] += (xyz[0]) * weight;
        chargecenter[1] += (xyz[1]) * weight;
        chargecenter[2] += (xyz[2]) * weight;
        totalweight += weight;
      } // if collection
    }   // hits
  }     // spacepoints
  return totalweight;
}

void FeatureHelper::reco_totalChargeCenter(std::vector<size_t> &shower_ids, std::vector<size_t> &track_ids, const art::Event &evt, std::string _pfp_producer,
                                           float &chargecenter_x, float &chargecenter_y, float &chargecenter_z)
{
  float chargecenter[3] = {0, 0, 0};
  float totalweight = 0;

  for (size_t shower_id : shower_ids)
  {
    // reco_calculateChargeCenter does not reset the chargecenter
    totalweight += reco_calculateChargeCenter(shower_id, evt, _pfp_producer, chargecenter);
  }
  for (size_t track_id : track_ids)
  {
    // reco_calculateChargeCenter does not reset the chargecenter
    totalweight += reco_calculateChargeCenter(track_id, evt, _pfp_producer, chargecenter);
  }

  chargecenter_x = chargecenter[0] / totalweight;
  chargecenter_y = chargecenter[1] / totalweight;
  chargecenter_z = chargecenter[2] / totalweight;
}

float FeatureHelper::reco_total_spacepoint_containment(std::vector<size_t> &shower_ids, std::vector<size_t> &track_ids, const art::Event &evt, std::string _pfp_producer)
{
  float in_fidvol_q_total = 0;
  float out_fidvol_q_total = 0;
  float in_fidvol_q;
  float out_fidvol_q;
  for (size_t shower_id : shower_ids)
  {
    // reco_spacepoint_containment does reset the counters
    reco_spacepoint_containment(shower_id, evt, _pfp_producer, in_fidvol_q, out_fidvol_q);
    in_fidvol_q_total += in_fidvol_q;
    out_fidvol_q_total += out_fidvol_q;
  }
  for (size_t track_id : track_ids)
  {
    // reco_spacepoint_containment does reset the counters
    reco_spacepoint_containment(track_id, evt, _pfp_producer, in_fidvol_q, out_fidvol_q);
    in_fidvol_q_total += in_fidvol_q;
    out_fidvol_q_total += out_fidvol_q;
  }
  if (in_fidvol_q_total > 0)
  {
    return in_fidvol_q_total / (in_fidvol_q_total + out_fidvol_q_total);
  }
  else
  {
    return -1;
  }
}

void FeatureHelper::reco_maxangle(std::vector<double> &_shower_dir_x, std::vector<double> &_shower_dir_y, std::vector<double> &_shower_dir_z,
                                  std::vector<double> &_track_dir_x, std::vector<double> &_track_dir_y, std::vector<double> &_track_dir_z,
                                  std::vector<double> &_track_maxangle, std::vector<double> &_shower_maxangle)
{
  uint n_tracks = _track_dir_x.size();
  uint n_showers = _shower_dir_x.size();
  _track_maxangle.resize(n_tracks);
  _shower_maxangle.resize(n_showers);

  double cosine = 1;
  double temp_cos = 1;

  //std::cout << "n_tracks " << n_tracks << " n_showers " << n_showers << std::endl;
  // First for tracks
  for (unsigned int i = 0; i < n_tracks; ++i)
  {
    double dir1[3] = {_track_dir_x[i], _track_dir_y[i], _track_dir_z[i]};
    for (unsigned int j = 0; j < n_tracks; ++j)
    {
      std::cout << i << " " << j << std::endl;
      double dir2[3] = {_track_dir_x[j], _track_dir_y[j], _track_dir_z[j]};
      temp_cos = dot_product(dir1, dir2);
      if (temp_cos < cosine)
      {
        //std::cout << "track_maxangle: temp_cos_track" << temp_cos << std::endl;
        cosine = temp_cos;
      }
    }
    for (unsigned int j = 0; j < n_showers; ++j)
    {
      double dir2[3] = {_shower_dir_x[j], _shower_dir_y[j], _shower_dir_z[j]};
      temp_cos = dot_product(dir1, dir2);
      if (temp_cos < cosine)
      {
        //std::cout << "track_maxangle: temp_cos_shower" << temp_cos << std::endl;
        cosine = temp_cos;
      }
    }
    //std::cout << "track_maxangle: " << cosine << std::endl;
    _track_maxangle[i] = cosine;
  }

  //Now the same for showers
  for (unsigned int i = 0; i < n_showers; ++i)
  {
    double dir1[3] = {_shower_dir_x[i], _shower_dir_y[i], _shower_dir_z[i]};
    for (unsigned int j = 0; j < n_tracks; ++j)
    {
      double dir2[3] = {_track_dir_x[j], _track_dir_y[j], _track_dir_z[j]};
      temp_cos = dot_product(dir1, dir2);
      if (temp_cos < cosine)
      {
        //std::cout << "shower_maxangle: temp_cos_track" << temp_cos << std::endl;
        cosine = temp_cos;
      }
    }
    for (unsigned int j = 0; j < n_showers; ++j)
    {
      double dir2[3] = {_shower_dir_x[j], _shower_dir_y[j], _shower_dir_z[j]};
      temp_cos = dot_product(dir1, dir2);
      if (temp_cos < cosine)
      {
        //std::cout << "shower_maxangle: temp_cos_shower" << temp_cos << std::endl;
        cosine = temp_cos;
      }
    }
    //std::cout << "shower_maxangle: " << cosine << std::endl;
    _shower_maxangle[i] = cosine;
  }
}

void FeatureHelper::reco_hierarchy(std::vector<size_t> &_nu_object_ids, const art::Event &evt, std::string _pfp_producer,
                                   std::vector<int> &object_daughter, std::vector<int> &object_is_daughter)
{
  uint n_objects = _nu_object_ids.size();
  object_daughter.resize(n_objects);
  object_is_daughter.resize(n_objects);
  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);

  for (uint i = 0; i < n_objects; ++i)
  {
    size_t pf_id = _nu_object_ids[i];
    recob::PFParticle const &pfparticle = pfparticle_handle->at(pf_id);

    int parent_pdg = pfparticle_handle->at(pfparticle.Parent()).PdgCode();
    const std::vector<size_t> daughter_ids = pfparticle.Daughters();
    if (parent_pdg == 11)
    {
      object_is_daughter[i] = 1;
    }
    else if (parent_pdg == 13)
    {
      object_is_daughter[i] = 2;
    }
    else
    {
      object_is_daughter[i] = 0;
    }

    if (daughter_ids.size() == 0)
    {
      object_daughter[i] = 0;
    }
    else if (daughter_ids.size() > 1)
    {
      object_daughter[i] = 3;
    }
    else
    {
      int daughter_pdg = pfparticle_handle->at(daughter_ids[0]).PdgCode();
      if (daughter_pdg == 11)
      {
        object_daughter[i] = 1;
      }
      else if (daughter_pdg == 13)
      {
        object_daughter[i] = 2;
      }
      else
      {
        object_daughter[i] = 0;
        std::cout << "[FeatureHelper] "
                  << "Object has daughter which is not a track/shower. ERROR!" << std::endl;
      }
    }
  }
}

void FeatureHelper::reco_spacepoint_containment(size_t pf_id, const art::Event &evt, std::string _pfp_producer,
                                                float &in_fidvol_q, float &out_fidvol_q)
{
  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);
  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);
  std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(pf_id);

  in_fidvol_q = 0;
  out_fidvol_q = 0;
  bool sps_is_fiducial = false;

  for (auto &_sps : spcpnts)
  {
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
    auto xyz = _sps->XYZ();
    sps_is_fiducial = geoHelper.isFiducial(xyz);

    for (auto &hit : hits)
    {
      if (hit->View() == geo::kZ)
      {
        float integral = (float)hit->Integral();
        if (sps_is_fiducial)
        {
          in_fidvol_q += integral;
        }
        else
        {
          out_fidvol_q += integral;
        }
      }
    }
  }
}

float FeatureHelper::reco_spacepoint_dqdx(size_t object_id, const art::Event &evt, std::string _pfp_producer,
                                          std::vector<float> &sps_dqdx_distance, std::vector<float> &sps_dqdx_integral)
{
  float chargecenter[3] = {0, 0, 0};
  float totalweight = reco_calculateChargeCenter(object_id, evt, _pfp_producer, chargecenter);
  if (totalweight > 0)
  {
    chargecenter[0] /= totalweight;
    chargecenter[1] /= totalweight;
    chargecenter[2] /= totalweight;
    auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
    auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);
    art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);
    art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, _pfp_producer);

    auto const &vertex_obj = vertex_per_pfpart.at(object_id);
    double vertex[3];
    vertex_obj->XYZ(vertex);

    float norm[3] = {chargecenter[0] - (float)vertex[0],
                     chargecenter[1] - (float)vertex[1],
                     chargecenter[2] - (float)vertex[2]};

    sps_dqdx_distance.clear();
    sps_dqdx_integral.clear();

    // Get the associated spacepoints:
    std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(object_id);

    float maxdist = 0;
    // Loop over the spacepoints and get the associated hits:
    for (auto &_sps : spcpnts)
    {
      float totalint = 0;
      float xyz[3] = {(float)_sps->XYZ()[0] - (float)vertex[0],
                      (float)_sps->XYZ()[1] - (float)vertex[1],
                      (float)_sps->XYZ()[2] - (float)vertex[2]};
      std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
      // Add the hits to the weighted average, if they are collection hits:
      for (auto &hit : hits)
      {
        if (hit->View() == geo::kZ)
        {
          // Collection hits only
          totalint += hit->Integral();
        } // if collection
      }   // hits
      if (totalint == 0)
      {
        // no collection charge in this point so do not save it.
        continue;
      }

      float this_distance = dot_product2(norm, xyz);
      //std::cout << "Spacepoint at distance " << this_distance << " totalint " << totalint << std::endl;
      sps_dqdx_distance.push_back(this_distance);
      sps_dqdx_integral.push_back(totalint);
      if (this_distance > maxdist)
      {
        maxdist = this_distance;
      }
    } // spacepoints
    return maxdist;
  }
  else
  {
    return -1;
  }
}

void FeatureHelper::reco_spacepoint_ratios(std::vector<size_t> _nu_object_ids, const art::Event &evt, std::string _pfp_producer,
                                           std::vector<float> &object_containment, std::vector<float> &object_dqdx_ratio)
{
  object_containment.clear();
  object_dqdx_ratio.clear();

  for (size_t pf_id : _nu_object_ids)
  {
    // For the containment:
    float in_fidvol_q_this;
    float out_fidvol_q_this;
    reco_spacepoint_containment(pf_id, evt, _pfp_producer, in_fidvol_q_this, out_fidvol_q_this);
    if (in_fidvol_q_this > 0)
    {
      std::cout << "object_containment: " << in_fidvol_q_this / (in_fidvol_q_this + out_fidvol_q_this) << std::endl;
      object_containment.push_back(in_fidvol_q_this / (in_fidvol_q_this + out_fidvol_q_this));
    }
    else
    {
      std::cout << "[FeatureHelper] PFP_id: " << pf_id << " had no spacepoint collection charge inside fidvol." << std::endl;
      object_containment.push_back(0.0);
    }

    // For the dqdx ratio:
    std::vector<float> sps_dqdx_distance;
    std::vector<float> sps_dqdx_integral;
    float maxdist = reco_spacepoint_dqdx(pf_id, evt, _pfp_producer, sps_dqdx_distance, sps_dqdx_integral);
    std::cout << "maxdist " << maxdist << std::endl;
    if (maxdist == -1)
    {
      object_dqdx_ratio.push_back(0.0);
      std::cout << "[FeatureHelper] PFP_id: " << pf_id << " had no spacepoint collection charge." << std::endl;
    }
    else
    {
      size_t n_points = sps_dqdx_distance.size();
      float part_1 = 0;
      float part_2 = 0;
      for (size_t i = 0; i < n_points; ++i)
      {
        if (sps_dqdx_distance[i] < (maxdist / 2))
        {
          part_1 += sps_dqdx_integral[i];
        }
        else
        {
          part_2 += sps_dqdx_integral[i];
        }
      }
      std::cout << "object_dqdx_ratio: part 1 " << part_1 << " part 2 " << part_2 << " ratio " << part_1 / part_2 << std::endl;
      object_dqdx_ratio.push_back(part_1 / part_2);
    }
  }
}

void FeatureHelper::reco_flash_info(std::vector<int> &_flash_passed, std::vector<double> &_flash_PE, std::vector<double> &_flash_time,
                                    double &_flash_PE_max, double &_flash_time_max)
{
  for (const int fl : _flash_passed)
  {
    if (fl != -1)
    {
      _flash_PE_max = _flash_PE[fl];
      _flash_time_max = _flash_time[fl];
    }
  }
}

void FeatureHelper::true_match_daughters(const art::Event &evt, std::string _pfp_producer,
                                         std::vector<size_t> &_nu_shower_ids, std::vector<size_t> &_nu_track_ids,
                                         std::vector<int> &_matched_showers, std::vector<int> &_matched_tracks)
{
  // First for showers
  for (uint ish = 0; ish < _matched_showers.size(); ++ish)
  {
    if (_matched_showers[ish] == std::numeric_limits<int>::lowest())
    {
      _matched_showers[ish] = daughterMatchHelper(evt, _pfp_producer,
                                                  _nu_shower_ids, _nu_track_ids,
                                                  _matched_showers, _matched_tracks,
                                                  _nu_shower_ids[ish]);
    }
  }
  // Now for tracks
  for (uint itr = 0; itr < _matched_tracks.size(); ++itr)
  {
    if (_matched_tracks[itr] == std::numeric_limits<int>::lowest())
    {
      _matched_tracks[itr] = daughterMatchHelper(evt, _pfp_producer,
                                                 _nu_shower_ids, _nu_track_ids,
                                                 _matched_showers, _matched_tracks,
                                                 _nu_track_ids[itr]);
    }
  }
}

void FeatureHelper::reco_track_containment(std::vector<double> &_track_end_x, std::vector<double> &_track_end_y, std::vector<double> &_track_end_z,
                                           std::vector<int> &_track_containment)
{
  uint n_tracks = _track_end_x.size();
  _track_containment.resize(n_tracks);
  for (uint i = 0; i < n_tracks; ++i)
  {
    _track_containment[i] = isContained(_track_end_x[i], _track_end_y[i], _track_end_z[i], track_containment_border);
  }
}

std::vector<double> FeatureHelper::reco_vtxdistance(double &vx, double &vy, double &vz,
                                                    std::vector<double> &start_x, std::vector<double> &start_y, std::vector<double> &start_z)
{
  uint n_objects = start_x.size();
  std::vector<double> object_dist(n_objects);
  for (uint i = 0; i < n_objects; ++i)
  {
    object_dist[i] = sqrt(((vx - start_x[i]) * (vx - start_x[i])) + ((vy - start_y[i]) * (vy - start_y[i])) + ((vz - start_z[i]) * (vz - start_z[i])));
  }
  return object_dist;
}

void FeatureHelper::reco_dedx(std::vector<std::vector<double>> &_object_dEdx_hits,
                              std::vector<std::vector<double>> &_object_dEdx,
                              std::vector<std::vector<float>> &_object_dQdx_cali,
                              std::vector<int> &_object_dedx_hits_w,
                              std::vector<float> &_object_dedx_w,
                              std::vector<float> &_object_dedx_best_w)
{
  uint n_objects = _object_dEdx.size();
  _object_dedx_hits_w.resize(n_objects);
  _object_dedx_w.resize(n_objects);
  _object_dedx_best_w.resize(n_objects);

  for (uint i = 0; i < n_objects; ++i)
  {
    _object_dedx_hits_w[i] = _object_dEdx_hits[i].size();
    _object_dedx_w[i] = _object_dEdx[i][2] * _object_dQdx_cali[i][2]; // collection plane

    double dedx_avg_collection = 0;
    // Protect against no cluster hits on collection plane:
    if (_object_dedx_hits_w[i] > 0)
    {
      dedx_avg_collection = std::accumulate(_object_dEdx_hits[i].begin(), _object_dEdx_hits[i].end(), 0.0) / _object_dedx_hits_w[i] * _object_dQdx_cali[i][2];
    }

    double dedx_choices[4] = {_object_dEdx[i][0] * _object_dQdx_cali[i][0] - target_electron_dedx,
                              _object_dEdx[i][1] * _object_dQdx_cali[i][1] - target_electron_dedx,
                              _object_dEdx[i][2] * _object_dQdx_cali[i][2] - target_electron_dedx,
                              dedx_avg_collection - target_electron_dedx};
    double min = dedx_choices[0];
    for (uint j = 1; j < 4; ++j)
    {
      if (abs(min) > abs(dedx_choices[j]))
      {
        min = dedx_choices[j];
        std::cout << "dedx_choices[j]" << dedx_choices[j] << std::endl;
      }
    }
    _object_dedx_best_w[i] = (float)min;
  }
}

void FeatureHelper::reco_energy(std::vector<std::vector<double>> &_object_energy_hits,
                                std::vector<std::vector<float>> &_object_energy_cali,
                                std::vector<std::vector<int>> &_object_nhits_cluster,
                                std::vector<std::vector<int>> &_object_nhits_spacepoint,
                                std::vector<float> &_object_energy_w,
                                std::vector<float> &_object_hitsratio_w,
                                std::vector<int> &_object_hits_w)
{
  uint n_objects = _object_energy_hits.size();
  _object_energy_w.resize(n_objects);
  _object_hitsratio_w.resize(n_objects);
  _object_hits_w.resize(n_objects);
  for (uint i = 0; i < n_objects; ++i)
  {
    _object_energy_w[i] = _object_energy_hits[i][2] * _object_energy_cali[i][2];
    _object_hits_w[i] = _object_nhits_cluster[i][2];

    _object_hitsratio_w[i] = 0;
    // Protect against no cluster hits on collection plane:
    if ( _object_hits_w[i] > 0)
    {
      _object_hitsratio_w[i] = _object_nhits_spacepoint[i][2] / _object_nhits_cluster[i][2];
    }
  }
}

} // namespace lee

#endif
