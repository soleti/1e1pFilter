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
  double nominator = (dir1[0] * dir2[0]) + (dir1[1] * dir2[1]) + (dir1[2] * dir2[2]);
  double denominator = sqrt(dir1[0] * dir1[0]) + (dir1[1] * dir1[1]) + (dir1[2] * dir1[2]) * sqrt(dir2[0] * dir2[0]) + (dir2[1] * dir2[1]) + (dir2[2] * dir2[2]);
  return nominator / denominator;
}

float dot_product2(float dir1[], float dir2[])
{
  float nominator = (dir1[0] * dir2[0]) + (dir1[1] * dir2[1]) + (dir1[2] * dir2[2]);
  float denominator = sqrt(dir1[0] * dir1[0]) + (dir1[1] * dir1[1]) + (dir1[2] * dir1[2]);
  return nominator / denominator;
}

// Method to calculate the centre of charge position for vector of pf_ids
void calculateChargeCenter(std::vector<size_t> &_nu_object_ids, const art::Event &evt, std::string _pfp_producer,
                           float chargecenter[])
{
  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);
  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);

  float totalweight = 0;

  for (auto &_i_pfp : _nu_object_ids)
  {
    // Get the associated spacepoints:
    std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(_i_pfp);

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
  }       // pfparticles
  // Normalize;
  chargecenter[0] /= totalweight;
  chargecenter[1] /= totalweight;
  chargecenter[2] /= totalweight;
}

//////////////////////////////////////////////////////
int FeatureHelper::true_thresholds_1eX(int &_true_nu_is_fiducial,
                                       std::vector<int> &_nu_daughters_pdg, std::vector<double> &_nu_daughters_E)
{
  if (_true_nu_is_fiducial)
  {
    for (unsigned int i = 0; i < _nu_daughters_pdg.size(); ++i)
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

  // First for tracks
  for (unsigned int i = 0; i < n_tracks; ++i)
  {
    double dir1[3] = {_track_dir_x[i], _track_dir_y[i], _track_dir_z[i]};
    for (unsigned int j = 0; i < n_tracks; ++i)
    {
      double dir2[3] = {_track_dir_x[j], _track_dir_y[j], _track_dir_z[j]};
      temp_cos = dot_product(dir1, dir2);
      if (temp_cos < cosine)
      {
        cosine = temp_cos;
      }
    }
    for (unsigned int j = 0; i < n_showers; ++i)
    {
      double dir2[3] = {_shower_dir_x[j], _shower_dir_y[j], _shower_dir_z[j]};
      temp_cos = dot_product(dir1, dir2);
      if (temp_cos < cosine)
      {
        cosine = temp_cos;
      }
    }
    _track_maxangle[i] = cosine;
  }

  //Now the same for showers
  for (unsigned int i = 0; i < n_showers; ++i)
  {
    double dir1[3] = {_shower_dir_x[i], _shower_dir_y[i], _shower_dir_z[i]};
    for (unsigned int j = 0; i < n_tracks; ++i)
    {
      double dir2[3] = {_track_dir_x[j], _track_dir_y[j], _track_dir_z[j]};
      temp_cos = dot_product(dir1, dir2);
      if (temp_cos < cosine)
      {
        cosine = temp_cos;
      }
    }
    for (unsigned int j = 0; i < n_showers; ++i)
    {
      double dir2[3] = {_shower_dir_x[j], _shower_dir_y[j], _shower_dir_z[j]};
      temp_cos = dot_product(dir1, dir2);
      if (temp_cos < cosine)
      {
        cosine = temp_cos;
      }
    }
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

  for (unsigned int i = 0; i < n_objects; ++i)
  {
    size_t pf_id = _nu_object_ids[i];
    recob::PFParticle const &pfparticle = pfparticle_handle->at(pf_id);

    int parent_pdg = pfparticle_handle->at(pfparticle.Parent()).PdgCode();
    const std::vector<size_t> daughter_ids = pfparticle.Daughters();
    if (parent_pdg == 11)
    {
      object_is_daughter[i] = 1;
    }
    else if (parent_pdg == 11)
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
                                                double &in_fidvol_q, double &out_fidvol_q)
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
        (sps_is_fiducial) ? in_fidvol_q += hit->Integral() : out_fidvol_q += hit->Integral();
      }
    }
  }
}

void FeatureHelper::reco_spacepoint_dqdx(double vertex[], std::vector<size_t> &_nu_object_ids, const art::Event &evt, std::string _pfp_producer,
                                         std::vector<float> &sps_dqdx_distance, std::vector<float> &sps_dqdx_integral)
{
  float chargecenter[3] = {0, 0, 0};
  calculateChargeCenter(_nu_object_ids, evt, _pfp_producer, chargecenter);

  float norm[3] = { chargecenter[0] - (float)vertex[0],
                    chargecenter[1] - (float)vertex[1],
                    chargecenter[2] - (float)vertex[2] };

  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);
  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);

  sps_dqdx_distance.clear();
  sps_dqdx_integral.clear();

  for (auto &_i_pfp : _nu_object_ids)
  {
    // Get the associated spacepoints:
    std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(_i_pfp);

    // Loop over the spacepoints and get the associated hits:
    for (auto &_sps : spcpnts)
    {
      float totalint = 0;
      float xyz[3] = { (float)_sps->XYZ()[0],(float)_sps->XYZ()[1],(float)_sps->XYZ()[2]};
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
      sps_dqdx_distance.push_back(dot_product2(norm, xyz));
      sps_dqdx_integral.push_back(totalint);
    } // spacepoints
  }   // pfparticles
}

} // namespace lee

#endif
