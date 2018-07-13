#ifndef ENERGYHELPER_CXX
#define ENERGYHELPER_CXX

#include "EnergyHelper.h"

namespace lee
{

EnergyHelper::EnergyHelper(bool is_data=false) {
    _is_data = is_data;
}


void EnergyHelper::track_dQdx(std::vector<art::Ptr<anab::Calorimetry>> *calos,
                              std::vector<double> &dqdx,
                              std::vector<double> &dedx) 
{

  int planenum = -1;
  for (auto c : *calos) {
    if (!c)
      continue; 

    planenum = c->PlaneID().Plane;
    std::vector<double> dqdxs = c->dQdx();
    std::vector<double> dedxs = c->dEdx();

    if (dqdxs.size() <= 1) {
      continue;
    }
    
    size_t n = dqdxs.size() / 2;
    std::nth_element(dqdxs.begin(), dqdxs.begin() + n, dqdxs.end());
    dqdx[planenum] = dqdxs[n];
    if (dedxs.size() == dqdxs.size()) {
      dedx[planenum] = dedxs[n];
    }
  }

}

void EnergyHelper::cluster_residuals(std::vector<art::Ptr<recob::Cluster>> *clusters,
                                     art::FindManyP<recob::Hit> *hits_per_cluster,
                                     double &mean_v,
                                     double &std_v)
{

  std::vector<double> distances;

  for (auto _cl: *clusters) {

    if (_cl->Plane().Plane != 2) continue;

    TVector3 start_cluster(_cl->StartWire() * _wire_spacing, _from_tick_to_ns * _drift * _cl->StartTick(), 0);
    TVector3 end_cluster(_cl->EndWire() * _wire_spacing, _from_tick_to_ns * _drift * _cl->EndTick(), 0);
    TVector3 line(start_cluster-end_cluster);
    
    std::vector< art::Ptr<recob::Hit> > hits = hits_per_cluster->at(_cl.key());

    for (auto &hit : hits)
    {
      double w = hit->WireID().Wire * _wire_spacing;
      double t = _from_tick_to_ns * _drift * hit->PeakTime();
      TVector3 hit_v(w, t, 0);
      TVector3 num(line.Cross(start_cluster - hit_v));
      double side = (hit_v[0] - start_cluster[0])*(end_cluster[1] - start_cluster[1]) - (hit_v[1] - start_cluster[1])*(end_cluster[0] - start_cluster[0]);
      int sign_side = (side > 0) - (side < 0);
      distances.push_back(sign_side * num.Mag()/line.Mag());
    }

  }

  double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
  double mean = sum / distances.size();

  std::vector<double> diff(distances.size());
  std::transform(distances.begin(), distances.end(), diff.begin(), [mean](double x) { return x - mean; });

  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdv = std::sqrt(sq_sum / distances.size());

  if (!isnan(mean) && !isnan(stdv)) {
    mean_v = mean;
    std_v = stdv;
  }

}

void EnergyHelper::energy_from_hits(std::vector<art::Ptr<recob::Cluster>> *clusters,
                                    art::FindManyP<recob::Hit> *hits_per_cluster,
                                    std::vector<int>    &nHits,
                                    std::vector<double> &pfenergy)
{

  std::vector<double> gain;

  if (_is_data) {
      gain = _data_gain;
  } else {
      gain = _mc_gain;
  }

  nHits.resize(3);
  pfenergy.resize(3);

  for (auto _cl: *clusters) {
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster->at(_cl.key());

    for (auto &hit : hits)
    {
      auto plane_nr = hit->View();
      if (plane_nr > 2 || plane_nr < 0)
        continue;

      pfenergy[plane_nr] += hit->Integral() * gain[plane_nr] * _work_function / _recombination_factor / 1000; // convert MeV to GeV
      nHits[plane_nr]++;
    }
  }
}

void EnergyHelper::PCA(std::vector<art::Ptr<recob::Cluster>> *clusters,
                       art::FindManyP<recob::Hit> *hits_per_cluster,
                       std::vector<std::vector<double>> &pca_planes)
{

  for (auto _cl: *clusters)
  {
    TPrincipal fPrincipal(2, "D");

    std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster->at(_cl.key());
    for (auto &hit : hits)
    {
      
      double data[2];
      double w = hit->WireID().Wire * _wire_spacing;
      double t = _from_tick_to_ns * _drift * hit->PeakTime();
      data[0] = w;
      data[1] = t;
      fPrincipal.AddRow(data);
    }

    fPrincipal.MakePrincipals();
    pca_planes[0][_cl->Plane().Plane] = (*fPrincipal.GetEigenValues())[0];
    pca_planes[1][_cl->Plane().Plane] = (*fPrincipal.GetEigenValues())[1];
  }

}

void EnergyHelper::get_cali(
    std::vector<art::Ptr<recob::SpacePoint>> *spcpnts,
    art::FindManyP<recob::Hit> *hits_per_spcpnts,
    std::vector<double> &cali_corr)
{
  cali_corr.resize(3);
  std::vector<float> total_charge(3, 0);

  for (auto _sps : *spcpnts)
  {
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts->at(_sps.key());
    const double *xyz = _sps->XYZ();

    for (auto &hit : hits)
    {
      auto plane_nr = hit->View();
      if (plane_nr > 2 || plane_nr < 0)
        continue;
      total_charge[plane_nr] += hit->Integral();
      float yzcorrection = _energy_calib_provider.YZdqdxCorrection(plane_nr, xyz[1], xyz[2]);
      float xcorrection = _energy_calib_provider.XdqdxCorrection(plane_nr, xyz[0]);
      if (!yzcorrection)
        yzcorrection = 1.0;
      if (!xcorrection)
        xcorrection = 1.0;
      cali_corr[plane_nr] += yzcorrection * xcorrection * hit->Integral();
    }
  }

  for (unsigned short i = 0; i < 3; ++i)
  {
    if (total_charge[i] > 0)
    {
      cali_corr[i] /= total_charge[i];
    }
    else
    {
      cali_corr[i] = 1;
    }
  }
}


double EnergyHelper::PID(const std::vector<anab::ParticleID> *pids,
                         int trackID,
                         std::string AlgName,
                         anab::kVariableType VariableType,
                         int pdgCode)
{
    anab::ParticleID selected_pid;
    bool there_is_pid = false;

    for (auto& pid: *pids) {
      if ((int)pid.PlaneID().Plane == trackID) {
          selected_pid = pid;
          there_is_pid = true;
          break;
      }
    }

    if (!there_is_pid) {
        return std::numeric_limits<double>::lowest();
    }
    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = selected_pid.ParticleIDAlgScores();

    for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
    {
      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
      int planeid = AlgScore.fPlaneID.Plane;

      if (planeid < 0 || planeid > 2)
      {
        std::cout << "[EnergyHelper::PID] No information for planeid " << planeid << std::endl;
        continue;
      }

      if (AlgScore.fAlgName == AlgName)
      {
        if (anab::kVariableType(AlgScore.fVariableType) == VariableType)
        {
          if (AlgScore.fAssumedPdg == pdgCode) {
              double alg_value = AlgScore.fValue;
              return alg_value;
          }
        }
      }
    }
    return std::numeric_limits<double>::lowest();
}

void EnergyHelper::dQdx(size_t pfp_id,
                        const art::Event &evt,
                        std::vector<double> &dqdx,
                        std::vector<double> &dqdx_cali,
                        std::vector<double> &dqdx_hits,
                        std::vector<double> &pitches,
                        double m_dQdxRectangleLength,
                        double m_dQdxRectangleWidth,
                        std::string _pfp_producer)
{

  double tolerance = 0.001;

  std::vector<double> _gain;
  if (evt.isRealData())
    _gain = _data_gain;
  else
    _gain = _mc_gain;


  // art::ServiceHandle<geo::Geometry> geom;

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  TVector3 pfp_dir;

  // Field needed for calibration factor
  double x_start, y_start, z_start;
  double x_middle, y_middle, z_middle;
  double x_end, y_end, z_end;
  float start_corr, middle_corr, end_corr;

  //For a shower
  if (pfparticle_handle->at(pfp_id).PdgCode() == 11)
  {

    art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    auto const &shower_obj = shower_per_pfpart.at(pfp_id);

    try
    {
      pfp_dir.SetX(shower_obj->Direction().X());
      pfp_dir.SetY(shower_obj->Direction().Y());
      pfp_dir.SetZ(shower_obj->Direction().Z());

      x_start = shower_obj->ShowerStart().X();
      y_start = shower_obj->ShowerStart().Y();
      z_start = shower_obj->ShowerStart().Z();
    }
    catch (...)
    {
      return;
    }
  }
  // For a track
  else
  {

    art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    auto const &track_obj = track_per_pfpart.at(pfp_id);

    pfp_dir.SetX(track_obj->StartDirection().X());
    pfp_dir.SetY(track_obj->StartDirection().Y());
    pfp_dir.SetZ(track_obj->StartDirection().Z());

    x_start = track_obj->Start().X();
    y_start = track_obj->Start().Y();
    z_start = track_obj->Start().Z();
  }
  pfp_dir.SetMag(2.); //Go 2cm along the direction of the object.
  x_middle = x_start + pfp_dir.X();
  y_middle = y_start + pfp_dir.Y();
  z_middle = z_start + pfp_dir.Z();
  x_end = x_middle + pfp_dir.X();
  y_end = y_middle + pfp_dir.Y();
  z_end = z_middle + pfp_dir.Z();
  pfp_dir.SetMag(1.); //Normalise again for safety (not needed).

  for (int plane_nr = 0; plane_nr < 3; ++plane_nr)
  {
    float yzcorrection_start = _energy_calib_provider.YZdqdxCorrection(plane_nr, y_start, z_start);
    float xcorrection_start = _energy_calib_provider.XdqdxCorrection(plane_nr, x_start);
    if (!yzcorrection_start)
      yzcorrection_start = 1.0;
    if (!xcorrection_start)
      xcorrection_start = 1.0;
    start_corr = yzcorrection_start * xcorrection_start;

    float yzcorrection_middle = _energy_calib_provider.YZdqdxCorrection(plane_nr, y_middle, z_middle);
    float xcorrection_middle = _energy_calib_provider.XdqdxCorrection(plane_nr, x_middle);
    if (!yzcorrection_middle)
      yzcorrection_middle = 1.0;
    if (!xcorrection_middle)
      xcorrection_middle = 1.0;
    middle_corr = yzcorrection_middle * xcorrection_middle;

    float yzcorrection_end = _energy_calib_provider.YZdqdxCorrection(plane_nr, y_end, z_end);
    float xcorrection_end = _energy_calib_provider.XdqdxCorrection(plane_nr, x_end);
    if (!yzcorrection_end)
      yzcorrection_end = 1.0;
    if (!xcorrection_end)
      xcorrection_end = 1.0;
    end_corr = yzcorrection_end * xcorrection_end;
    //std::cout << "[EnergyHelper] dqdx_cali " << start_corr << middle_corr << end_corr << std::endl;
    dqdx_cali[plane_nr] = (start_corr + middle_corr + end_corr) / 3;
  }

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_cluster(cluster_handle, evt, _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pfp_id);

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    std::vector<art::Ptr<recob::Hit>> hits =
        hits_per_cluster.at(clusters[icl].key());

    std::vector<double> cluster_axis;
    std::vector<double> cluster_start;
    std::vector<double> cluster_end;
    double start_x = _detprop->ConvertTicksToX(clusters[icl]->StartTick(), clusters[icl]->Plane());
    double end_x = _detprop->ConvertTicksToX(clusters[icl]->EndTick(), clusters[icl]->Plane());
    // std::cout << "start " << start_x << ", end: " << end_x << std::endl;
    if (pfp_dir.Z() >= 0)
    {
      std::reverse(hits.begin(), hits.end());
      cluster_axis = {cos(clusters[icl]->StartAngle()),
                      sin(clusters[icl]->StartAngle())};

      cluster_start = {clusters[icl]->StartWire() * _wire_spacing - tolerance * cos(clusters[icl]->StartAngle()),
                       start_x - tolerance * sin(clusters[icl]->StartAngle())};
      cluster_end = {clusters[icl]->EndWire() * _wire_spacing, end_x};
    }
    else
    {
      cluster_axis = {-1. * cos(clusters[icl]->StartAngle()),
                      -1. * sin(clusters[icl]->StartAngle())};
      cluster_start = {clusters[icl]->EndWire() * _wire_spacing + tolerance * cos(clusters[icl]->StartAngle()),
                       end_x + tolerance * sin(clusters[icl]->StartAngle())};
      cluster_end = {clusters[icl]->StartWire() * _wire_spacing, start_x};
    }

    double cluster_length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) +
                                 pow(cluster_end[1] - cluster_start[1], 2));
    if (cluster_length <= 0)
      continue;

    double pitch = geo_helper.getPitch(pfp_dir, clusters[icl]->Plane().Plane);

    // Build rectangle 4 x 1 cm around the cluster axis
    std::vector<std::vector<double>> points;
    geo_helper.buildRectangle(m_dQdxRectangleLength, m_dQdxRectangleWidth,
                             cluster_start, cluster_axis, points);

    std::vector<double> dqdxs;

    for (auto &hit : hits)
    {
      std::vector<double> hit_pos = {hit->WireID().Wire * _wire_spacing,
                                     _detprop->ConvertTicksToX(hit->PeakTime(), clusters[icl]->Plane())};

      bool is_within = geo_helper.isInside(hit_pos, points);

      if (is_within)
      {
        double q = hit->Integral() * _gain[clusters[icl]->Plane().Plane];
        dqdxs.push_back(q / pitch);
        dqdx_hits.push_back(q / pitch);
        pitches.push_back(pitch);
      }
    }

    // Get the median
    if (dqdxs.size() > 0)
    {
      std::nth_element(dqdxs.begin(), dqdxs.begin() + dqdxs.size() / 2, dqdxs.end());
      dqdx[clusters[icl]->Plane().Plane] = dqdxs[dqdxs.size() / 2];
    }
  }
}

void EnergyHelper::dEdx_from_dQdx(std::vector<double> &dedx,
                                std::vector<double> &dqdx)
{
  for (size_t i = 0; i < dqdx.size(); i++)
  {
    if (dqdx[i] > 0)
      dedx[i] = dqdx[i] * (_work_function) / _recombination_factor;
  }
}

} // namespace lee

#endif
