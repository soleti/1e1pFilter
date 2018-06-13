#ifndef ENERGYHELPER_CXX
#define ENERGYHELPER_CXX

#include "EnergyHelper.h"

namespace lee
{

EnergyHelper::EnergyHelper() {
  detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

void EnergyHelper::cluster_residuals(const art::Event &evt,
                                   std::string _pfp_producer,
                                   size_t pfp_id,
                                   double &mean_v,
                                   double &std_v)
{
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pfp_id);

  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt, _pfp_producer);

  double drift = detprop->DriftVelocity() * 1e-3;
  double fromTickToNs = 4.8 / detprop->ReadOutWindowSize() * 1e6;
  double wireSpacing = 0.3;
  std::vector<double> distances;

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    if (clusters[icl]->Plane().Plane != 2) continue;

    TVector3 start_cluster(clusters[icl]->StartWire() * wireSpacing, fromTickToNs * drift * clusters[icl]->StartTick(), 0);
    TVector3 end_cluster(clusters[icl]->EndWire() * wireSpacing, fromTickToNs * drift * clusters[icl]->EndTick(), 0);
    TVector3 line(start_cluster-end_cluster);
    
    std::vector< art::Ptr<recob::Hit> > hits = hits_per_clusters.at(clusters[icl].key());


    for (auto &hit : hits)
    {
      double w = hit->WireID().Wire * wireSpacing;
      double t = fromTickToNs * drift * hit->PeakTime();
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

void EnergyHelper::trackResiduals(const art::Event &e,
                                  std::string _pfp_producer,
                                  art::Ptr<recob::Track> candidate_track,
                                  double &mean,
                                  double &std)
{
  lar_pandora::TrackVector allPfParticleTracks;
  lar_pandora::TracksToHits trackToHitsMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, trackToHitsMap);

  std::vector<TVector3> hit_v;   // a vec of hits from coll plane
  std::vector<TVector3> track_v; // a vec of hits from coll plane

  // Collect hits
  auto iter = trackToHitsMap.find(candidate_track);
  if (iter != trackToHitsMap.end())
  {
    std::vector<art::Ptr<recob::Hit>> hits = iter->second;
    for (auto hit : hits)
    {
      if (hit->View() == 2)
      {
        TVector3 h(hit->WireID().Wire, hit->PeakTime(), 0);
        //std::cout << "emplacing hit with wire " << h.X() << ", and time " << h.Y() << std::endl;
        hit_v.emplace_back(h);
      }
    }
  }

  // Collect track points
  for (size_t i = 0; i < candidate_track->NumberTrajectoryPoints(); i++)
  {
    try
    {
      if (candidate_track->HasValidPoint(i))
      {
        TVector3 trk_pt = candidate_track->LocationAtPoint(i);
        double wire = geo->NearestWire(trk_pt, 2);
        double time = detprop->ConvertXToTicks(trk_pt.X(), geo::PlaneID(0, 0, 2));
        TVector3 p(wire, time, 0.);
        //std::cout << "emplacing track point on wire " << p.X() << ", and time " << p.Y() << std::endl;
        track_v.emplace_back(p);
      }
    }
    catch (...)
    {
      continue;
    }
  }

  ubana::TrackQuality _track_quality;
  _track_quality.SetTrackPoints(track_v);
  _track_quality.SetHitCollection(hit_v);
  std::pair<double, double> residual_mean_std = _track_quality.GetResiduals();

  if (!isnan(residual_mean_std.first) && !isnan(residual_mean_std.second)) {
    mean = residual_mean_std.first;
    std = residual_mean_std.second;
  }

}



void EnergyHelper::energy_cali(recob::PFParticle const &pfparticle,
                                    std::vector<double> &pfenergy,
                                    const art::Event &evt,
                                    std::string _pfp_producer)
{
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);

  std::vector<double> _gain;
  if (evt.isRealData()){
    _gain = _data_gain;
  }
  else{
    _gain = _mc_gain;
  }

  pfenergy.resize(3);

  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);

  std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(pfparticle.Self());
  std::vector<double> cali_corr(3, 0);
  std::vector<double> total_charge(3, 0);

  for (auto &_sps : spcpnts)
  {
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
    const double *xyz = _sps->XYZ();

    for (auto &hit : hits)
    {
      auto plane_nr = hit->View();
      if (plane_nr > 2 || plane_nr < 0)
        continue;
      total_charge[plane_nr] += hit->Integral();
      float yzcorrection = energyCalibProvider.YZdqdxCorrection(plane_nr, xyz[1], xyz[2]);
      float xcorrection = energyCalibProvider.XdqdxCorrection(plane_nr, xyz[0]);
      if (!yzcorrection)
        yzcorrection = 1.0;
      if (!xcorrection)
        xcorrection = 1.0;
      pfenergy[plane_nr] += yzcorrection * xcorrection * hit->Integral() * _gain[plane_nr] * work_function / recombination_factor / 1000; // convert MeV to GeV
    }
  }

}

void EnergyHelper::energyFromHits(recob::PFParticle const &pfparticle,
                                    std::vector<int>    &nHits,
                                    std::vector<double> &pfenergy,
                                    const art::Event &evt,
                                    std::string _pfp_producer)
{
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> cluster_pfparticle_ass(pfparticle_handle, evt,
                                                        _pfp_producer);
  std::vector<art::Ptr<recob::Cluster>> clusters =
      cluster_pfparticle_ass.at(pfparticle.Self());

      std::vector<double> _gain;
      if (evt.isRealData()){
        _gain = _data_gain;
      }
      else{
        _gain = _mc_gain;
      }

  nHits.resize(3);
  pfenergy.resize(3);

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    art::FindManyP<recob::Hit> hit_cluster_ass(cluster_handle, evt,
                                               _pfp_producer);
    std::vector<art::Ptr<recob::Hit>> hits =
        hit_cluster_ass.at(clusters[icl].key());

    for (size_t ihit = 0; ihit < hits.size(); ++ihit)
    {
      auto plane_nr = hits[ihit]->View();
      if (plane_nr > 2 || plane_nr < 0)
        continue;

      pfenergy[plane_nr] += hits[ihit]->Integral() * _gain[plane_nr] * work_function / recombination_factor /1000; // convert MeV to GeV
      nHits[plane_nr] ++;
    }
  }
}

double EnergyHelper::trackEnergy_dedx(const art::Ptr<recob::Track> &track,
                                      const art::Event &evt,
                                      std::string _pfp_producer,
                                      std::string _calo_producer)
{
  auto const &track_handle = evt.getValidHandle<std::vector<recob::Track>>(_pfp_producer);
  art::FindManyP<anab::Calorimetry> calo_track_ass(track_handle, evt, _calo_producer);

  const std::vector<art::Ptr<anab::Calorimetry>> calos = calo_track_ass.at(track->ID());

  double E = 0;


  for (size_t ical = 0; ical < calos.size(); ++ical)
  {
    if (E != 0)
      continue;
    if (!calos[ical])
      continue;

    if (!calos[ical]->PlaneID().isValid)
      continue;

    int planenum = calos[ical]->PlaneID().Plane;

    if (planenum < 0 || planenum > 2)
      continue;
    if (planenum != 2)
      continue; // Use informartion from collection plane only
    E=calos[ical]->KineticEnergy() / 1000;  // convert to GeV
  }
  return E;
}

void EnergyHelper::nHits(size_t pfp_id,
                         const art::Event &evt,
                         std::vector<int> &nHits,
                         std::string _pfp_producer)
{

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pfp_id);
  nHits.resize(3);

  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt, _pfp_producer);

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    nHits[clusters[icl]->Plane().Plane] = hits_per_clusters.at(clusters[icl].key()).size();
  }
}

void EnergyHelper::PCA(size_t pfp_id,
                       const art::Event &evt,
                       std::vector<std::vector<double>> &pca_planes,
                       std::string _pfp_producer)
{

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pfp_id);

  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt, _pfp_producer);

  double drift = detprop->DriftVelocity() * 1e-3;
  double fromTickToNs = 4.8 / detprop->ReadOutWindowSize() * 1e6;
  double wireSpacing = 0.3;

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    TPrincipal fPrincipal(2, "D");

    std::vector<art::Ptr<recob::Hit>> hits = hits_per_clusters.at(clusters[icl].key());

    for (auto &hit : hits)
    {
      double data[2];
      double w = hit->WireID().Wire * wireSpacing;
      double t = fromTickToNs * drift * hit->PeakTime();
      data[0] = w;
      data[1] = t;
      fPrincipal.AddRow(data);
    }

    fPrincipal.MakePrincipals();
    pca_planes[clusters[icl]->Plane().Plane][0] = (*fPrincipal.GetEigenValues())[0];
    pca_planes[clusters[icl]->Plane().Plane][1] = (*fPrincipal.GetEigenValues())[1];
  }
}

void EnergyHelper::dQdx_cali(size_t pfp_id,
                        const art::Event &evt,
                        std::vector<double> &dqdx,
                        std::vector<double> &dqdx_hits,
                        std::vector<double> &pitches,
                        double m_dQdxRectangleLength,
                        double m_dQdxRectangleWidth,
                        std::string _pfp_producer)
{
  double wireSpacing = 0.3;
  double tolerance = 0.001;
  std::vector<double> _gain;
  if (evt.isRealData())
    _gain = _data_gain;
  else
    _gain = _mc_gain;

  detinfo::DetectorProperties const *detprop =
      lar::providerFrom<detinfo::DetectorPropertiesService>();
  // art::ServiceHandle<geo::Geometry> geom;
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  TVector3 pfp_dir;

  //For a shower
  if (pfparticle_handle->at(pfp_id).PdgCode() == 11)
  {
    try
    {
      art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);
      auto const &shower_obj = shower_per_pfpart.at(pfp_id);
      pfp_dir.SetX(shower_obj->Direction().X());
      pfp_dir.SetY(shower_obj->Direction().Y());
      pfp_dir.SetZ(shower_obj->Direction().Z());
    }
    catch (...)
    {
      return;
    }
  }
  // For a track
  else
  {
    try
    {
    art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    auto const &track_obj = track_per_pfpart.at(pfp_id);

    pfp_dir.SetX(track_obj->StartDirection().X());
    pfp_dir.SetY(track_obj->StartDirection().Y());
    pfp_dir.SetZ(track_obj->StartDirection().Z());
    }
    catch (...)
    {
      return;
    }
  }

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt,
                                                     _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt,
                                               _pfp_producer);
  std::vector<art::Ptr<recob::Cluster>> clusters =
      clusters_per_pfpart.at(pfp_id);
  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    std::vector<art::Ptr<recob::Hit>> hits =
        hits_per_clusters.at(clusters[icl].key());

    std::vector<double> cluster_axis;
    std::vector<double> cluster_start;
    std::vector<double> cluster_end;

    double start_x = detprop->ConvertTicksToX(clusters[icl]->StartTick(), clusters[icl]->Plane());
    double end_x = detprop->ConvertTicksToX(clusters[icl]->EndTick(), clusters[icl]->Plane());
    // std::cout << "start " << start_x << ", end: " << end_x << std::endl;
    if (pfp_dir.Z() >= 0)
    {
      std::reverse(hits.begin(), hits.end());
      cluster_axis = {cos(clusters[icl]->StartAngle()),
                      sin(clusters[icl]->StartAngle())};

      cluster_start = {clusters[icl]->StartWire() * wireSpacing - tolerance*cos(clusters[icl]->StartAngle()),
                       start_x - tolerance*sin(clusters[icl]->StartAngle())};
      cluster_end = {clusters[icl]->EndWire() * wireSpacing, end_x};
    }
    else
    {
      cluster_axis = {-1.*cos(clusters[icl]->StartAngle()),
                      -1.*sin(clusters[icl]->StartAngle())};
      cluster_start = {clusters[icl]->EndWire() * wireSpacing + tolerance*cos(clusters[icl]->StartAngle()),
                       end_x + tolerance*sin(clusters[icl]->StartAngle())};
      cluster_end = {clusters[icl]->StartWire() * wireSpacing, start_x};
    }

    double cluster_length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) +
                                 pow(cluster_end[1] - cluster_start[1], 2));
    if (cluster_length <= 0)
      continue;

    double pitch =
        geoHelper.getPitch(pfp_dir, clusters[icl]->Plane().Plane);

    // Build rectangle 4 x 1 cm around the cluster axis
    std::vector<std::vector<double>> points;
    geoHelper.buildRectangle(m_dQdxRectangleLength, m_dQdxRectangleWidth,
                             cluster_start, cluster_axis, points);

    std::vector<double> dqdxs;
    art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);
    
    std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(pfp_id);

    for (auto &_sps : spcpnts)
    {
      std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
      const double *xyz = _sps->XYZ();
      for (auto &hit : hits)
      {
        auto plane_nr = hit->View();
        if (plane_nr > 2 || plane_nr < 0)
          continue;



        std::vector<double> hit_pos = {hit->WireID().Wire * wireSpacing,
                                       detprop->ConvertTicksToX(hit->PeakTime(), clusters[icl]->Plane())};

        bool is_within = geoHelper.isInside(hit_pos, points);

        if (is_within)
        {
          float yzcorrection = energyCalibProvider.YZdqdxCorrection(plane_nr, xyz[1], xyz[2]);
          float xcorrection = energyCalibProvider.XdqdxCorrection(plane_nr, xyz[0]);
          double q = yzcorrection * xcorrection * hit->Integral() * _gain[clusters[icl]->Plane().Plane];
          if (!yzcorrection)
            yzcorrection = 1.0;
          if (!xcorrection)
            xcorrection = 1.0;
          dqdxs.push_back(q / pitch);
          if (plane_nr == 2)
          {
            dqdx_hits.push_back(q / pitch);
            // dqdx_wires.push_back(hit->WireID().Wire);
          }
        }
      }

    }
    
    pitches[clusters[icl]->Plane().Plane] = pitch;

    // Get the median
    if (dqdxs.size() > 0)
    {
      std::nth_element(dqdxs.begin(), dqdxs.begin() + dqdxs.size()/2, dqdxs.end());
      dqdx[clusters[icl]->Plane().Plane] = dqdxs[dqdxs.size()/2];
    }
  }
}

void EnergyHelper::dQdx(size_t pfp_id,
                        const art::Event &evt,
                        std::vector<double> &dqdx,
                        std::vector<double> &dqdx_hits,
                        std::vector<double> &pitches,
                        double m_dQdxRectangleLength,
                        double m_dQdxRectangleWidth,
                        std::string _pfp_producer)
{
  double wireSpacing = 0.3;
  double tolerance = 0.001;
  std::vector<double> _gain;
  if (evt.isRealData())
    _gain = _data_gain;
  else
    _gain = _mc_gain;

  detinfo::DetectorProperties const *detprop =
      lar::providerFrom<detinfo::DetectorPropertiesService>();
  // art::ServiceHandle<geo::Geometry> geom;

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  TVector3 pfp_dir;

  //For a shower
  if (pfparticle_handle->at(pfp_id).PdgCode() == 11)
  {
    try
    {
      art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);
      auto const &shower_obj = shower_per_pfpart.at(pfp_id);
      pfp_dir.SetX(shower_obj->Direction().X());
      pfp_dir.SetY(shower_obj->Direction().Y());
      pfp_dir.SetZ(shower_obj->Direction().Z());
    }
    catch (...)
    {
      return;
    }
  }
  // For a track
  else
  {
    try
    {
    art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    auto const &track_obj = track_per_pfpart.at(pfp_id);

    pfp_dir.SetX(track_obj->StartDirection().X());
    pfp_dir.SetY(track_obj->StartDirection().Y());
    pfp_dir.SetZ(track_obj->StartDirection().Z());
    }
    catch (...)
    {
      return;
    }
  }

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt,
                                                     _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt,
                                               _pfp_producer);
  std::vector<art::Ptr<recob::Cluster>> clusters =
      clusters_per_pfpart.at(pfp_id);
  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    std::vector<art::Ptr<recob::Hit>> hits =
        hits_per_clusters.at(clusters[icl].key());

    std::vector<double> cluster_axis;
    std::vector<double> cluster_start;
    std::vector<double> cluster_end;

    double start_x = detprop->ConvertTicksToX(clusters[icl]->StartTick(), clusters[icl]->Plane());
    double end_x = detprop->ConvertTicksToX(clusters[icl]->EndTick(), clusters[icl]->Plane());
    // std::cout << "start " << start_x << ", end: " << end_x << std::endl;
    if (pfp_dir.Z() >= 0)
    {
      std::reverse(hits.begin(), hits.end());
      cluster_axis = {cos(clusters[icl]->StartAngle()),
                      sin(clusters[icl]->StartAngle())};

      cluster_start = {clusters[icl]->StartWire() * wireSpacing - tolerance*cos(clusters[icl]->StartAngle()),
                       start_x - tolerance*sin(clusters[icl]->StartAngle())};
      cluster_end = {clusters[icl]->EndWire() * wireSpacing, end_x};
    }
    else
    {
      cluster_axis = {-1.*cos(clusters[icl]->StartAngle()),
                      -1.*sin(clusters[icl]->StartAngle())};
      cluster_start = {clusters[icl]->EndWire() * wireSpacing + tolerance*cos(clusters[icl]->StartAngle()),
                       end_x + tolerance*sin(clusters[icl]->StartAngle())};
      cluster_end = {clusters[icl]->StartWire() * wireSpacing, start_x};
    }

    double cluster_length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) +
                                 pow(cluster_end[1] - cluster_start[1], 2));
    if (cluster_length <= 0)
      continue;

    double pitch =
        geoHelper.getPitch(pfp_dir, clusters[icl]->Plane().Plane);

    // Build rectangle 4 x 1 cm around the cluster axis
    std::vector<std::vector<double>> points;
    geoHelper.buildRectangle(m_dQdxRectangleLength, m_dQdxRectangleWidth,
                             cluster_start, cluster_axis, points);

    std::vector<double> dqdxs;

    for (auto &hit : hits)
    {
      std::vector<double> hit_pos = {hit->WireID().Wire * wireSpacing,
                                     detprop->ConvertTicksToX(hit->PeakTime(), clusters[icl]->Plane())};

      bool is_within = geoHelper.isInside(hit_pos, points);

      if (is_within)
      {
        double q = hit->Integral() * _gain[clusters[icl]->Plane().Plane];
        dqdxs.push_back(q / pitch);
        if (clusters[icl]->Plane().Plane == 2)
        {
          dqdx_hits.push_back(q / pitch);
          // dqdx_wires.push_back(hit->WireID().Wire);
        }
      }
    }

    pitches[clusters[icl]->Plane().Plane] = pitch;

    // Get the median
    if (dqdxs.size() > 0)
    {
      std::nth_element(dqdxs.begin(), dqdxs.begin() + dqdxs.size()/2, dqdxs.end());
      dqdx[clusters[icl]->Plane().Plane] = dqdxs[dqdxs.size()/2];
    }
  }
}

void EnergyHelper::dEdxFromdQdx(std::vector<double> &dedx,
                                std::vector<double> &dqdx)
{
  for (size_t i = 0; i < dqdx.size(); i++)
  {
    if (dqdx[i] > 0)
      dedx[i] = dqdx[i] * (work_function) / recombination_factor;
  }
}

} // namespace lee

#endif
