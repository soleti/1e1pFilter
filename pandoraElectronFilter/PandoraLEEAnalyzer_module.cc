////////////////////////////////////////////////////////////////////////
// Class:       PandoraLEEAnalyzer
// Module Type: analyzer
// File:        PandoraLEEAnalyzer_module.cc
//
// Generated at Thu Jun 23 00:24:52 2016 by Lorena Escudero Sanchez using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "PandoraLEEAnalyzer.h"


lee::PandoraLEEAnalyzer::PandoraLEEAnalyzer(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)  // ,
// More initializers here.
{

  //create output tree
  art::ServiceHandle<art::TFileService> tfs;
  // myTFile = new TFile("PandoraLEEAnalyzerOutput.root", "RECREATE");
  myTTree = tfs->make<TTree>("pandoratree", "PandoraAnalysis Tree");

  myPOTTTree = tfs->make<TTree>("pot", "POT Tree");

  myTTree->Branch("category",  &_category, "category/i");
  myTTree->Branch("E",  &_energy, "E/d");
  myTTree->Branch("n_tracks",  &_n_tracks, "n_tracks/i");
  myTTree->Branch("n_showers",  &_n_showers, "n_showers/i");
  myTTree->Branch("vx",  &_vx, "vx/d");
  myTTree->Branch("vy",  &_vy, "vy/d");
  myTTree->Branch("vz",  &_vz, "vz/d");

  myTTree->Branch("true_vx",  &_true_vx, "true_vx/d");
  myTTree->Branch("true_vy",  &_true_vy, "true_vy/d");
  myTTree->Branch("true_vz",  &_true_vz, "true_vz/d");

  myTTree->Branch("true_vx_sce",  &_true_vx_sce, "true_vx_sce/d");
  myTTree->Branch("true_vy_sce",  &_true_vy_sce, "true_vy_sce/d");
  myTTree->Branch("true_vz_sce",  &_true_vz_sce, "true_vz_sce/d");

  myTTree->Branch("nu_E",  &_nu_energy, "nu_E/d");
  myTTree->Branch("passed",  &_event_passed, "passed/I");
  myTTree->Branch("n_candidates", &_n_candidates, "n_candidates/i");
  myTTree->Branch("n_true_nu", &_n_true_nu, "n_true_nu/i");
  myTTree->Branch("distance", &_distance, "distance/d");
  myTTree->Branch("true_nu_is_fiducial", &_true_nu_is_fiducial, "true_nu_is_fiducial/I");

  myTTree->Branch("n_matched", &_n_matched, "n_matched/i");
  myTTree->Branch("nu_matched_tracks", &_nu_matched_tracks, "nu_matched_tracks/i");
  myTTree->Branch("nu_matched_showers", &_nu_matched_showers, "nu_matched_showers/i");

  myTTree->Branch("nu_daughters_pdg",  "std::vector< int >", &_nu_daughters_pdg);
  myTTree->Branch("nu_daughters_E",  "std::vector< double >", &_nu_daughters_E);

  myTTree->Branch("nu_daughters_vx",  "std::vector< double >", &_nu_daughters_vx);
  myTTree->Branch("nu_daughters_vy",  "std::vector< double >", &_nu_daughters_vy);
  myTTree->Branch("nu_daughters_vz",  "std::vector< double >", &_nu_daughters_vz);

  myTTree->Branch("nu_daughters_endx",  "std::vector< double >", &_nu_daughters_endx);
  myTTree->Branch("nu_daughters_endy",  "std::vector< double >", &_nu_daughters_endy);
  myTTree->Branch("nu_daughters_endz",  "std::vector< double >", &_nu_daughters_endz);

  myTTree->Branch("nu_daughters_px",  "std::vector< double >", &_nu_daughters_px);
  myTTree->Branch("nu_daughters_py",  "std::vector< double >", &_nu_daughters_py);
  myTTree->Branch("nu_daughters_pz",  "std::vector< double >", &_nu_daughters_pz);

  myTTree->Branch("event", &_event, "event/i");
  myTTree->Branch("run", &_run, "run/i");
  myTTree->Branch("subrun", &_subrun, "subrun/i");

  myTTree->Branch("bnbweight", &_bnbweight, "bnbweight/d");


  myTTree->Branch("flash_passed", &_flash_passed, "flash_passed/i");
  myTTree->Branch("track_passed", &_track_passed, "track_passed/i");
  myTTree->Branch("shower_passed", &_shower_passed, "shower_passed/i");


  myTTree->Branch("shower_dir_x",  "std::vector< double >", &_shower_dir_x);
  myTTree->Branch("shower_dir_y",  "std::vector< double >", &_shower_dir_y);
  myTTree->Branch("shower_dir_z",  "std::vector< double >", &_shower_dir_z);

  myTTree->Branch("shower_start_x",  "std::vector< double >", &_shower_start_x);
  myTTree->Branch("shower_start_y",  "std::vector< double >", &_shower_start_y);
  myTTree->Branch("shower_start_z",  "std::vector< double >", &_shower_start_z);

  myTTree->Branch("shower_theta",  "std::vector< double >", &_shower_theta);
  myTTree->Branch("shower_phi",  "std::vector< double >", &_shower_phi);

  myTTree->Branch("shower_energy",  "std::vector< double >", &_shower_energy);
  myTTree->Branch("track_energy",  "std::vector< double >", &_track_energy);

  myTTree->Branch("track_dir_x",  "std::vector< double >", &_track_dir_x);
  myTTree->Branch("track_dir_y",  "std::vector< double >", &_track_dir_y);
  myTTree->Branch("track_dir_z",  "std::vector< double >", &_track_dir_z);

  myTTree->Branch("track_start_x",  "std::vector< double >", &_track_start_x);
  myTTree->Branch("track_start_y",  "std::vector< double >", &_track_start_y);
  myTTree->Branch("track_start_z",  "std::vector< double >", &_track_start_z);

  myTTree->Branch("track_is_fiducial",  "std::vector< int >", &_track_is_fiducial);
  myTTree->Branch("shower_is_fiducial",  "std::vector< int >", &_shower_is_fiducial);


  myTTree->Branch("track_end_x",  "std::vector< double >", &_track_end_x);
  myTTree->Branch("track_end_y",  "std::vector< double >", &_track_end_y);
  myTTree->Branch("track_end_z",  "std::vector< double >", &_track_end_z);

  myTTree->Branch("track_theta",  "std::vector< double >", &_track_theta);
  myTTree->Branch("track_phi",  "std::vector< double >", &_track_phi);

  myTTree->Branch("track_len",  "std::vector< double >", &_track_length);
  myTTree->Branch("track_id",  "std::vector< double >", &_track_id);

  myTTree->Branch("nu_pdg",  &_nu_pdg, "nu_pdg/i");

  myPOTTTree->Branch("run", &_run_sr, "run/i");
  myPOTTTree->Branch("subrun", &_subrun_sr, "subrun/i");
  myPOTTTree->Branch("pot", &_pot, "pot/d");

  myTTree->Branch("predict_p","std::vector< double >",&_predict_p);
  myTTree->Branch("predict_mu","std::vector< double >",&_predict_mu);
  myTTree->Branch("predict_pi","std::vector< double >",&_predict_pi);
  myTTree->Branch("predict_em","std::vector< double >",&_predict_em);
  myTTree->Branch("predict_cos","std::vector< double >",&_predict_cos);

  myTTree->Branch("interaction_type",&_interaction_type, "interaction_type/i");

  myTTree->Branch("shower_dQdx","std::vector< std::vector< double > >",&_shower_dQdx);
  myTTree->Branch("shower_dEdx","std::vector< std::vector< double > >",&_shower_dEdx);

  myTTree->Branch("shower_open_angle","std::vector< double >",&_shower_open_angle);

  myTTree->Branch("matched_tracks","std::vector< int > _matched_tracks",&_matched_tracks);
  myTTree->Branch("matched_showers","std::vector< int > _matched_tracks",&_matched_showers);



  this->reconfigure(pset);

}

lee::PandoraLEEAnalyzer::~PandoraLEEAnalyzer()
{
  std::cout << "[PandoraLEE] " << "End!" << std::endl;
}



art::Ptr<recob::Shower> lee::PandoraLEEAnalyzer::get_most_energetic_shower(std::vector< art::Ptr<recob::Shower> > &showers) {
  art::Ptr<recob::Shower> most_energetic_shower;

  double max_energy = std::numeric_limits<double>::lowest();
  for (auto const& shower : showers) {
    if (shower->Energy()[shower->best_plane()] > max_energy) {
      most_energetic_shower = shower;
      max_energy = shower->Energy()[shower->best_plane()];
    }
  }
  return most_energetic_shower;
}

art::Ptr<recob::Track> lee::PandoraLEEAnalyzer::get_longest_track(std::vector< art::Ptr<recob::Track> > &tracks) {
  art::Ptr<recob::Track> longest_track;

  double max_length = std::numeric_limits<double>::lowest();
  for (auto const& track : tracks) {
    try {
      if (track->Length() > max_length) {
        longest_track = track;
        max_length = track->Length();
      }
    } catch (...) {
      std::cout << "[PandoraLEE] " << "Error getting longest track " << track << std::endl;
    }
  }
  return longest_track;
}

void lee::PandoraLEEAnalyzer::get_daughter_tracks( std::vector < size_t > pf_ids, const art::Event & evt, std::vector< art::Ptr<recob::Track> > &tracks) {

  auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( _pfp_producer );

  art::FindOneP< recob::Track > track_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  for (auto const& pf_id : pf_ids) {
    try {
      auto const& track_obj = track_per_pfpart.at(pf_id);
      tracks.push_back(track_obj);
    } catch (...) {
      std::cout << "[PandoraLEE] " << "Error getting the track" << std::endl;
    }
  }
}

void lee::PandoraLEEAnalyzer::get_daughter_showers(std::vector < size_t > pf_ids, const art::Event & evt, std::vector< art::Ptr<recob::Shower> > &showers) {

  auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( _pfp_producer );

  art::FindOneP< recob::Shower > shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  for (auto const& pf_id : pf_ids) {
    try {

      auto const& shower_obj = shower_per_pfpart.at(pf_id);
      showers.push_back(shower_obj);
    } catch (...) {
      std::cout << "[PandoraLEE] " << "Error getting the shower" << std::endl;
    }
  }

}

void lee::PandoraLEEAnalyzer::dEdxFromdQdx(std::vector< double > &dedx, std::vector< double > &dqdx) {

  double work_function = 23;
  double recombination_factor = 0.62;

  for (size_t i = 0; i < dqdx.size(); i++) {
    if (dqdx[i] > 0) dedx[i] = dqdx[i]*(work_function/1e6)/recombination_factor;
    std::cout << "[dEdx] " << i << " " << dedx[i] << std::endl;
  }

}

double lee::PandoraLEEAnalyzer::getPitch(const TVector3& direction, const int& pl) const
{
  // prepare a direction vector for the plane
  TVector3 wireDir = {0., 0., 0.};
  // the direction of the plane is the vector uniting two consecutive wires
  // such that this vector is perpendicular to both wires
  // basically this is the vector perpendicular to the wire length direction,
  // and still in the wire-plane direction
  if (pl == 0)
  wireDir = {0., -sqrt(3) / 2., 1 / 2.};
  else if (pl == 1)
  wireDir = {0.,  sqrt(3) / 2., 1 / 2.};
  else if (pl == 2)
  wireDir = {0., 0., 1.};

  // cosine between shower direction and plane direction gives the factor
  // by which to divide 0.3, the minimum wire-spacing
  double minWireSpacing = 0.3;
  double cos = wireDir.Dot(direction);
  if (cos < 0) cos *= -1;
  cos /= (wireDir.Mag() * direction.Mag());
  // if cosine is 0 the direction is perpendicular and the wire-spacing is infinite
  if (cos == 0)
  return std::numeric_limits<double>::max();
  double pitch = minWireSpacing / cos;
  return pitch;
}

bool lee::PandoraLEEAnalyzer::withinRectangle(std::vector< std::vector < double > > &points, std::vector< double > &point) {
  double m1_top = (points[1][1]-points[0][1])/(points[1][0]-points[0][0]);
  double m2_top = (points[1][1]-points[3][1])/(points[1][0]-points[3][0]);
  double m1_bottom = (points[3][1]-points[2][1])/(points[3][0]-points[2][0]);
  double m2_bottom = (points[0][1]-points[2][1])/(points[0][0]-points[2][0]);

  bool is_below_m1_top = point[1] < m1_top*(point[0]-points[1][0]) + points[1][1];
  bool is_below_m2_top = point[1] < m2_top*(point[0]-points[1][0]) + points[1][1];
  bool is_above_m1_bottom = point[1] > m1_bottom*(point[0]-points[2][0]) + points[2][1];
  bool is_above_m2_bottom = point[1] > m2_bottom*(point[0]-points[2][0]) + points[2][1];

  return is_below_m1_top and is_below_m2_top and is_above_m2_bottom and is_above_m1_bottom;
}

void lee::PandoraLEEAnalyzer::buildRectangle(double length, double width, std::vector< double > &start, std::vector< double > &axis, std::vector< std::vector < double > > &points) {
  double perp_axis[2] = {-axis[1],axis[0]};

  std::vector< double > p1 = {start[0]+perp_axis[0]*width/2, start[1]+perp_axis[1]*width/2};
  std::vector< double > p2 = {p1[0]+axis[0]*length,p1[1]+axis[1]*length};
  std::vector< double > p3 = {start[0]-perp_axis[0]*width/2,start[1]-perp_axis[1]*width/2};
  std::vector< double > p4 = {p3[0]+axis[0]*length,p3[1]+axis[1]*length};

  points.insert(points.end(), { p1, p2, p3, p4});
}

void lee::PandoraLEEAnalyzer::dQdx(size_t pfp_id, const art::Event & evt, std::vector< double > &dqdx) {



  detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;

  auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( _pfp_producer );
  auto const& cluster_handle = evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindOneP< recob::Shower > shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  auto const& shower_obj = shower_per_pfpart.at(pfp_id);

  TVector3 shower_dir(shower_obj->Direction().X(), shower_obj->Direction().Y(), shower_obj->Direction().Z());

  art::FindManyP<recob::Cluster > clusters_per_pfpart ( pfparticle_handle, evt, _pfp_producer );
  art::FindManyP<recob::Hit > hits_per_clusters ( cluster_handle, evt, _pfp_producer );

  std::vector<art::Ptr < recob::Cluster > > clusters = clusters_per_pfpart.at(pfp_id);

  double drift = detprop->DriftVelocity()*1e-3;
  std::cout << drift << std::endl;
  std::cout << "[dQdx] Clusters size " << clusters.size() << std::endl;

  for (size_t icl = 0; icl < clusters.size(); icl++) {
    std::vector<art::Ptr<recob::Hit> > hits = hits_per_clusters.at(clusters[icl].key());

    std::cout << "[dQdx] " << "Cluster " << icl << std::endl;

    std::cout << "[dQdx] " << "Wire coordinate " << clusters[icl]->StartWire() << std::endl;
    std::cout << "[dQdx] " << "Tick coordinate " << clusters[icl]->StartTick() << std::endl;

    // TODO Use variable from detector properties!
    // To get the time in ns -> 4.8 ms / 9600 ticks * 1e6 = 500
    // 0.3 wire spacing

    double fromTickToNs = 4.8 / detprop->ReadOutWindowSize() * 1e6;
    double wireSpacing = 0.3;

    std::vector< double > cluster_start = {clusters[icl]->StartWire()*wireSpacing, drift*clusters[icl]->StartTick()*fromTickToNs};
    std::vector< double > cluster_end = {clusters[icl]->EndWire()*wireSpacing, drift*clusters[icl]->EndTick()*fromTickToNs};

    double cluster_length = sqrt(pow(cluster_end[0]-cluster_start[0],2)+pow(cluster_end[1]-cluster_start[1],2));
    if (cluster_length <= 0) continue;

    std::vector< double > cluster_axis = {(cluster_end[0]-cluster_start[0])/cluster_length,(cluster_end[1]-cluster_start[1])/cluster_length};

    // Build rectangle 4 x 1 cm around the cluster axis
    std::vector< std::vector < double > > points;
    buildRectangle(m_dQdxRectangleLength, m_dQdxRectangleWidth, cluster_start, cluster_axis, points);

    std::vector<double> dqdxs;

    for (auto & hit: hits) {

      // std::cout << "[PandoraLEE] Hit wire ID " << hit->WireID().Wire << std::endl;
      // std::cout << "Hit peak time " << hit->PeakTime() << std::endl;

      std::vector< double > hit_pos = {hit->WireID().Wire*wireSpacing,fromTickToNs*drift*hit->PeakTime()};

      double pitch = getPitch(shower_dir, clusters[icl]->Plane().Plane);

      // Hit within the rectangle
      if (withinRectangle(points, hit_pos) and pitch > 0) {
        double q = hit->Integral()*_gain;
        dqdxs.push_back(q/pitch);
      }
    }

    // Get the median
    size_t n = dqdxs.size() / 2;
    std::nth_element(dqdxs.begin(), dqdxs.begin()+n, dqdxs.end());
    if (n > 0) {
      std::cout << "[dQdx] Plane dQdx " << clusters[icl]->Plane().Plane << " " << dqdxs[n] << " " << dqdxs[n]*(23./1e6)/0.62 << std::endl;

      dqdx[clusters[icl]->Plane().Plane] = dqdxs[n];
    }

  }
}

// CORRECT SHOWER DIRECTION THAT SOMETIMES IS INVERTED
int lee::PandoraLEEAnalyzer::correct_direction(size_t pfp_id, const art::Event & evt) {




  auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( _pfp_producer );

  art::FindManyP<recob::SpacePoint > spcpnts_per_pfpart ( pfparticle_handle, evt, _pfp_producer );
  std::vector<art::Ptr < recob::SpacePoint > > spcpnts = spcpnts_per_pfpart.at(pfp_id);

  int direction = 1;

  if (spcpnts.size() > 0) {
    art::FindOneP< recob::Vertex > vertex_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    auto const& vertex_obj = vertex_per_pfpart.at(pfp_id);
    double vertex_xyz[3];
    vertex_obj->XYZ(vertex_xyz);
    TVector3 start_vec(vertex_xyz[0],vertex_xyz[1],vertex_xyz[2]);

    art::FindOneP< recob::Shower > shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    auto const& shower_obj = shower_per_pfpart.at(pfp_id);
    TVector3 shower_vec(shower_obj->Direction().X(),shower_obj->Direction().Y(),shower_obj->Direction().Z());

    TVector3 avg_spcpnt = average_position(spcpnts);

    TVector3 a;
    TVector3 b;
    a = shower_vec;
    b = avg_spcpnt - start_vec;
    double costheta = a.Dot(b);
    direction = costheta >= 0 ? 1 : -1;
  }

  return direction;
}

TVector3 lee::PandoraLEEAnalyzer::average_position(std::vector<art::Ptr < recob::SpacePoint > > &spcpnts) {
  double avg_x = 0;
  double avg_y = 0;
  double avg_z = 0;

  for (auto &spcpnt: spcpnts) {
    auto spcpnt_xyz = spcpnt->XYZ();
    avg_x += spcpnt_xyz[0];
    avg_y += spcpnt_xyz[1];
    avg_z += spcpnt_xyz[2];
  }

  return TVector3(avg_x/spcpnts.size(), avg_y/spcpnts.size(), avg_z/spcpnts.size());

}


double lee::PandoraLEEAnalyzer::showerEnergy(const art::Ptr<recob::Shower>& shower, const art::Event & evt)
{
  auto const& shower_handle = evt.getValidHandle< std::vector< recob::Shower > >( _pfp_producer );
  auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( _pfp_producer );
  auto const& cluster_handle = evt.getValidHandle< std::vector< recob::Cluster > >( _pfp_producer );

  art::FindOneP<recob::PFParticle> pfparticle_shower_ass(shower_handle, evt, "pandoraNu");
  const art::Ptr<recob::PFParticle> pfparticle = pfparticle_shower_ass.at(shower.key());

  art::FindManyP<recob::Cluster> cluster_pfparticle_ass(pfparticle_handle, evt, "pandoraNu");
  std::vector<art::Ptr < recob::Cluster > > clusters = cluster_pfparticle_ass.at(pfparticle.key());

  double energy[3] = {0,0,0};
  int n_hits[3] = {0,0,0};
  for (size_t icl = 0; icl < clusters.size(); icl++) {
    art::FindManyP<recob::Hit> hit_cluster_ass(cluster_handle, evt, "pandoraNu");
    std::vector<art::Ptr < recob::Hit > > hits = hit_cluster_ass.at(clusters[icl].key());

    for (size_t ihit = 0; ihit < hits.size(); ++ihit)
    {

      if (hits[ihit]->View() > 2 || hits[ihit]->View() < 0) continue;
      // 23 work function (23 eV/e- in the argon)
      // 0.62 recombination factor
      // TODO move to other function or use detector properties
      energy[hits[ihit]->View()] += hits[ihit]->Integral()*_gain*(23./1e6)/0.62;
      n_hits[hits[ihit]->View()]++;
    }
  }

  const int n = sizeof(n_hits) / sizeof(int);

  return energy[std::distance(n_hits, std::max_element(n_hits, n_hits + n))]/1000; // convert to GeV
}


double lee::PandoraLEEAnalyzer::trackEnergy(const art::Ptr<recob::Track>& track, const art::Event & evt)
{
  auto const& track_handle = evt.getValidHandle< std::vector< recob::Track > >( _pfp_producer );
  art::FindManyP<anab::Calorimetry> calo_track_ass(track_handle, evt, "pandoraNucalo");
  const std::vector<art::Ptr<anab::Calorimetry>> calos = calo_track_ass.at(track->ID());
  double E = 0;
  double Eapprox = 0;

  for (size_t ical = 0; ical < calos.size(); ++ical)
  {
    if (E != 0) continue;
    if (!calos[ical]) continue;
    if (!calos[ical]->PlaneID().isValid) continue;
    int planenum = calos[ical]->PlaneID().Plane;
    if (planenum < 0 || planenum > 2) continue;
    if (planenum != 2) continue;                           // Use informartion from collection plane only

    // Understand if the calo module flipped the track
    //double dqdx_start = (calos[ical]->dQdx())[0] + (calos[ical]->dQdx())[1] + (calos[ical]->dQdx())[2];
    //double dqdx_end   = (calos[ical]->dQdx())[calos[ical]->dQdx().size()-1] + (calos[ical]->dQdx())[calos[ical]->dQdx().size()-2] + (calos[ical]->dQdx())[calos[ical]->dQdx().size()-3];
    //bool caloFlippedTrack = dqdx_start < dqdx_end;

    double mean = 0;
    double dedx = 0;
    double prevresrange = 0;

    if (calos[ical]->ResidualRange()[0] > track->Length() / 2)
    {
      prevresrange = track->Length();
    }

    double currentresrange = 0;

    for (size_t iTrkHit = 0; iTrkHit < calos[ical]->dEdx().size(); ++iTrkHit)
    {
      dedx = calos[ical]->dEdx()[iTrkHit];
      currentresrange = calos[ical]->ResidualRange()[iTrkHit];
      if (dedx > 0 && dedx < 10)
      {
        //std::cout << "[PandoraLEE] " << dedx << "\t" << currentresrange << "\t"<< prevresrange<<std::endl;
        mean += dedx;
        E += dedx * abs(prevresrange - currentresrange);
        prevresrange = currentresrange;
      }
    }
    //std::cout << "[PandoraLEE] " << "Length: " << track->Length() << "and Energy approximation is " << mean/calos[ical]->dEdx().size()*track->Length()<< "MeV"<<std::endl;
    Eapprox = mean / calos[ical]->dEdx().size() * track->Length();
  }
  return Eapprox / 1000; // convert to GeV
}


void lee::PandoraLEEAnalyzer::measure_energy(size_t ipf, const art::Event & evt, double & energy) {


  auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( _pfp_producer );
  auto const& pfparticles(*pfparticle_handle);

  art::FindManyP<recob::Shower > showers_per_pfparticle ( pfparticle_handle, evt, _pfp_producer );
  std::vector<art::Ptr<recob::Shower>> showers = showers_per_pfparticle.at(ipf);

  for (size_t ish = 0; ish < showers.size(); ish++) {
    if (showerEnergy(showers[ish], evt) > 0) {
      energy += showerEnergy(showers[ish], evt);
    }
  }

  art::FindManyP<recob::Track > tracks_per_pfparticle ( pfparticle_handle, evt, _pfp_producer );
  std::vector<art::Ptr<recob::Track>> tracks = tracks_per_pfparticle.at(ipf);

  for (size_t itr = 0; itr < tracks.size(); itr++) {
    if (trackEnergy(tracks[itr], evt) > 0) {
      energy += trackEnergy(tracks[itr], evt);
    }
  }

  for (auto const& pfdaughter : pfparticles[ipf].Daughters()) {
    measure_energy(pfdaughter, evt, energy);
  }
}

size_t lee::PandoraLEEAnalyzer::choose_candidate(std::vector<size_t> & candidates, const art::Event & evt) {

  size_t chosen_candidate = 0;
  double most_z = -1;
  double longest_track_dir;

  for (auto const& ic : candidates) {
    std::vector<art::Ptr<recob::Track>> nu_tracks;
    std::cout << "[PandoraLEE] Candidate " << ic << std::endl;
    std::vector< size_t > pfp_tracks_id = fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary().at(ic);
    get_daughter_tracks(pfp_tracks_id, evt, nu_tracks);
    longest_track_dir = get_longest_track(nu_tracks)->StartDirection().Z();

    if (longest_track_dir > most_z) {
      chosen_candidate = ic;
      most_z = longest_track_dir;
    }

  }

  return chosen_candidate;

}

void lee::PandoraLEEAnalyzer::endSubRun(const art::SubRun& sr)
{

  _run_sr = sr.run();
  _subrun_sr = sr.subRun();

  art::Handle< sumdata::POTSummary > potListHandle;
  if (!m_isData) {
    if (sr.getByLabel("generator", potListHandle))
    _pot = potListHandle->totpot;
    else
    _pot = 0.;
  } else {
    if (sr.getByLabel("beamdata",  "bnbETOR860", potListHandle))
    _pot = potListHandle->totpot;
    else
    _pot = 0.;
  }

  myPOTTTree->Fill();

}
void lee::PandoraLEEAnalyzer::clear() {
  _interaction_type = std::numeric_limits<int>::lowest();
  _shower_open_angle.clear();
  _shower_dir_x.clear();
  _shower_dir_y.clear();
  _shower_dir_z.clear();
  _shower_dQdx.clear();
  _shower_dEdx.clear();
  _matched_tracks.clear();
  _matched_showers.clear();


  _shower_start_x.clear();
  _shower_start_y.clear();
  _shower_start_z.clear();

  _track_is_fiducial.clear();
  _shower_is_fiducial.clear();

  _track_start_x.clear();
  _track_start_y.clear();
  _track_start_z.clear();

  _track_end_x.clear();
  _track_end_y.clear();
  _track_end_z.clear();

  _predict_p.clear();
  _predict_mu.clear();
  _predict_em.clear();
  _predict_pi.clear();
  _predict_cos.clear();

  _shower_theta.clear();
  _shower_phi.clear();

  _track_dir_x.clear();
  _track_dir_y.clear();
  _track_dir_z.clear();

  _track_theta.clear();
  _track_phi.clear();

  _shower_energy.clear();
  _track_energy.clear();

  _track_length.clear();
  _track_id.clear();

  _nu_pdg = 0;

  _flash_passed = 0;
  _track_passed = 0;
  _shower_passed = 0;
  _energy = std::numeric_limits<double>::lowest();

  _true_nu_is_fiducial = std::numeric_limits<int>::lowest();
  _nu_energy = std::numeric_limits<double>::lowest();

  _n_tracks = std::numeric_limits<int>::lowest();
  _n_showers = std::numeric_limits<int>::lowest();

  _vx = std::numeric_limits<double>::lowest();
  _vy = std::numeric_limits<double>::lowest();
  _vz = std::numeric_limits<double>::lowest();

  _true_vx = std::numeric_limits<double>::lowest();
  _true_vy = std::numeric_limits<double>::lowest();
  _true_vz = std::numeric_limits<double>::lowest();

  _true_vx_sce = std::numeric_limits<double>::lowest();
  _true_vy_sce = std::numeric_limits<double>::lowest();
  _true_vz_sce = std::numeric_limits<double>::lowest();

  _nu_matched_tracks = std::numeric_limits<int>::lowest();
  _nu_matched_showers = std::numeric_limits<int>::lowest();

  _category = std::numeric_limits<int>::lowest();
  _run = std::numeric_limits<int>::lowest();
  _subrun = std::numeric_limits<int>::lowest();
  _event = std::numeric_limits<int>::lowest();
  _n_candidates = std::numeric_limits<int>::lowest();
  _n_true_nu = std::numeric_limits<int>::lowest();
  _run_sr = std::numeric_limits<int>::lowest();
  _subrun_sr = std::numeric_limits<int>::lowest();
  _n_matched = std::numeric_limits<int>::lowest();
  _pot = std::numeric_limits<double>::lowest();
  _event_passed = 0;
  _distance = std::numeric_limits<double>::lowest();

  _nu_daughters_E.clear();
  _nu_daughters_pdg.clear();

  _nu_daughters_px.clear();
  _nu_daughters_py.clear();
  _nu_daughters_pz.clear();

  _nu_daughters_vx.clear();
  _nu_daughters_vy.clear();
  _nu_daughters_vz.clear();

  _nu_daughters_endx.clear();
  _nu_daughters_endy.clear();
  _nu_daughters_endz.clear();

  _bnbweight = std::numeric_limits<int>::lowest();

}

void lee::PandoraLEEAnalyzer::analyze(art::Event const & evt)
{
  clear();


  _run = evt.run();
  _subrun = evt.subRun();
  _event = evt.id().event();

  std::cout << "[PandoraLEE] " << "RUN " << _run << " SUBRUN " << _subrun << " EVENT " << _event << std::endl;

  std::vector<size_t> nu_candidates;

  _event_passed = int(fElectronEventSelectionAlg.eventSelected(evt));

  _category = 0;
  std::vector< double > true_neutrino_vertex(3);

  if (!evt.isRealData()) {
    _gain = 200;

    // nu_e flux must be corrected by event weight

    art::InputTag eventweight_tag("eventweight");

    auto const & eventweights_handle = evt.getValidHandle<std::vector <evwgh::MCEventWeight> >(eventweight_tag);
    auto const & eventweights(*eventweights_handle);

    if (eventweights.size() > 0) {

      for(auto last : eventweights.at(0).fWeight) {
        if(last.first.find("bnbcorrection") != std::string::npos) {

          if(!std::isfinite(last.second.at(0))) {
            _bnbweight = 1;
          } else {
            _bnbweight = last.second.at(0);
          }

        }

      }

    } else { _bnbweight = 1;}


    auto const& generator_handle = evt.getValidHandle< std::vector< simb::MCTruth > >( _mctruthLabel );
    auto const& generator(*generator_handle);
    _n_true_nu = generator.size();
    _true_nu_is_fiducial = 0;
    std::vector<simb::MCParticle> nu_mcparticles;

    if (generator.size() > 0) {
      _nu_pdg = generator[0].GetNeutrino().Nu().PdgCode();
      _nu_energy = generator[0].GetNeutrino().Nu().E();
      int ccnc = generator[0].GetNeutrino().CCNC();
      if (ccnc == 1) {
        _category = k_nc;
      }

      true_neutrino_vertex[0] = generator[0].GetNeutrino().Nu().Vx();
      true_neutrino_vertex[1] = generator[0].GetNeutrino().Nu().Vy();
      true_neutrino_vertex[2] = generator[0].GetNeutrino().Nu().Vz();
      _true_vx = true_neutrino_vertex[0];
      _true_vy = true_neutrino_vertex[1];
      _true_vz = true_neutrino_vertex[2];
      _true_nu_is_fiducial = int(fElectronEventSelectionAlg.is_fiducial(true_neutrino_vertex));

      _interaction_type = generator[0].GetNeutrino().InteractionType();

      std::string _env = std::getenv("MRB_INSTALL");
      _env = _env + "/pandoraElectronFilter/v00_01_00/slf6.x86_64.e10.prof/lib/";

      SpaceChargeMicroBooNE SCE =
      SpaceChargeMicroBooNE(_env + "SCEoffsets_MicroBooNE_E273.root");

      _true_vx_sce = _true_vx-SCE.GetPosOffsets(_true_vx, _true_vy, _true_vz)[0]+0.7;
      _true_vy_sce = _true_vy+SCE.GetPosOffsets(_true_vx, _true_vy, _true_vz)[1];
      _true_vz_sce = _true_vz+SCE.GetPosOffsets(_true_vx, _true_vy, _true_vz)[2];

      if (fElectronEventSelectionAlg.is_dirt(true_neutrino_vertex)) {
        _category = k_dirt;
      }

      for (int i = 0; i < generator[0].NParticles(); i++) {
        if (generator[0].Origin() == 1) {
          nu_mcparticles.push_back(generator[0].GetParticle(i));
        }
      }
    } else {
      _category = k_cosmic;
      _nu_energy = std::numeric_limits<double>::lowest();
    }

    for (auto& mcparticle : nu_mcparticles) {
      if (mcparticle.Process() == "primary" and mcparticle.T() != 0 and mcparticle.StatusCode() == 1) {

        _nu_daughters_E.push_back(mcparticle.E());
        _nu_daughters_pdg.push_back(mcparticle.PdgCode());

        _nu_daughters_px.push_back(mcparticle.Px());
        _nu_daughters_py.push_back(mcparticle.Py());
        _nu_daughters_pz.push_back(mcparticle.Pz());

        _nu_daughters_vx.push_back(mcparticle.Vx());
        _nu_daughters_vy.push_back(mcparticle.Vy());
        _nu_daughters_vz.push_back(mcparticle.Vz());

        _nu_daughters_endx.push_back(mcparticle.EndX());
        _nu_daughters_endy.push_back(mcparticle.EndY());
        _nu_daughters_endz.push_back(mcparticle.EndZ());

      }
    }

    if (_category != k_cosmic && _category != k_dirt && _category != k_nc) {
      if (abs(_nu_pdg) == 12) {
        _category = k_nu_e;
      }
      if (abs(_nu_pdg) == 14) {
        _category = k_nu_mu;
      }
    }
  } else {
    _gain = 240;
    _category = k_data;
  }

  std::cout << "[PandoraLEE] " << "True neutrino PDG " << _nu_pdg << std::endl;
  std::cout << "[PandoraLEE] " << "Nu energy " << _nu_energy << std::endl;

  _energy = std::numeric_limits<double>::lowest();

  for (auto & i_primary : fElectronEventSelectionAlg.get_primary_indexes() ) {
    try {

      if (fElectronEventSelectionAlg.get_op_flash_indexes().at(i_primary) != -1) {
        _flash_passed = 1;
      }
      if (fElectronEventSelectionAlg.get_n_showers().at(i_primary) != 0) {
        _shower_passed = 1;
      }
      if (fElectronEventSelectionAlg.get_n_tracks().at(i_primary) != 0) {
        _track_passed = 1;
      }
    } catch (...) {
      std::cout << "[PandoraLEE] Error getting passed events" << std::endl;
    }
  }

  if (_event_passed) {

    _n_candidates = nu_candidates.size();
    for (auto & inu : fElectronEventSelectionAlg.get_primary_indexes()) {
      if (fElectronEventSelectionAlg.get_neutrino_candidate_passed().at(inu)) {
        nu_candidates.push_back(inu);
      }
    }

    std::cout << "[PandoraLEE] " << "EVENT PASSED" << std::endl;
    auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( _pfp_producer );
    //auto const& pfparticles(*pfparticle_handle);


    size_t ipf_candidate = choose_candidate(nu_candidates, evt);
    std::cout << "[PandoraLEE] " << "Neutrino candidate " << ipf_candidate << std::endl;

    //size_t ipf_candidate = 0;
    _energy = 0;

    measure_energy(ipf_candidate, evt, _energy);

    art::FindOneP< recob::Vertex > vertex_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    auto const& vertex_obj = vertex_per_pfpart.at(ipf_candidate);


    double reco_neutrino_vertex[3];
    vertex_obj->XYZ(reco_neutrino_vertex);
    _vx = reco_neutrino_vertex[0];
    _vy = reco_neutrino_vertex[1];
    _vz = reco_neutrino_vertex[2];

    if (_category != k_data) {
      TVector3 v_reco_vertex(_vx, _vy, _vz);
      TVector3 sce_true_vertex(_true_vx_sce, _true_vy_sce, _true_vz_sce);

      _distance = fElectronEventSelectionAlg.distance(v_reco_vertex, sce_true_vertex);

    }

    std::vector< size_t > pfp_tracks_id;

    try {
      pfp_tracks_id =  fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary().at(ipf_candidate);
      _n_tracks = fElectronEventSelectionAlg.get_n_tracks().at(ipf_candidate);
    } catch (...) {
      std::cout << "[PandoraLEE] Error getting associated tracks to neutrino candidate" << std::endl;
    }


    for (auto &pf_id: pfp_tracks_id) {
      _matched_tracks.push_back(std::numeric_limits<int>::lowest());

      art::FindOneP< recob::Track > track_per_pfpart(pfparticle_handle, evt, _pfp_producer);
      auto const& track_obj = track_per_pfpart.at(pf_id);

      auto const& trackVecHandle = evt.getValidHandle< std::vector< recob::Track > >( _pfp_producer );

      art::FindManyP<anab::CosmicTag> dtAssns(trackVecHandle, evt, "decisiontreeid");

      std::vector< art::Ptr<anab::CosmicTag> > dtVec = dtAssns.at(track_obj.key());

      for(auto const& dttag : dtVec)
      {
        if(dttag->CosmicType() == TAGID_P)  _predict_p.push_back(dttag->CosmicScore());
        else if(dttag->CosmicType() == TAGID_MU) _predict_mu.push_back(dttag->CosmicScore());
        else if(dttag->CosmicType() == TAGID_PI) _predict_pi.push_back(dttag->CosmicScore());
        else if(dttag->CosmicType() == TAGID_EM) _predict_em.push_back(dttag->CosmicScore());
        else if(dttag->CosmicType() == TAGID_CS) _predict_cos.push_back(dttag->CosmicScore());
      }

      std::vector<double> start_point;
      std::vector<double> end_point;

      start_point.push_back(track_obj->Start().X());
      start_point.push_back(track_obj->Start().Y());
      start_point.push_back(track_obj->Start().Z());

      end_point.push_back(track_obj->End().X());
      end_point.push_back(track_obj->End().Y());
      end_point.push_back(track_obj->End().Z());

      _track_is_fiducial.push_back(int(fElectronEventSelectionAlg.is_fiducial(start_point) && fElectronEventSelectionAlg.is_fiducial(end_point)));

      _track_energy.push_back(trackEnergy(track_obj,evt));
      std::cout << "[PandoraLEE] " << "Track energy " << trackEnergy(track_obj,evt) << std::endl;
      _track_length.push_back(track_obj->Length());
      _track_id.push_back(track_obj->ID());
      _track_dir_x.push_back(track_obj->StartDirection().X());
      _track_dir_y.push_back(track_obj->StartDirection().Y());
      _track_dir_z.push_back(track_obj->StartDirection().Z());

      _track_start_x.push_back(track_obj->Start().X());
      _track_start_y.push_back(track_obj->Start().Y());
      _track_start_z.push_back(track_obj->Start().Z());

      _track_end_x.push_back(track_obj->End().X());
      _track_end_y.push_back(track_obj->End().Y());
      _track_end_z.push_back(track_obj->End().Z());

      _track_theta.push_back(track_obj->Theta());
      _track_phi.push_back(track_obj->Phi());

    }


    std::vector< size_t > pfp_showers_id;

    try {
      pfp_showers_id = fElectronEventSelectionAlg.get_pfp_id_showers_from_primary().at(ipf_candidate);
      _n_showers = fElectronEventSelectionAlg.get_n_showers().at(ipf_candidate);
    } catch (...) {
      std::cout << "[PandoraLEE] Error getting associated showers to neutrino candidate" << std::endl;
    }

    for (auto &pf_id: pfp_showers_id) {

      std::vector<double> dqdx(3,std::numeric_limits<double>::lowest());
      std::vector<double> dedx(3,std::numeric_limits<double>::lowest());

      dQdx(pf_id, evt, dqdx);
      dEdxFromdQdx(dedx, dqdx);

      _matched_showers.push_back(std::numeric_limits<int>::lowest());

      _shower_dQdx.push_back(dqdx);
      _shower_dEdx.push_back(dedx);

      int direction = correct_direction(pf_id, evt);
      //std::cout << "[PandoraLEE] " << "Correct direction " << direction << std::endl;
      art::FindOneP< recob::Shower > shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);
      auto const& shower_obj = shower_per_pfpart.at(pf_id);

      TVector3 correct_dir(direction*shower_obj->Direction().X(),direction*shower_obj->Direction().Y(),direction*shower_obj->Direction().Z());

      _shower_dir_x.push_back(correct_dir.X());
      _shower_dir_y.push_back(correct_dir.Y());
      _shower_dir_z.push_back(correct_dir.Z());

      _shower_open_angle.push_back(shower_obj->OpenAngle());

      std::vector<double> correct_dir_v;
      correct_dir_v.resize(3);
      correct_dir_v[0] = correct_dir.X();
      correct_dir_v[1] = correct_dir.Y();
      correct_dir_v[2] = correct_dir.Z();
      double shower_length = shower_obj->Length();

      std::vector<double> start_point;
      std::vector<double> end_point;
      start_point.resize(3);
      end_point.resize(3);
      for (int ix = 0; ix < 3; ix++)
      {
        start_point[ix] = shower_obj->ShowerStart()[ix];
        end_point[ix] = shower_obj->ShowerStart()[ix] + shower_length * correct_dir_v[ix];
      }

      _shower_is_fiducial.push_back(int(fElectronEventSelectionAlg.is_fiducial(start_point) && fElectronEventSelectionAlg.is_fiducial(end_point)));

      _shower_start_x.push_back(shower_obj->ShowerStart().X());
      _shower_start_y.push_back(shower_obj->ShowerStart().Y());
      _shower_start_z.push_back(shower_obj->ShowerStart().Z());

      _shower_phi.push_back(correct_dir.Phi());
      _shower_theta.push_back(correct_dir.Theta());

      _shower_energy.push_back(showerEnergy(shower_obj, evt));
      std::cout << "[PandoraLEE] Shower energy " << showerEnergy(shower_obj, evt) << std::endl;
      for (size_t i = 0; i < 3; i++) {
        std::cout << "[PandoraLEE] Shower obj energy Plane " << i << " " << shower_obj->Energy()[i] << std::endl;
      }
    }

    lar_pandora::MCParticlesToPFParticles matchedParticles;    // This is a map: MCParticle to matched PFParticle
    lar_pandora::MCParticlesToHits        matchedParticleHits;

    // --- Do the matching
    fElectronEventSelectionAlg.GetRecoToTrueMatches(evt, _pfp_producer, _spacepointLabel, _geantModuleLabel, _hitfinderLabel, matchedParticles, matchedParticleHits);

    art::ServiceHandle<cheat::BackTracker> bt;

    std::vector< art::Ptr<recob::PFParticle> > neutrino_pf;
    std::vector< art::Ptr<recob::PFParticle> > cosmic_pf;

    std::vector<int> neutrino_pdg;
    std::vector<int> cosmic_pdg;

    for (lar_pandora::MCParticlesToPFParticles::const_iterator iter1 = matchedParticles.begin(), iterEnd1 = matchedParticles.end();
    iter1 != iterEnd1; ++iter1) {

      art::Ptr<simb::MCParticle>  mc_par = iter1->first;   // The MCParticle
      art::Ptr<recob::PFParticle> pf_par = iter1->second;  // The matched PFParticle

      const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());

      if (mc_truth->Origin() == simb::kBeamNeutrino) {
        //std::cout << "[PandoraLEE] " << "Matched neutrino" << std::endl;
        //std::cout << "[PandoraLEE] " << "Pf PDG: " << pf_par->PdgCode() << " MC PDG: " << mc_par->PdgCode() << std::endl;
        neutrino_pf.push_back(pf_par);
        neutrino_pdg.push_back(mc_par->PdgCode());
      }

      if (mc_truth->Origin() == simb::kCosmicRay) {
        cosmic_pf.push_back(pf_par);
        cosmic_pdg.push_back(mc_par->PdgCode());
      }

    }

    _n_matched = neutrino_pf.size();

    _nu_matched_showers = 0;
    _nu_matched_tracks = 0;

    bool shower_cr_found = false;

    for (size_t ish = 0; ish < pfp_showers_id.size(); ish++) {
      for (size_t ipf = 0; ipf < cosmic_pf.size(); ipf++ ) {
        if (pfp_showers_id[ish] == cosmic_pf[ipf].key()) {
          shower_cr_found = true;
          _matched_showers[ish] = cosmic_pdg[ipf];
        }
      }

      for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++ ) {
        if (pfp_showers_id[ish] == neutrino_pf[ipf].key()) {
          _nu_matched_showers++;
          _matched_showers[ish] = neutrino_pdg[ipf];
        }
      }

      std::cout << "[PandoraLEE] " << "Shower PFP " << pfp_showers_id[ish] << std::endl;
      std::cout << "[PandoraLEE] " << "Neutrino? " << _nu_matched_showers << std::endl;
      std::cout << "[PandoraLEE] " << "Cosmic? " << shower_cr_found << std::endl;
      if (!shower_cr_found && _nu_matched_showers == 0  && _category != k_dirt) {
        _category = k_other;
        std::cout << "[PandoraLEE] " << "***NOT NEUTRINO NOR COSMIC***" << std::endl;
      }
      if (shower_cr_found && _nu_matched_showers > 0) {
        _category = k_mixed;
        std::cout << "[PandoraLEE] " << "***MIXED COSMIC/NEUTRINO***" << std::endl;
      }
    }

    bool track_cr_found = false;

    for (size_t itr = 0; itr < pfp_tracks_id.size(); itr++) {

      for (size_t ipf = 0; ipf < cosmic_pf.size(); ipf++ ) {
        if (pfp_tracks_id[itr] == cosmic_pf[ipf].key()) {
          track_cr_found = true;
          _matched_tracks[itr] = cosmic_pdg[ipf];
        }
      }

      for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++ ) {
        if (pfp_tracks_id[itr] == neutrino_pf[ipf].key()) {
          _nu_matched_tracks++;
          _matched_tracks[itr] = neutrino_pdg[ipf];
        }
      }


      std::cout << "[PandoraLEE] " << "Track PFP " << pfp_tracks_id[itr] << std::endl;
      std::cout << "[PandoraLEE] " << "Neutrino? " << _nu_matched_tracks << std::endl;
      std::cout << "[PandoraLEE] " << "Cosmic? " << track_cr_found << std::endl;
      if (!track_cr_found && _nu_matched_tracks == 0 && _category != k_dirt) {
        _category = k_other;
        std::cout << "[PandoraLEE] " << "***NOT NEUTRINO NOR COSMIC***" << std::endl;
      }
      if(track_cr_found && _nu_matched_tracks > 0) {
        _category = k_mixed;
        std::cout << "[PandoraLEE] " << "***MIXED COSMIC/NEUTRINO***" << std::endl;
      }
    }

    if ((track_cr_found || shower_cr_found) && _category != k_mixed) _category = k_cosmic;
    std::cout << "[PandoraLEE] " << "Category " << _category << std::endl;

  } else {
    std::cout << "[PandoraLEE] " << "EVENT NOT PASSED" << std::endl;
  }

  myTTree->Fill();
  std::cout << "[PandoraLEE] " << "END ANALYZER" << std::endl;

} // end analyze function

//------------------------------------------------------------------------------------------------------------------------------------


void lee::PandoraLEEAnalyzer::reconfigure(fhicl::ParameterSet const & pset)
{

  //TODO: add an external fcl file to change configuration
  //add what you want to read, and default values of your labels etc. example:
  //  m_particleLabel = pset.get<std::string>("PFParticleModule","pandoraNu");
  fElectronEventSelectionAlg.reconfigure(pset.get<fhicl::ParameterSet>("ElectronSelectionAlg"));

  m_printDebug = pset.get<bool>("PrintDebug", false);

  m_fidvolXstart = pset.get<double>("fidvolXstart", 10);
  m_fidvolXend = pset.get<double>("fidvolXstart", 10);

  m_fidvolYstart = pset.get<double>("fidvolYstart", 20);
  m_fidvolYend = pset.get<double>("fidvolYend", 20);

  m_fidvolZstart = pset.get<double>("fidvolZstart", 10);
  m_fidvolZend = pset.get<double>("fidvolZend", 50);

  m_isData = pset.get<bool>("isData", false);
  m_dQdxRectangleWidth = pset.get<double>("dQdxRectangleWidth", 1);
  m_dQdxRectangleLength = pset.get<double>("dQdxRectangleLength", 4);

}

//---------------------------------------------------------------------------------------------------------------------------
//add other functions here

DEFINE_ART_MODULE(lee::PandoraLEEAnalyzer)
