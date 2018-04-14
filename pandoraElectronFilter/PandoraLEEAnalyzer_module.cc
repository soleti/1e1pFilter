////////////////////////////////////////////////////////////////////////
// Class:       PandoraLEEAnalyzer
// Module Type: analyzer
// File:        PandoraLEEAnalyzer_module.cc
//
// by Roberto and Wouter using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "PandoraLEEAnalyzer.h"

lee::PandoraLEEAnalyzer::PandoraLEEAnalyzer(fhicl::ParameterSet const &pset)
    : EDAnalyzer(pset) // ,
// More initializers here.
{

    // create output tree
    art::ServiceHandle<art::TFileService> tfs;
    // myTFile = new TFile("PandoraLEEAnalyzerOutput.root", "RECREATE");
    myTTree = tfs->make<TTree>("pandoratree", "PandoraAnalysis Tree");

    myPOTTTree = tfs->make<TTree>("pot", "POT Tree");

    myPOTTTree->Branch("run", &_run_sr, "run/i");
    myPOTTTree->Branch("subrun", &_subrun_sr, "subrun/i");
    myPOTTTree->Branch("pot", &_pot, "pot/d");

    myTTree->Branch("category", &_category, "category/i");

    myTTree->Branch("n_tracks", &_n_tracks, "n_tracks/i");
    myTTree->Branch("n_showers", &_n_showers, "n_showers/i");

    myTTree->Branch("ccnc", &_ccnc, "ccnc/i");
    myTTree->Branch("qsqr", &_qsqr, "qsqr/d");
    myTTree->Branch("theta", &_theta, "theta/d");
    myTTree->Branch("vx", &_vx, "vx/d");
    myTTree->Branch("vy", &_vy, "vy/d");
    myTTree->Branch("vz", &_vz, "vz/d");
    myTTree->Branch("fiducial", &_fiducial, "fiducial/i");

    myTTree->Branch("true_vx", &_true_vx, "true_vx/d");
    myTTree->Branch("true_vy", &_true_vy, "true_vy/d");
    myTTree->Branch("true_vz", &_true_vz, "true_vz/d");
    myTTree->Branch("true_nu_fiducial", &_true_nu_is_fiducial, "true_nu_fiducial/i");

    myTTree->Branch("true_shower_x_sce", "std::vector< double >", &_true_shower_x_sce);
    myTTree->Branch("true_shower_y_sce", "std::vector< double >", &_true_shower_y_sce);
    myTTree->Branch("true_shower_z_sce", "std::vector< double >", &_true_shower_z_sce);
    myTTree->Branch("true_shower_pdg", "std::vector< int >", &_true_shower_pdg);
    myTTree->Branch("true_shower_depE", "std::vector< double >", &_true_shower_depE);

    myTTree->Branch("true_vx_sce", &_true_vx_sce, "true_vx_sce/d");
    myTTree->Branch("true_vy_sce", &_true_vy_sce, "true_vy_sce/d");
    myTTree->Branch("true_vz_sce", &_true_vz_sce, "true_vz_sce/d");
    myTTree->Branch("lepton_E", &_lepton_E, "lepton_E/d");
    myTTree->Branch("lepton_theta", &_lepton_theta, "lepton_theta/d");

    myTTree->Branch("nu_E", &_nu_energy, "nu_E/d");
    myTTree->Branch("passed", &_event_passed, "passed/I");
    myTTree->Branch("numu_passed", &_numu_passed, "numu_passed/I");
    myTTree->Branch("numu_cuts", &_numu_cuts, "numu_cuts/I");

    myTTree->Branch("n_true_nu", &_n_true_nu, "n_true_nu/i");
    myTTree->Branch("distance", &_distance, "distance/d");

    myTTree->Branch("n_matched", &_n_matched, "n_matched/i");
    myTTree->Branch("nu_matched_tracks", &_nu_matched_tracks,
                    "nu_matched_tracks/i");
    myTTree->Branch("nu_matched_showers", &_nu_matched_showers,
                    "nu_matched_showers/i");

    myTTree->Branch("nu_daughters_pdg", "std::vector< int >", &_nu_daughters_pdg);
    myTTree->Branch("nu_daughters_E", "std::vector< double >", &_nu_daughters_E);

    myTTree->Branch("nu_daughters_vx", "std::vector< double >",
                    &_nu_daughters_vx);
    myTTree->Branch("nu_daughters_vy", "std::vector< double >",
                    &_nu_daughters_vy);
    myTTree->Branch("nu_daughters_vz", "std::vector< double >",
                    &_nu_daughters_vz);

    myTTree->Branch("nu_daughters_endx", "std::vector< double >",
                    &_nu_daughters_endx);
    myTTree->Branch("nu_daughters_endy", "std::vector< double >",
                    &_nu_daughters_endy);
    myTTree->Branch("nu_daughters_endz", "std::vector< double >",
                    &_nu_daughters_endz);

    myTTree->Branch("nu_daughters_px", "std::vector< double >",
                    &_nu_daughters_px);
    myTTree->Branch("nu_daughters_py", "std::vector< double >",
                    &_nu_daughters_py);
    myTTree->Branch("nu_daughters_pz", "std::vector< double >",
                    &_nu_daughters_pz);

    myTTree->Branch("nu_track_ids", "std::vector< size_t >",
                    &_nu_track_ids);
    myTTree->Branch("nu_shower_ids", "std::vector< size_t >",
                    &_nu_shower_ids);

    myTTree->Branch("nu_shower_daughters", "std::vector< std::vector< int > >",
                    &_nu_shower_daughters);
    myTTree->Branch("nu_track_daughters", "std::vector< std::vector< int > >",
                    &_nu_track_daughters);

    myTTree->Branch("event", &_event, "event/i");
    myTTree->Branch("run", &_run, "run/i");
    myTTree->Branch("subrun", &_subrun, "subrun/i");

    myTTree->Branch("bnbweight", &_bnbweight, "bnbweight/d");

    myTTree->Branch("chosen_candidate", &_chosen_candidate, "chosen_candidate/i");
    myTTree->Branch("candidate_pdg", &_candidate_pdg, "candidate_pdg/i");
    myTTree->Branch("n_primaries", &_n_primaries, "n_primaries/i");

    myTTree->Branch("primary_indexes", "std::vector< int >", &_primary_indexes);
    myTTree->Branch("number_tracks", "std::vector< int >", &_number_tracks);
    myTTree->Branch("number_showers", "std::vector< int >", &_number_showers);

    myTTree->Branch("flash_time", "std::vector< double >", &_flash_time);
    myTTree->Branch("flash_PE", "std::vector< double >", &_flash_PE);
    myTTree->Branch("flash_x", &_flash_x, "flash_x/d");
    myTTree->Branch("TPC_x", &_TPC_x, "TPC_x/d");

    myTTree->Branch("flash_passed", "std::vector< int >", &_flash_passed);
    myTTree->Branch("track_passed", "std::vector< int >", &_track_passed);
    myTTree->Branch("shower_passed", "std::vector< int >", &_shower_passed);

    myTTree->Branch("shower_dir_x", "std::vector< double >", &_shower_dir_x);
    myTTree->Branch("shower_dir_y", "std::vector< double >", &_shower_dir_y);
    myTTree->Branch("shower_dir_z", "std::vector< double >", &_shower_dir_z);

    myTTree->Branch("shower_start_x", "std::vector< double >", &_shower_start_x);
    myTTree->Branch("shower_start_y", "std::vector< double >", &_shower_start_y);
    myTTree->Branch("shower_start_z", "std::vector< double >", &_shower_start_z);

    myTTree->Branch("shower_theta", "std::vector< double >", &_shower_theta);
    myTTree->Branch("shower_phi", "std::vector< double >", &_shower_phi);

    myTTree->Branch("shower_energy_hits", "std::vector< std::vector< double > >", &_shower_energy_hits);
    myTTree->Branch("shower_energy_cali", "std::vector< std::vector< float > >", &_shower_energy_cali);
    myTTree->Branch("shower_energy_product", "std::vector< std::vector< double > >", &_shower_energy_product);
    myTTree->Branch("track_energy_dedx", "std::vector< double >", &_track_energy_dedx);
    myTTree->Branch("track_energy_hits", "std::vector< std::vector< double > >", &_track_energy_hits);
    myTTree->Branch("track_energy_cali", "std::vector< std::vector< float > >", &_track_energy_cali);

    myTTree->Branch("track_dir_x", "std::vector< double >", &_track_dir_x);
    myTTree->Branch("track_dir_y", "std::vector< double >", &_track_dir_y);
    myTTree->Branch("track_dir_z", "std::vector< double >", &_track_dir_z);

    myTTree->Branch("track_start_x", "std::vector< double >", &_track_start_x);
    myTTree->Branch("track_start_y", "std::vector< double >", &_track_start_y);
    myTTree->Branch("track_start_z", "std::vector< double >", &_track_start_z);

    myTTree->Branch("track_end_x", "std::vector< double >", &_track_end_x);
    myTTree->Branch("track_end_y", "std::vector< double >", &_track_end_y);
    myTTree->Branch("track_end_z", "std::vector< double >", &_track_end_z);

    myTTree->Branch("track_theta", "std::vector< double >", &_track_theta);
    myTTree->Branch("track_phi", "std::vector< double >", &_track_phi);

    myTTree->Branch("track_len", "std::vector< double >", &_track_length);
    myTTree->Branch("track_id", "std::vector< double >", &_track_id);

    myTTree->Branch("track_pidchi", "std::vector< double >", &_track_pidchi);
    myTTree->Branch("track_pidchipr", "std::vector< double >", &_track_pidchipr);
    myTTree->Branch("track_pidchika", "std::vector< double >", &_track_pidchika);
    myTTree->Branch("track_pidchipi", "std::vector< double >", &_track_pidchipi);
    myTTree->Branch("track_pidchimu", "std::vector< double >", &_track_pidchimu);
    myTTree->Branch("track_pida", "std::vector< double >", &_track_pida);
    myTTree->Branch("track_res_mean", "std::vector< double >", &_track_res_mean);
    myTTree->Branch("track_res_std", "std::vector< double >", &_track_res_std);

    myTTree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/i");

    myTTree->Branch("predict_p", "std::vector< double >", &_predict_p);
    myTTree->Branch("predict_mu", "std::vector< double >", &_predict_mu);
    myTTree->Branch("predict_pi", "std::vector< double >", &_predict_pi);
    myTTree->Branch("predict_em", "std::vector< double >", &_predict_em);
    myTTree->Branch("predict_cos", "std::vector< double >", &_predict_cos);

    myTTree->Branch("interaction_type", &_interaction_type, "interaction_type/i");

    myTTree->Branch("shower_dQdx", "std::vector< std::vector< double > >", &_shower_dQdx);
    myTTree->Branch("shower_dQdx_cali", "std::vector< std::vector< float > >", &_shower_dQdx_cali);
    myTTree->Branch("shower_dEdx", "std::vector< std::vector< double > >", &_shower_dEdx);
    myTTree->Branch("shower_dQdx_wires", "std::vector< std::vector< int > >", &_shower_dQdx_wires);

    myTTree->Branch("track_dQdx", "std::vector< std::vector< double > >", &_track_dQdx);
    myTTree->Branch("track_dQdx_cali", "std::vector< std::vector< float > >", &_track_dQdx_cali);
    myTTree->Branch("track_dEdx", "std::vector< std::vector< double > >", &_track_dEdx);
    myTTree->Branch("track_dQdx_wires", "std::vector< std::vector< int > >", &_track_dQdx_wires);

    myTTree->Branch("shower_open_angle", "std::vector< double >",
                    &_shower_open_angle);

    myTTree->Branch("shower_length", "std::vector< double >",
                    &_shower_length);

    myTTree->Branch("shower_dQdx_hits", "std::vector< std::vector< double > >",
                    &_shower_dQdx_hits);
    myTTree->Branch("shower_dEdx_hits", "std::vector< std::vector< double > >",
                    &_shower_dEdx_hits);

    myTTree->Branch("track_dQdx_hits", "std::vector< std::vector< double > >",
                    &_track_dQdx_hits);
    myTTree->Branch("track_dEdx_hits", "std::vector< std::vector< double > >",
                    &_track_dEdx_hits);

    myTTree->Branch("matched_tracks", "std::vector< int >",
                    &_matched_tracks);
    myTTree->Branch("matched_tracks_energy", "std::vector< double >",
                    &_matched_tracks_energy);
    myTTree->Branch("matched_tracks_process", "std::vector< std::string >",
                    &_matched_tracks_process);

    myTTree->Branch("matched_showers", "std::vector< int >",
                    &_matched_showers);
    myTTree->Branch("matched_showers_process", "std::vector< std::string >",
                    &_matched_showers_process);
    myTTree->Branch("matched_showers_energy", "std::vector< double >",
                    &_matched_showers_energy);

    myTTree->Branch("shower_pca", "std::vector< double >",
                    &_shower_pca);

    myTTree->Branch("shower_nhits_cluster", "std::vector< std::vector<int> >",
                    &_shower_nhits_cluster);
    myTTree->Branch("shower_nhits_spacepoint", "std::vector< std::vector<int> >",
                    &_shower_nhits_spacepoint);

    myTTree->Branch("track_pca", "std::vector< double >",
                    &_track_pca);

    myTTree->Branch("track_nhits_cluster", "std::vector< std::vector<int> >",
                    &_track_nhits_cluster);
    myTTree->Branch("track_nhits_spacepoint", "std::vector< std::vector<int> >",
                    &_track_nhits_spacepoint);

    // Do not store all the hits/spacepoints any longer.
    myTTree->Branch("shower_sp_x", "std::vector< float >", &_shower_sp_x);
    myTTree->Branch("shower_sp_y", "std::vector< float >", &_shower_sp_y);
    myTTree->Branch("shower_sp_z", "std::vector< float >", &_shower_sp_z);
    myTTree->Branch("shower_sp_int", "std::vector< float >", &_shower_sp_int);

    // Features from the feature helper:
    myTTree->Branch("true_1eX_signal", &_true_1eX_signal, "true_1eX_signal/i");
    myTTree->Branch("track_bdt_precut", &_track_bdt_precut, "track_bdt_precut/i");

    myTTree->Branch("shower_vtxdistance", "std::vector< double >", &_shower_vtxdistance);
    myTTree->Branch("track_vtxdistance", "std::vector< double >", &_track_vtxdistance);

    myTTree->Branch("track_maxangle", "std::vector< double >", &_track_maxangle);
    myTTree->Branch("shower_maxangle", "std::vector< double >", &_shower_maxangle);

    myTTree->Branch("track_daughter", "std::vector< int >", &_track_daughter);
    myTTree->Branch("track_is_daughter", "std::vector< int >", &_track_is_daughter);
    myTTree->Branch("shower_daughter", "std::vector< int >", &_shower_daughter);
    myTTree->Branch("shower_is_daughter", "std::vector< int >", &_shower_is_daughter);

    myTTree->Branch("shower_cle", "std::vector< int >", &_shower_cle);
    myTTree->Branch("track_cle", "std::vector< int >", &_track_cle);

    myTTree->Branch("shower_fidvol_ratio", "std::vector<float>", &_shower_fidvol_ratio);
    myTTree->Branch("shower_spacepoint_dqdx_ratio", "std::vector<float>", &_shower_spacepoint_dqdx_ratio);
    myTTree->Branch("track_spacepoint_dqdx_ratio", "std::vector<float>", &_track_spacepoint_dqdx_ratio);
    myTTree->Branch("track_containment", "std::vector<int>", &_track_containment);
 
    myTTree->Branch("flash_PE_max", &_flash_PE_max, "flash_PE_max/d");
    myTTree->Branch("flash_time_max", &_flash_time_max, "flash_time_max/d");

    myTTree->Branch("chargecenter_x", &_chargecenter_x, "chargecenter_x/F");
    myTTree->Branch("chargecenter_y", &_chargecenter_y, "chargecenter_y/F");
    myTTree->Branch("chargecenter_z", &_chargecenter_z, "chargecenter_z/F");
    myTTree->Branch("total_spacepoint_containment", &_total_spacepoint_containment, "total_spacepoint_containment/F");

    myTTree->Branch("shower_dedx_hits_w", "std::vector<int>", &_shower_dedx_hits_w);
    myTTree->Branch("shower_dedx_w", "std::vector<float>", &_shower_dedx_w);
    myTTree->Branch("shower_dedx_best_w", "std::vector<float>", &_shower_dedx_best_w);
    myTTree->Branch("track_dedx_hits_w", "std::vector<int>", &_track_dedx_hits_w);
    myTTree->Branch("track_dedx_w", "std::vector<float>", &_track_dedx_w);
    myTTree->Branch("track_dedx_best_w", "std::vector<float>", &_track_dedx_best_w);

    myTTree->Branch("shower_hits_w", "std::vector<int>", &_shower_hits_w);
    myTTree->Branch("shower_energy_w", "std::vector<float>", &_shower_energy_w);
    myTTree->Branch("shower_hitsratio_w", "std::vector<float>", &_shower_hitsratio_w);
    myTTree->Branch("track_hits_w", "std::vector<int>", &_track_hits_w);
    myTTree->Branch("track_energy_w", "std::vector<float>", &_track_energy_w);
    myTTree->Branch("track_hitsratio_w", "std::vector<float>", &_track_hitsratio_w);

    this->reconfigure(pset);
}

lee::PandoraLEEAnalyzer::~PandoraLEEAnalyzer()
{
    std::cout << "[PandoraLEE] "
              << "End!" << std::endl;
}

art::Ptr<recob::Shower> lee::PandoraLEEAnalyzer::get_most_energetic_shower(
    std::vector<art::Ptr<recob::Shower>> &showers)
{
    art::Ptr<recob::Shower> most_energetic_shower;

    double max_energy = std::numeric_limits<double>::lowest();
    for (auto const &shower : showers)
    {
        if (shower->Energy()[shower->best_plane()] > max_energy)
        {
            most_energetic_shower = shower;
            max_energy = shower->Energy()[shower->best_plane()];
        }
    }
    return most_energetic_shower;
}

void lee::PandoraLEEAnalyzer::endSubRun(const art::SubRun &sr)
{

    _run_sr = sr.run();
    _subrun_sr = sr.subRun();

    art::Handle<sumdata::POTSummary> potListHandle;
    if (!m_isData)
    {
        if (sr.getByLabel("generator", potListHandle))
            _pot = potListHandle->totpot;
        else
            _pot = 0.;
    }
    else
    {
        if (sr.getByLabel("beamdata", "bnbETOR860", potListHandle))
            _pot = potListHandle->totpot;
        else
            _pot = 0.;
    }

    myPOTTTree->Fill();
}
void lee::PandoraLEEAnalyzer::clear()
{
    _track_res_mean.clear();
    _track_res_std.clear();
    _qsqr = std::numeric_limits<double>::lowest();
    _theta = std::numeric_limits<double>::lowest();

    _ccnc = std::numeric_limits<int>::lowest();

    _track_pida.clear();
    _track_pidchipr.clear();
    _track_pidchika.clear();
    _track_pidchipi.clear();
    _track_pidchimu.clear();
    _track_pidchi.clear();

    _interaction_type = std::numeric_limits<int>::lowest();

    _shower_pca.clear();
    _track_pca.clear();

    _shower_nhits_cluster.clear();
    _shower_nhits_spacepoint.clear();

    _track_nhits_cluster.clear();
    _track_nhits_spacepoint.clear();

    _matched_tracks.clear();
    _matched_tracks_process.clear();
    _matched_tracks_energy.clear();

    _matched_showers.clear();
    _matched_showers_process.clear();
    _matched_showers_energy.clear();

    _shower_dQdx_hits.clear();
    _shower_dEdx_hits.clear();
    _shower_dQdx.clear();
    _shower_dQdx_cali.clear();
    _shower_dQdx_wires.clear();
    _shower_dEdx.clear();

    _track_dQdx_hits.clear();
    _track_dEdx_hits.clear();
    _track_dQdx.clear();
    _track_dQdx_cali.clear();
    _track_dQdx_wires.clear();
    _track_dEdx.clear();

    _shower_sp_x.clear();
    _shower_sp_y.clear();
    _shower_sp_z.clear();
    _shower_sp_int.clear();

    _shower_open_angle.clear();
    _shower_length.clear();
    _shower_dir_x.clear();
    _shower_dir_y.clear();
    _shower_dir_z.clear();

    _nu_track_ids.clear();
    _nu_shower_ids.clear();
    _nu_track_daughters.clear();
    _nu_shower_daughters.clear();

    _shower_start_x.clear();
    _shower_start_y.clear();
    _shower_start_z.clear();

    _track_start_x.clear();
    _track_start_y.clear();
    _track_start_z.clear();

    _true_shower_x_sce.clear();
    _true_shower_y_sce.clear();
    _true_shower_z_sce.clear();
    _true_shower_pdg.clear();
    _true_shower_depE.clear();

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

    _shower_energy_hits.clear();
    _shower_energy_cali.clear();
    _shower_energy_product.clear();
    _track_energy_dedx.clear();
    _track_energy_hits.clear();
    _track_energy_cali.clear();

    _track_length.clear();
    _track_id.clear();

    _nu_pdg = 0;

    _flash_passed.clear();
    _track_passed.clear();
    _shower_passed.clear();
    _primary_indexes.clear();
    _number_tracks.clear();
    _number_showers.clear();

    _flash_PE.clear();
    _flash_time.clear();

    _chosen_candidate = std::numeric_limits<int>::lowest();
    _candidate_pdg = std::numeric_limits<int>::lowest();
    _n_primaries = 0;

    _nu_energy = std::numeric_limits<double>::lowest();

    _n_tracks = std::numeric_limits<int>::lowest();
    _n_showers = std::numeric_limits<int>::lowest();

    _vx = std::numeric_limits<double>::lowest();
    _vy = std::numeric_limits<double>::lowest();
    _vz = std::numeric_limits<double>::lowest();
    _fiducial = std::numeric_limits<int>::lowest();

    _true_vx = std::numeric_limits<double>::lowest();
    _true_vy = std::numeric_limits<double>::lowest();
    _true_vz = std::numeric_limits<double>::lowest();
    _true_nu_is_fiducial = std::numeric_limits<int>::lowest();
    _lepton_E = std::numeric_limits<double>::lowest();
    _lepton_theta = std::numeric_limits<double>::lowest();

    _true_vx_sce = std::numeric_limits<double>::lowest();
    _true_vy_sce = std::numeric_limits<double>::lowest();
    _true_vz_sce = std::numeric_limits<double>::lowest();

    _nu_matched_tracks = std::numeric_limits<int>::lowest();
    _nu_matched_showers = std::numeric_limits<int>::lowest();

    _category = std::numeric_limits<int>::lowest();
    _run = std::numeric_limits<int>::lowest();
    _subrun = std::numeric_limits<int>::lowest();
    _event = std::numeric_limits<int>::lowest();
    _n_true_nu = std::numeric_limits<int>::lowest();
    _run_sr = std::numeric_limits<int>::lowest();
    _subrun_sr = std::numeric_limits<int>::lowest();
    _n_matched = std::numeric_limits<int>::lowest();
    _pot = std::numeric_limits<double>::lowest();
    _event_passed = 0;
    _numu_cuts = 0;
    _numu_passed = 0;
    _distance = std::numeric_limits<double>::lowest();

    _flash_x = std::numeric_limits<double>::lowest();
    _TPC_x = std::numeric_limits<double>::lowest();

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

    // Features from the featurehelper
    _true_1eX_signal = 0;
    _track_bdt_precut = 0;
    _shower_maxangle.clear();
    _track_maxangle.clear();
    _shower_vtxdistance.clear();
    _track_vtxdistance.clear();
    _shower_cle.clear();
    _track_cle.clear();
    _track_daughter.clear();
    _track_is_daughter.clear();
    _shower_daughter.clear();
    _shower_is_daughter.clear();
    _shower_fidvol_ratio.clear();
    _shower_spacepoint_dqdx_ratio.clear();
    _track_spacepoint_dqdx_ratio.clear();
    _track_containment.clear();
    _flash_PE_max = std::numeric_limits<double>::lowest();
    _flash_time_max = std::numeric_limits<double>::lowest();
    _chargecenter_x = std::numeric_limits<float>::lowest();
    _chargecenter_y = std::numeric_limits<float>::lowest();
    _chargecenter_z = std::numeric_limits<float>::lowest();
    _total_spacepoint_containment = std::numeric_limits<float>::lowest();

    _shower_dedx_hits_w.clear();
    _shower_dedx_w.clear();
    _shower_dedx_best_w.clear();
    _track_dedx_hits_w.clear();
    _track_dedx_w.clear();
    _track_dedx_best_w.clear();
    _shower_energy_w.clear();
    _shower_hitsratio_w.clear();
    _shower_hits_w.clear();
    _track_energy_w.clear();
    _track_hitsratio_w.clear();
    _track_hits_w.clear();
}

void lee::PandoraLEEAnalyzer::categorizePFParticles(
    art::Event const &evt,
    std::vector<int> &neutrino_pdg,
    std::vector<std::string> &neutrino_process,
    std::vector<double> &neutrino_energy,
    std::vector<art::Ptr<recob::PFParticle>> &neutrino_pf,

    std::vector<int> &cosmic_pdg,
    std::vector<std::string> &cosmic_process,
    std::vector<double> &cosmic_energy,
    std::vector<art::Ptr<recob::PFParticle>> &cosmic_pf)
{

    lar_pandora::PFParticlesToMCParticles matchedParticles;

    std::cout << "[PandoraLEE] Before configure " << std::endl;

    pandoraHelper.Configure(evt, m_pfp_producer, m_spacepointLabel,
                            m_hitfinderLabel, _geantModuleLabel, _mcpHitAssLabel);

    std::cout << "[PandoraLEE] Before GetRecoToTrueMatches " << std::endl;

    pandoraHelper.GetRecoToTrueMatches(matchedParticles);
    std::cout << "[PandoraLEE] Matched size " << matchedParticles.size() << std::endl;

    // art::ServiceHandle<cheat::BackTracker> bt;

    for (lar_pandora::PFParticlesToMCParticles::const_iterator
             iter = matchedParticles.begin(),
             iterEnd = matchedParticles.end();
         iter != iterEnd; ++iter)
    {

        art::Ptr<simb::MCParticle> mc_par = iter->second; // The MCParticle
        art::Ptr<recob::PFParticle> pf_par = iter->first; // The matched PFParticle

        const auto mc_truth =
            pandoraHelper.TrackIDToMCTruth(evt, "largeant", mc_par->TrackId());

        if (mc_truth->Origin() == simb::kBeamNeutrino)
        {
            neutrino_pf.push_back(pf_par);
            neutrino_pdg.push_back(mc_par->PdgCode());
            if (mc_par->EndProcess() != NULL)
            {
                neutrino_process.push_back(mc_par->EndProcess());
            }
            neutrino_energy.push_back(mc_par->E());
        }

        if (mc_truth->Origin() == simb::kCosmicRay || mc_truth->Origin() == simb::kUnknown)
        {
            cosmic_pf.push_back(pf_par);
            cosmic_pdg.push_back(mc_par->PdgCode());
            if (mc_par->EndProcess() != NULL)
            {
                cosmic_process.push_back(mc_par->EndProcess());
            }
            cosmic_energy.push_back(mc_par->E());
        }
    }
}

// Stupid function to cast arrays of size_t to arrays of integers because root tree branches do not like long unsigned ints seemingly
std::vector<int> lee::PandoraLEEAnalyzer::vectorCast(std::vector<long unsigned int> vec_size_t)
{

    std::vector<int> vec_int;
    for (size_t id : vec_size_t)
    {
        vec_int.push_back(static_cast<int>(id));
    }
    return vec_int;
}

void lee::PandoraLEEAnalyzer::analyze(art::Event const &evt)
{
    clear();

    _run = evt.run();
    _subrun = evt.subRun();
    _event = evt.id().event();

    std::cout << "[PandoraLEE] "
              << "RUN " << _run << " SUBRUN " << _subrun << " EVENT " << _event
              << std::endl;

    std::vector<size_t> nu_candidates;

    _event_passed = int(fElectronEventSelectionAlg.eventSelected(evt));

    /*
    if (_event_passed && !evt.isRealData()) {
    art::Handle<std::vector<simb::GTruth>> gTruthHandle;
    evt.getByLabel("generator", gTruthHandle);
    if (!gTruthHandle.isValid())
      return;
    std::vector<art::Ptr<simb::GTruth>> gTruthVec;
    art::fill_ptr_vector(gTruthVec, gTruthHandle);
    if (gTruthVec.size() == 0)
    {
      std::cout << "\n[NUMUSEL] No GTruth Information" << std::endl;
      return;
    }

    art::Handle<std::vector<simb::MCFlux>> mcFluxHandle;
    evt.getByLabel("generator", mcFluxHandle);
    if (!mcFluxHandle.isValid())
      return;
    std::vector<art::Ptr<simb::MCFlux>> mcFluxVec;
    art::fill_ptr_vector(mcFluxVec, mcFluxHandle);
    if (mcFluxVec.size() == 0)
    {
      std::cout << "\n[NUMUSEL] No MCFlux Information" << std::endl;
      return;
    }

    art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
    evt.getByLabel("generator", mcTruthHandle);
    if (!mcTruthHandle.isValid())
      return;
    std::vector<art::Ptr<simb::MCTruth>> mcTruthVec;
    art::fill_ptr_vector(mcTruthVec, mcTruthHandle);
    if (mcTruthVec.size() == 0)
    {
      std::cout << "\n[NUMUSEL] No MCTruth Information" << std::endl;
      return;
    }

    const art::Ptr<simb::MCFlux> mcFlux = mcFluxVec.at(0);
    const art::Ptr<simb::GTruth> gTruth = gTruthVec.at(0);
    const art::Ptr<simb::MCTruth> mcTruth = mcTruthVec.at(0);

    ewutil.WriteTree(evt, mcFlux, mcTruth, gTruth);
    }
    */

    _flash_PE = fElectronEventSelectionAlg.get_flash_PE();
    _flash_time = fElectronEventSelectionAlg.get_flash_time();

    _flash_x = fElectronEventSelectionAlg.get_flash_x();
    _TPC_x = fElectronEventSelectionAlg.get_TPC_x();
    _category = 0;
    std::vector<double> true_neutrino_vertex(3);
    std::cout << "[PandoraLEEAnalyzer] Real data " << evt.isRealData() << std::endl;

    // Get info from Marco's analyzer
    art::Handle<std::vector<ubana::SelectionResult>> selection_h;
    evt.getByLabel("UBXSec", selection_h);

    if (!selection_h.isValid() || selection_h->empty())
    {
        std::cout << "[PandoraLEEAnalyzer] SelectionResult handle is not valid or empty." << std::endl;
    }

    std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
    art::fill_ptr_vector(selection_v, selection_h);

    if (selection_v.size() > 0)
    {
        _numu_passed = int(selection_v.at(0)->GetSelectionStatus());
        if (selection_v.at(0)->GetSelectionStatus())
        {
            std::cout << "[PandoraLEEAnalyzer] Event is selected by UBXSec" << std::endl;
        }
        else
        {
            std::cout << "[PandoraLEEAnalyzer] Event is not selected by UBXSec" << std::endl;
            std::cout << "[PandoraLEEAnalyzer] Failure reason " << selection_v.at(0)->GetFailureReason() << std::endl;
        }
        std::map<std::string, bool> failure_map = selection_v.at(0)->GetCutFlowStatus();
        for (auto iter : failure_map)
        {
            std::cout << "[PandoraLEEAnalyzer] UBXSec Cut: " << iter.first << "  >>>  " << (iter.second ? "PASSED" : "NOT PASSED") << std::endl;
            if (iter.second)
            {
                _numu_cuts += 1;
            }
        }
    }

    if ((!evt.isRealData() || m_isOverlaidSample))
    {
        // nu_e flux must be corrected by event weight

        try
        {
            art::InputTag eventweight_tag("eventweight");
            auto const &eventweights_handle =
                evt.getValidHandle<std::vector<evwgh::MCEventWeight>>(eventweight_tag);
            auto const &eventweights(*eventweights_handle);

            if (eventweights.size() > 0)
            {

                for (auto last : eventweights.at(0).fWeight)
                {
                    if (last.first.find("bnbcorrection") != std::string::npos)
                    {

                        if (!std::isfinite(last.second.at(0)))
                        {
                            _bnbweight = 1;
                        }
                        else
                        {
                            _bnbweight = last.second.at(0);
                        }
                    }
                }
            }
            else
            {
                _bnbweight = 1;
            }
        }
        catch (...)
        {
            std::cout << "[PandoraLEEAnalyzer] No MCEventWeight data product" << std::endl;
            _bnbweight = 1;
        }

        auto const &generator_handle =
            evt.getValidHandle<std::vector<simb::MCTruth>>(_mctruthLabel);
        auto const &generator(*generator_handle);
        _n_true_nu = generator.size();

        bool there_is_a_neutrino = false;
        //std::cout << "[PandoraLEEAnalyzer] Generator size " << generator.size() << std::endl;
        if (generator.size() > 0)
        {
            for (auto &gen : generator)
            {
                //std::cout << "[PandoraLEEAnalyzer] Generator origin " << gen.Origin() << std::endl;

                if (gen.Origin() == simb::kBeamNeutrino)
                {
                    there_is_a_neutrino = true;
                    _nu_pdg = gen.GetNeutrino().Nu().PdgCode();

                    _nu_energy = gen.GetNeutrino().Nu().E();
                    _ccnc = gen.GetNeutrino().CCNC();
                    _qsqr = gen.GetNeutrino().QSqr();
                    _theta = gen.GetNeutrino().Theta();
                    _lepton_E = gen.GetNeutrino().Lepton().E();
                    _lepton_theta = gen.GetNeutrino().Lepton().Momentum().Theta();

                    if (_ccnc == simb::kNC)
                    {
                        _category = k_nc;
                    }

                    true_neutrino_vertex[0] = gen.GetNeutrino().Nu().Vx();
                    true_neutrino_vertex[1] = gen.GetNeutrino().Nu().Vy();
                    true_neutrino_vertex[2] = gen.GetNeutrino().Nu().Vz();
                    _true_vx = true_neutrino_vertex[0];
                    _true_vy = true_neutrino_vertex[1];
                    _true_vz = true_neutrino_vertex[2];
                    _true_nu_is_fiducial = geoHelper.isFiducial(true_neutrino_vertex);

                    _interaction_type = gen.GetNeutrino().InteractionType();

                    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
                    if (SCE->GetPosOffsets(_true_vx, _true_vy, _true_vz).size() == 3)
                    {

                        _true_vx_sce =
                            _true_vx - SCE->GetPosOffsets(_true_vx, _true_vy, _true_vz)[0] + 0.7;
                        _true_vy_sce =
                            _true_vy + SCE->GetPosOffsets(_true_vx, _true_vy, _true_vz)[1];
                        _true_vz_sce =
                            _true_vz + SCE->GetPosOffsets(_true_vx, _true_vy, _true_vz)[2];
                    }
                    else
                    {
                        std::cout << "[PandoraLEEAnalyzer] "
                                  << "Space Charge service offset size < 3" << std::endl;
                        continue;
                    }

                    if (!geoHelper.isActive(true_neutrino_vertex))
                    {
                        _category = k_dirt;
                    }
                }
            }
        }

        if (!there_is_a_neutrino)
            _category = k_cosmic;

        auto const &mcparticles_handle = evt.getValidHandle<std::vector<simb::MCParticle>>("largeant");
        auto const &mcparticles(*mcparticles_handle);

        for (auto &mcparticle : mcparticles)
        {
            if (!(mcparticle.Process() == "primary" &&
                  mcparticle.T() != 0 &&
                  mcparticle.StatusCode() == 1))
                continue;

            const auto mc_truth = pandoraHelper.TrackIDToMCTruth(evt, "largeant", mcparticle.TrackId());
            if (mc_truth->Origin() == simb::kBeamNeutrino)
            {
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

        //Insert block to save the start point of the MCshower object for all showers that have a neutrino as mother and a kbeamneutrino as origin
        auto const &mcshower_handle = evt.getValidHandle<std::vector<sim::MCShower>>("mcreco");
        for (size_t _i_mcs = 0; _i_mcs < mcshower_handle->size(); _i_mcs++)
        {
            int pdg_mother = mcshower_handle->at(_i_mcs).MotherPdgCode();
            int origin = mcshower_handle->at(_i_mcs).Origin();

            if ((pdg_mother == 22 || pdg_mother == 11) && origin == 1)
            {
                _true_shower_pdg.push_back(mcshower_handle->at(_i_mcs).AncestorPdgCode());
                _true_shower_depE.push_back(mcshower_handle->at(_i_mcs).DetProfile().E());

                double x_det = mcshower_handle->at(_i_mcs).Start().X();
                double y_det = mcshower_handle->at(_i_mcs).Start().Y();
                double z_det = mcshower_handle->at(_i_mcs).Start().Z();

                if (pdg_mother == 22)
                { //For photons take the end of the shower
                    x_det = mcshower_handle->at(_i_mcs).End().X();
                    y_det = mcshower_handle->at(_i_mcs).End().Y();
                    z_det = mcshower_handle->at(_i_mcs).End().Z();
                }

                std::vector<double> dqdx = mcshower_handle->at(_i_mcs).dQdx();

                auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
                _true_shower_x_sce.push_back(x_det - SCE->GetPosOffsets(x_det, y_det, z_det)[0] + 0.7);
                _true_shower_y_sce.push_back(y_det + SCE->GetPosOffsets(x_det, y_det, z_det)[1]);
                _true_shower_z_sce.push_back(z_det + SCE->GetPosOffsets(x_det, y_det, z_det)[2]);
            }
        }

        if (_category != k_cosmic && _category != k_dirt && _category != k_nc)
        {
            if (abs(_nu_pdg) == 12)
            {
                _category = k_nu_e;
            }
            if (abs(_nu_pdg) == 14)
            {
                _category = k_nu_mu;
            }
        }
        _true_1eX_signal = featureHelper.true_thresholds_1eX(_true_nu_is_fiducial, _nu_daughters_pdg, _nu_daughters_E);
    }
    else
    {
        _category = k_data;
    }

    //std::cout << "[PandoraLEE] "
    //          << "True neutrino PDG " << _nu_pdg << std::endl;
    //std::cout << "[PandoraLEE] "
    //          << "Nu energy " << _nu_energy << std::endl;

    std::vector<art::Ptr<recob::PFParticle>> neutrino_pf;
    std::vector<art::Ptr<recob::PFParticle>> cosmic_pf;

    std::vector<int> neutrino_pdg;
    std::vector<std::string> neutrino_process;
    std::vector<double> neutrino_energy;

    std::vector<int> cosmic_pdg;
    std::vector<std::string> cosmic_process;
    std::vector<double> cosmic_energy;

    if ((!evt.isRealData() || m_isOverlaidSample))
    {

        categorizePFParticles(evt,
                              neutrino_pdg, neutrino_process, neutrino_energy, neutrino_pf,
                              cosmic_pdg, cosmic_process, cosmic_energy, cosmic_pf);

        _n_matched = neutrino_pf.size();
    }

    size_t ipf_candidate = std::numeric_limits<size_t>::lowest();
    if (_event_passed)
    {
        for (auto &inu : fElectronEventSelectionAlg.get_primary_indexes())
        {
            if (fElectronEventSelectionAlg.get_neutrino_candidate_passed().at(inu))
            {
                ipf_candidate = inu;
                break;
            }
        }

        std::cout << "[PandoraLEE] "
                  << "EVENT PASSED" << std::endl;
        auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);

        recob::PFParticle const &pfpneutrino = pfparticle_handle->at(ipf_candidate);

        _chosen_candidate = ipf_candidate;
        _candidate_pdg = pfpneutrino.PdgCode();

        art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
        auto const &vertex_obj = vertex_per_pfpart.at(ipf_candidate);

        double reco_neutrino_vertex[3];
        vertex_obj->XYZ(reco_neutrino_vertex);
        _vx = reco_neutrino_vertex[0];
        _vy = reco_neutrino_vertex[1];
        _vz = reco_neutrino_vertex[2];
        _fiducial = geoHelper.isFiducial(reco_neutrino_vertex);

        if (_category != k_data)
        {
            TVector3 v_reco_vertex(_vx, _vy, _vz);
            TVector3 sce_true_vertex(_true_vx_sce, _true_vy_sce, _true_vz_sce);

            _distance = geoHelper.distance(v_reco_vertex, sce_true_vertex);
        }

        // Load all the initializations of associations
        auto const &trackVecHandle = evt.getValidHandle<std::vector<recob::Track>>(m_pfp_producer);
        auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(m_pfp_producer);
        auto const &cluster_handle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pfp_producer);

        art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
        art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
        art::FindMany<anab::ParticleID> fmpid(trackVecHandle, evt, m_pid_producer);

        art::FindManyP<anab::CosmicTag> dtAssns(trackVecHandle, evt, "decisiontreeid");

        art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
        art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt, m_pfp_producer);

        art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
        art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, m_pfp_producer);

        _n_tracks = fElectronEventSelectionAlg.get_n_tracks().at(ipf_candidate);

        if (_n_tracks > 0)
        {
            _nu_track_ids = fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary().at(ipf_candidate);
            for (auto &pf_id : _nu_track_ids)
            {

                recob::PFParticle const &pfparticle = pfparticle_handle->at(pf_id);
                _nu_track_daughters.push_back(vectorCast(pfparticle.Daughters()));

                std::vector<double> dqdx(3, std::numeric_limits<double>::lowest());
                // Currently not used, reserved for calibrated variant.
                std::vector<float> dqdx_cali(3, std::numeric_limits<float>::lowest());
                std::vector<double> dedx(3, std::numeric_limits<double>::lowest());
                std::vector<double> dqdx_hits_track;
                std::vector<int> dqdx_wires_track;

                energyHelper.dQdx(pf_id, evt, dqdx, dqdx_cali, dqdx_hits_track, dqdx_wires_track, m_dQdxRectangleLength, m_dQdxRectangleWidth, m_pfp_producer);
                _track_dQdx_cali.push_back(dqdx_cali);
                _track_dQdx_hits.push_back(dqdx_hits_track);
                _track_dQdx_wires.push_back(dqdx_wires_track);

                std::vector<double> dedx_hits_track(dqdx_hits_track.size(), std::numeric_limits<double>::lowest());

                energyHelper.dEdxFromdQdx(dedx, dqdx);

                energyHelper.dEdxFromdQdx(dedx_hits_track, dqdx_hits_track);
                _track_dEdx_hits.push_back(dedx_hits_track);

                _track_dQdx.push_back(dqdx);
                _track_dEdx.push_back(dedx);

                art::Ptr<recob::Track> const &track_obj = track_per_pfpart.at(pf_id);

                if (track_obj.isNull())
                {
                    std::cout << "[LEE Analyzer] track is Null" << std::endl;
                    continue;
                }

                double mean = std::numeric_limits<double>::lowest();
                double stdev = std::numeric_limits<double>::lowest();

                energyHelper.trackResiduals(evt, m_pfp_producer, track_obj, mean, stdev);
                _track_res_mean.push_back(mean);
                _track_res_std.push_back(stdev);

                std::vector<const anab::ParticleID *> pids = fmpid.at(track_obj.key());

                for (size_t ipid = 0; ipid < pids.size(); ++ipid)
                {
                    if (!pids[ipid] || !pids[ipid]->PlaneID().isValid)
                    {
                        _track_pida.push_back(std::numeric_limits<double>::lowest());
                        _track_pidchipr.push_back(std::numeric_limits<double>::lowest());
                        _track_pidchika.push_back(std::numeric_limits<double>::lowest());
                        _track_pidchipi.push_back(std::numeric_limits<double>::lowest());
                        _track_pidchimu.push_back(std::numeric_limits<double>::lowest());
                        _track_pidchi.push_back(std::numeric_limits<double>::lowest());
                        break;
                    }

                    int planenum = pids[ipid]->PlaneID().Plane;
                    if (planenum != 2)
                        continue;

                    _track_pida.push_back(pids[ipid]->PIDA());
                    _track_pidchipr.push_back(pids[ipid]->Chi2Proton());
                    _track_pidchika.push_back(pids[ipid]->Chi2Kaon());
                    _track_pidchipi.push_back(pids[ipid]->Chi2Pion());
                    _track_pidchimu.push_back(pids[ipid]->Chi2Muon());
                    _track_pidchi.push_back(pids[ipid]->MinChi2());
                }
                _matched_tracks.push_back(std::numeric_limits<int>::lowest());
                _matched_tracks_process.push_back("");
                _matched_tracks_energy.push_back(std::numeric_limits<double>::lowest());

                std::vector<art::Ptr<anab::CosmicTag>> dtVec = dtAssns.at(track_obj.key());
                for (auto const &dttag : dtVec)
                {
                    if (dttag->CosmicType() == TAGID_P)
                    {
                        _predict_p.push_back(dttag->CosmicScore());
                    }
                    else if (dttag->CosmicType() == TAGID_MU)
                        _predict_mu.push_back(dttag->CosmicScore());
                    else if (dttag->CosmicType() == TAGID_PI)
                        _predict_pi.push_back(dttag->CosmicScore());
                    else if (dttag->CosmicType() == TAGID_EM)
                        _predict_em.push_back(dttag->CosmicScore());
                    else if (dttag->CosmicType() == TAGID_CS)
                        _predict_cos.push_back(dttag->CosmicScore());
                }

                std::vector<double> start_point;
                std::vector<double> end_point;

                start_point.push_back(track_obj->Start().X());
                start_point.push_back(track_obj->Start().Y());
                start_point.push_back(track_obj->Start().Z());

                end_point.push_back(track_obj->End().X());
                end_point.push_back(track_obj->End().Y());
                end_point.push_back(track_obj->End().Z());

                // needed for amount of hits and charge weightd calibration constant
                std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(pf_id);
                std::vector<int> nhits_spacepoint(3, 0);
                std::vector<float> cali_corr(3, 0);
                std::vector<float> total_charge(3, 0);

                for (auto &_sps : spcpnts)
                {
                    std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
                    const double *xyz = _sps->XYZ();

                    for (auto &hit : hits)
                    {
                        auto plane_nr = hit->View();
                        if (plane_nr > 2 || plane_nr < 0)
                            continue;
                        nhits_spacepoint[plane_nr]++;
                        total_charge[plane_nr] += hit->Integral();
                        float yzcorrection = energyCalibProvider.YZdqdxCorrection(plane_nr, xyz[1], xyz[2]);
                        float xcorrection = energyCalibProvider.XdqdxCorrection(plane_nr, xyz[0]);
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
                    //Else the value is zero
                }
                _track_energy_cali.push_back(cali_corr);
                _track_nhits_spacepoint.push_back(nhits_spacepoint);

                std::vector<double> this_energy;
                std::vector<int> this_nhits;

                energyHelper.energyFromHits(pfparticle, this_nhits, this_energy, evt, m_pfp_producer);

                _track_energy_hits.push_back(this_energy);
                // Alternative way to calculate the energy using dedx.
                _track_energy_dedx.push_back(energyHelper.trackEnergy_dedx(track_obj, evt, m_pfp_producer));

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

                std::vector<std::vector<double>> pca;
                pca.resize(3, std::vector<double>(2));

                energyHelper.PCA(pf_id, evt, pca, m_pfp_producer);

                double weighted_pca = 0;
                int total_hits = 0;
                for (size_t i = 0; i < this_nhits.size(); i++)
                {
                    weighted_pca += this_nhits[i] * pca[i][0];
                    total_hits += this_nhits[i];
                }

                weighted_pca /= total_hits;
                _track_pca.push_back(weighted_pca);
                _track_nhits_cluster.push_back(this_nhits);

                //std::cout << "[PCA] Track " << pca[2][0] << " " << pca[2][1] << std::endl;
            }
        }

        try
        {
            _nu_shower_ids = fElectronEventSelectionAlg.get_pfp_id_showers_from_primary().at(
                ipf_candidate);
            _n_showers = fElectronEventSelectionAlg.get_n_showers().at(ipf_candidate);
        }
        catch (...)
        {
            std::cout << "[PandoraLEE] Error getting associated showers to neutrino "
                         "candidate"
                      << std::endl;
        }

        for (auto &pf_id : _nu_shower_ids)
        {

            auto const &shower_obj = shower_per_pfpart.at(pf_id);
            if (shower_obj.isNull())
            {
                std::cout << "[PandoraLEEAnalyzer]  aranShower pointer " << pf_id << " is null, exiting" << std::endl;
                continue;
            }

            recob::PFParticle const &pfparticle = pfparticle_handle->at(pf_id);
            _nu_shower_daughters.push_back(vectorCast(pfparticle.Daughters()));

            std::vector<double> dqdx(3, std::numeric_limits<double>::lowest());
            // Currently not used, reserved for calibrated variant.
            std::vector<float> dqdx_cali(3, std::numeric_limits<float>::lowest());
            std::vector<double> dedx(3, std::numeric_limits<double>::lowest());
            std::vector<double> dqdx_hits_shower;
            std::vector<int> dqdx_wires_shower;

            _matched_showers.push_back(std::numeric_limits<int>::lowest());
            _matched_showers_process.push_back("");

            _matched_showers_energy.push_back(std::numeric_limits<double>::lowest());

            energyHelper.dQdx(pf_id, evt, dqdx, dqdx_cali, dqdx_hits_shower, dqdx_wires_shower, m_dQdxRectangleLength,
                              m_dQdxRectangleWidth, m_pfp_producer);

            _shower_dQdx_hits.push_back(dqdx_hits_shower);
            _shower_dQdx_cali.push_back(dqdx_cali);
            _shower_dQdx_wires.push_back(dqdx_wires_shower);
            //std::cout << "[dQdx] nohits " << dqdx_hits_shower.size() << " " << dqdx[0] << " " << dqdx[1] << " " << dqdx[2] << std::endl;

            std::vector<double> dedx_hits_shower(dqdx_hits_shower.size(), std::numeric_limits<double>::lowest());

            energyHelper.dEdxFromdQdx(dedx, dqdx);
            energyHelper.dEdxFromdQdx(dedx_hits_shower, dqdx_hits_shower);
            _shower_dEdx_hits.push_back(dedx_hits_shower);

            _shower_dQdx.push_back(dqdx);
            _shower_dEdx.push_back(dedx);

            _shower_dir_x.push_back(shower_obj->Direction().X());
            _shower_dir_y.push_back(shower_obj->Direction().Y());
            _shower_dir_z.push_back(shower_obj->Direction().Z());

            _shower_open_angle.push_back(shower_obj->OpenAngle());
            double shower_length = shower_obj->Length();
            _shower_length.push_back(shower_length);

            std::vector<double> start_point;
            std::vector<double> end_point;
            start_point.resize(3);
            end_point.resize(3);
            for (int ix = 0; ix < 3; ix++)
            {
                start_point[ix] = shower_obj->ShowerStart()[ix];
                end_point[ix] =
                    shower_obj->ShowerStart()[ix] + shower_length;
            }

            _shower_start_x.push_back(shower_obj->ShowerStart().X());
            _shower_start_y.push_back(shower_obj->ShowerStart().Y());
            _shower_start_z.push_back(shower_obj->ShowerStart().Z());

            _shower_energy_product.push_back(shower_obj->Energy());

            _shower_phi.push_back(shower_obj->Direction().Phi());
            _shower_theta.push_back(shower_obj->Direction().Theta());

            // needed for amount of hits and charge weightd calibration constant
            std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(pf_id);
            std::vector<int> nhits_spacepoint(3, 0);
            std::vector<float> cali_corr(3, 0);
            std::vector<float> total_charge(3, 0);

            for (auto &_sps : spcpnts)
            {
                std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
                const double *xyz = _sps->XYZ();

                for (auto &hit : hits)
                {
                    auto plane_nr = hit->View();
                    if (plane_nr > 2 || plane_nr < 0)
                        continue;
                    nhits_spacepoint[plane_nr]++;
                    total_charge[plane_nr] += hit->Integral();
                    float yzcorrection = energyCalibProvider.YZdqdxCorrection(plane_nr, xyz[1], xyz[2]);
                    float xcorrection = energyCalibProvider.XdqdxCorrection(plane_nr, xyz[0]);
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
                //Else the value is zero
            }
            _shower_energy_cali.push_back(cali_corr);
            _shower_nhits_spacepoint.push_back(nhits_spacepoint);

            std::vector<double> this_energy;
            std::vector<int> this_nhits;

            energyHelper.energyFromHits(pfparticle, this_nhits, this_energy, evt, m_pfp_producer);
            _shower_energy_hits.push_back(this_energy);

            std::vector<std::vector<double>> pca;
            pca.resize(3, std::vector<double>(2));

            energyHelper.PCA(pf_id, evt, pca, m_pfp_producer);

            double weighted_pca = 0;
            int total_hits = 0;
            for (size_t i = 0; i < this_nhits.size(); i++)
            {
                weighted_pca += this_nhits[i] * pca[i][0];
                total_hits += this_nhits[i];
            }

            weighted_pca /= total_hits;
            _shower_pca.push_back(weighted_pca);
            _shower_nhits_cluster.push_back(this_nhits);

            //std::cout << "[PCA] Shower " << pca[2][0] << " " << pca[2][1] << std::endl;

            for (auto &_sps : spcpnts)
            {
                std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
                auto xyz = _sps->XYZ();
                float integral = 0;

                for (auto &hit : hits)
                {
                    if (hit->View() == geo::kZ)
                    {
                        integral += hit->Integral();
                    }
                }
                // Only save points with nonzero collection charge
                if (integral > 0)
                {
                    _shower_sp_int.push_back(integral);
                    _shower_sp_x.push_back(xyz[0]);
                    _shower_sp_y.push_back(xyz[1]);
                    _shower_sp_z.push_back(xyz[2]);
                }
            }
        }

        /* For each neutrino PFParticle checks how many daughter tracks and showers
    with true neutrino origin we have */
        _nu_matched_showers = 0;
        bool shower_cr_found = false;

        _nu_matched_tracks = 0;
        bool track_cr_found = false;

        for (auto &inu : fElectronEventSelectionAlg.get_primary_indexes())
        {
            _primary_indexes.push_back(inu);

            _flash_passed.push_back(fElectronEventSelectionAlg.get_op_flash_indexes().at(inu));
            _number_tracks.push_back(fElectronEventSelectionAlg.get_n_tracks().at(inu));
            _number_showers.push_back(fElectronEventSelectionAlg.get_n_showers().at(inu));

            std::vector<size_t> pfp_showers_id = fElectronEventSelectionAlg.get_pfp_id_showers_from_primary().at(inu);

            int pass_shower = 0;

            for (size_t ish = 0; ish < pfp_showers_id.size(); ish++)
            {

                for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++)
                {

                    if (pfp_showers_id[ish] == neutrino_pf[ipf].key())
                    {
                        pass_shower += 1;

                        if (inu == ipf_candidate)
                        {
                            _nu_matched_showers++;
                            _matched_showers[ish] = neutrino_pdg[ipf];
                            _matched_showers_process[ish] = neutrino_process[ipf];
                            _matched_showers_energy[ish] = neutrino_energy[ipf];

                            // Also fill the daugthers with this info:
                        }
                    }
                }

                if (inu == ipf_candidate)
                {
                    for (size_t ipf = 0; ipf < cosmic_pf.size(); ipf++)
                    {
                        if (pfp_showers_id[ish] == cosmic_pf[ipf].key())
                        {
                            shower_cr_found = true;
                            _matched_showers[ish] = cosmic_pdg[ipf];
                            _matched_showers_process[ish] = cosmic_process[ipf];
                            _matched_showers_energy[ish] = cosmic_energy[ipf];
                        }
                    }
                }
            }

            _shower_passed.push_back(pass_shower);

            if (_number_tracks.back() > 0)
            {
                std::vector<size_t> pfp_tracks_id = fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary().at(inu);

                int pass_track = 0;

                for (size_t itr = 0; itr < pfp_tracks_id.size(); itr++)
                {

                    for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++)
                    {
                        if (pfp_tracks_id[itr] == neutrino_pf[ipf].key())
                        {
                            pass_track += 1;
                            if (inu == ipf_candidate)
                            {
                                _nu_matched_tracks++;
                                _matched_tracks[itr] = neutrino_pdg[ipf];
                                _matched_tracks_process[itr] = neutrino_process[ipf];
                                _matched_tracks_energy[itr] = neutrino_energy[ipf];
                            }
                        }
                    }

                    if (inu == ipf_candidate)
                    {
                        for (size_t ipf = 0; ipf < cosmic_pf.size(); ipf++)
                        {
                            if (pfp_tracks_id[itr] == cosmic_pf[ipf].key())
                            {
                                track_cr_found = true;
                                _matched_tracks[itr] = cosmic_pdg[ipf];
                                _matched_tracks_process[itr] = cosmic_process[ipf];
                                _matched_tracks_energy[itr] = cosmic_energy[ipf];
                            }
                        }
                    }
                }

                _track_passed.push_back(pass_track);
            }
        }

        if (!track_cr_found && _nu_matched_tracks == 0 && !shower_cr_found && _nu_matched_showers == 0 && _category != k_dirt && _category != k_data)
        {
            _category = k_other;
            std::cout << "[PandoraLEE] "
                      << "***NOT NEUTRINO NOR COSMIC***" << std::endl;
        }

        // if ((track_cr_found && _nu_matched_tracks > 0)
        // || (shower_cr_found && _nu_matched_showers > 0)) {
        //   _category = k_mixed;
        //   std::cout << "[PandoraLEE] "
        //             << "***MIXED COSMIC/NEUTRINO***" << std::endl;
        // }

        if ((track_cr_found || shower_cr_found) &&
            (_nu_matched_tracks > 0 || _nu_matched_showers > 0))
        {
            _category = k_mixed;
        }

        if ((track_cr_found || shower_cr_found) && _category != k_mixed)
        {
            _category = k_cosmic;
        }

        std::cout << "[PandoraLEE] "
                  << "Category " << _category << std::endl;

        //Fill the FeatureHelper fields:
        featureHelper.reco_maxangle(_shower_dir_x, _shower_dir_y, _shower_dir_z,
                                    _track_dir_x, _track_dir_y, _track_dir_z,
                                    _track_maxangle, _shower_maxangle);

        featureHelper.reco_hierarchy(_nu_track_ids, evt, m_pfp_producer,
                                     _track_daughter, _track_is_daughter);

        featureHelper.reco_hierarchy(_nu_shower_ids, evt, m_pfp_producer,
                                     _shower_daughter, _shower_is_daughter);

        _shower_vtxdistance = featureHelper.reco_vtxdistance(_vx, _vy, _vz, _shower_start_x, _shower_start_y, _shower_start_z);
        _track_vtxdistance = featureHelper.reco_vtxdistance(_vx, _vy, _vz, _track_start_x, _track_start_y, _track_start_z);

        featureHelper.reco_spacepoint_ratios(_nu_shower_ids, evt, m_pfp_producer,
                                             _shower_fidvol_ratio, _shower_spacepoint_dqdx_ratio);

        std::vector<float> _track_fidvol_ratio; // Dummy variable, not saving this field.
        featureHelper.reco_spacepoint_ratios(_nu_track_ids, evt, m_pfp_producer,
                                             _track_fidvol_ratio, _track_spacepoint_dqdx_ratio);

        featureHelper.reco_flash_info(_flash_passed, _flash_PE, _flash_time,
                                      _flash_PE_max, _flash_time_max);

        _track_bdt_precut = featureHelper.reco_bdt_track_precut(_predict_mu, _predict_cos, _n_tracks);

        featureHelper.reco_track_containment(_track_end_x, _track_end_y, _track_end_z, _track_containment);

        featureHelper.reco_totalChargeCenter(_nu_shower_ids, _nu_track_ids, evt, m_pfp_producer,
                                             _chargecenter_x, _chargecenter_y, _chargecenter_z);

        _total_spacepoint_containment = featureHelper.reco_total_spacepoint_containment(_nu_shower_ids, _nu_track_ids, evt, m_pfp_producer);

        featureHelper.reco_dedx(_shower_dEdx_hits,
                                _shower_dEdx,
                                _shower_dQdx_cali,
                                _shower_dedx_hits_w,
                                _shower_dedx_w,
                                _shower_dedx_best_w);

        featureHelper.reco_dedx(_track_dEdx_hits,
                                _track_dEdx,
                                _track_dQdx_cali,
                                _track_dedx_hits_w,
                                _track_dedx_w,
                                _track_dedx_best_w);

        featureHelper.reco_energy(_shower_energy_hits,
                                  _shower_energy_cali,
                                  _shower_nhits_cluster,
                                  _shower_nhits_spacepoint,
                                  _shower_energy_w,
                                  _shower_hitsratio_w,
                                  _shower_hits_w);

        featureHelper.reco_energy(_track_energy_hits,
                                  _track_energy_cali,
                                  _track_nhits_cluster,
                                  _track_nhits_spacepoint,
                                  _track_energy_w,
                                  _track_hitsratio_w,
                                  _track_hits_w);

        // Should work on overlaidsample
        if (!evt.isRealData())
        {
            featureHelper.true_closest_electron_matched(_matched_showers, _matched_tracks,
                                                        _true_vx_sce, _true_vy_sce, _true_vz_sce,
                                                        _shower_start_x, _shower_start_y, _shower_start_z,
                                                        _track_start_x, _track_start_y, _track_start_z,
                                                        _shower_cle, _track_cle);

            featureHelper.true_match_daughters(evt, m_pfp_producer,
                                               _nu_shower_ids, _nu_track_ids,
                                               _matched_showers, _matched_tracks);
        }
    }
    else
    {
        std::cout << "[PandoraLEE] "
                  << "EVENT NOT PASSED" << std::endl;

        for (auto &inu : fElectronEventSelectionAlg.get_primary_indexes())
        {
            _primary_indexes.push_back(inu);

            _flash_passed.push_back(fElectronEventSelectionAlg.get_op_flash_indexes().at(inu));
            _number_tracks.push_back(fElectronEventSelectionAlg.get_n_tracks().at(inu));
            _number_showers.push_back(fElectronEventSelectionAlg.get_n_showers().at(inu));

            std::vector<size_t> pfp_showers_id = fElectronEventSelectionAlg.get_pfp_id_showers_from_primary().at(inu);

            int pass_shower = 0;

            for (size_t ish = 0; ish < pfp_showers_id.size(); ish++)
            {
                for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++)
                {
                    if (pfp_showers_id[ish] == neutrino_pf[ipf].key())
                    {
                        pass_shower += 1;
                    }
                }
            }

            _shower_passed.push_back(pass_shower);

            if (_number_tracks.back() > 0)
            {
                std::vector<size_t> pfp_tracks_id = fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary().at(inu);
                int pass_track = 0;

                for (size_t itr = 0; itr < pfp_tracks_id.size(); itr++)
                {

                    for (size_t ipf = 0; ipf < neutrino_pf.size(); ipf++)
                    {
                        if (pfp_tracks_id[itr] == neutrino_pf[ipf].key())
                        {
                            pass_track += 1;
                        }
                    }
                }

                _track_passed.push_back(pass_track);
            }
        }

        // Fill fields that do not depend on passed or not passed.
        _n_primaries = _primary_indexes.size();
    }

    myTTree->Fill();
    std::cout << "[PandoraLEE] "
              << "END ANALYZER" << std::endl;

} // end analyze function

//------------------------------------------------------------------------------------------------------------------------------------

void lee::PandoraLEEAnalyzer::reconfigure(fhicl::ParameterSet const &pset)
{
    // UGLY AND DANGEROUS!
    geoHelper.setFiducialVolumeCuts(10, 10, 20, 20, 10, 50);

    // add what you want to read, and default values of your labels etc. example:
    fElectronEventSelectionAlg.reconfigure(pset.get<fhicl::ParameterSet>("ElectronSelectionAlg"));

    m_hitfinderLabel = pset.get<std::string>("HitFinderLabel", "pandoraCosmicHitRemoval::McRecoStage2");
    m_pid_producer = pset.get<std::string>("ParticleIDModuleLabel", "pandoraNupid::McRecoStage2");
    m_pfp_producer = pset.get<std::string>("PFParticleLabel", "pandoraNu::McRecoStage2");
    m_spacepointLabel = pset.get<std::string>("SpacePointLabel", "pandoraNu::McRecoStage2");
    //m_spacepointLabel = pset.get<std::string>("SpacePointLabel", "pandoraNu::PandoraLEEAnalyzer");

    m_printDebug = pset.get<bool>("PrintDebug", false);

    m_isData = pset.get<bool>("isData", false);
    m_isCosmicInTime = pset.get<bool>("isCosmicInTime", false);
    m_isOverlaidSample = pset.get<bool>("isOverlaidSample", false);

    m_dQdxRectangleWidth = pset.get<double>("dQdxRectangleWidth", 1);
    m_dQdxRectangleLength = pset.get<double>("dQdxRectangleLength", 4);
}

//---------------------------------------------------------------------------------------------------------------------------
// add other functions here

DEFINE_ART_MODULE(lee::PandoraLEEAnalyzer)
