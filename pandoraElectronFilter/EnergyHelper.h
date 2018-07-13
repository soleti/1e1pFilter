////////////////////////////////////////////////////////////////////////
// Class:       EnergyHelper
// Module Type: filter
// File:        EnergyHelper.h
//
////////////////////////////////////////////////////////////////////////

#ifndef ENERGYHELPER_H
#define ENERGYHELPER_H

#include "HelperBase.h"

#include "GeometryHelper.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "TPrincipal.h"
#include "uboone/UBXSec/Algorithms/TrackQuality.h"
#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

namespace lee {

class EnergyHelper : public HelperBase {
public:
  explicit EnergyHelper(bool is_data);
  ~EnergyHelper() = default;

  /**
   * @brief      Measure the energy of a pfp_particle
   *
   * @param[in]  shower  The shower
   * @param[in]  evt     The art::event
   *
   * @return     Energy in units of GeV
   */
  void dQdx(size_t pfp_id,
            const art::Event &evt,
            std::vector<double> &dqdx,
            std::vector<double> &dqdx_cali,
            std::vector<double> &dqdx_hits,
            std::vector<double> &pitches,
            double m_dQdxRectangleLength,
            double m_dQdxRectangleWidth,
            std::string _pfp_producer);


  double PID(const std::vector<anab::ParticleID> *pids,
             int trackID,
             std::string AlgName,
             anab::kVariableType VariableType,
             int pdgCode);

  void dEdx_from_dQdx(std::vector<double> &dedx,
                    std::vector<double> &dqdx);

  void PCA(std::vector<art::Ptr<recob::Cluster>> *clusters,
           art::FindManyP<recob::Hit> *hits_per_clusters,
           std::vector<std::vector<double>> &pca_planes);

  void get_cali(std::vector<art::Ptr<recob::SpacePoint>> *spcpnts,
                art::FindManyP<recob::Hit> *hits_per_spcpnts,
                std::vector<double> &cali_corr);

  void  energy_from_hits(std::vector<art::Ptr<recob::Cluster>> *clusters,
                                    art::FindManyP<recob::Hit> *hits_per_cluster,
                                    std::vector<int>    &nHits,
                                    std::vector<double> &pfenergy);

  void cluster_residuals(std::vector<art::Ptr<recob::Cluster>> *clusters,
                         art::FindManyP<recob::Hit> *hits_per_cluster,
                         double &mean_v,
                         double &std_v);

  void track_dQdx(std::vector<art::Ptr<anab::Calorimetry>> *calos,
                  std::vector<double> &dqdx,
                  std::vector<double> &dedx);

  private:
    std::vector<double> _data_gain = {236.41, 228.83, 242.72}; // DocDB 14754
    std::vector<double> _mc_gain = {193.05, 196.85, 196.85};   // Plane 0, plane 1, plane 2
    bool _is_data = false;
    const lariov::TPCEnergyCalibProvider &_energy_calib_provider = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
    const detinfo::DetectorProperties *_detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    const double _drift = _detprop->DriftVelocity() * 1e-3;
    const double _readout_window = 4.8;
    const double _from_tick_to_ns = _readout_window / _detprop->ReadOutWindowSize() * 1e6;
    const double _wire_spacing = 0.3;
    const double _work_function = 23 / 1e6;
    const double _recombination_factor = 0.62;

    GeometryHelper geo_helper;
};
} // namespace lee

#endif
