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

#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibService.h"

#include "uboone/UBXSec/Algorithms/TrackQuality.h"

namespace lee
{

class EnergyHelper : public HelperBase
{
public:
  explicit EnergyHelper();
  ~EnergyHelper() = default;

  /**
   * @brief      Measure the energy of a track using dedx
   *
   * @param[in]  track  The track
   * @param[in]  evt    The art::event
   *
   * @return     Energy in units of GeV
   */

  double trackEnergy_dedx(const art::Ptr<recob::Track> &track,
                          const art::Event &evt,
                          std::string _pfp_producer);

  /**
   * @brief      Measure the energy of a pfp_particle
   *
   * @param[in]  shower  The shower
   * @param[in]  evt     The art::event
   *
   * @return     Energy in units of GeV
   */
  void energyFromHits(recob::PFParticle const &pfparticle,
                      std::vector<int> &nHits,
                      std::vector<double> &pfenergy,
                      const art::Event &evt,
                      std::string _pfp_producer);

  void dQdx(size_t pfp_id, const art::Event &evt,
            std::vector<double> &dqdx,
            std::vector<float> &dqdx_cali,
            std::vector<double> &dqdx_hits,
            std::vector<int> &dqdx_wires,
            double m_dQdxRectangleLength, double m_dQdxRectangleWidth,
            std::string _pfp_producer);

  void dEdxFromdQdx(std::vector<double> &dedx,
                    std::vector<double> &dqdx);

  void PCA(size_t pfp_id,
           const art::Event &evt,
           std::vector<std::vector<double>> &pca_planes,
           std::string _pfp_producer);

  void nHits(size_t pfp_id,
             const art::Event &evt,
             std::vector<int> &nHits,
             std::string _pfp_producer);

  void trackResiduals(const art::Event &e,
                      std::string _pfp_producer,
                      art::Ptr<recob::Track> candidate_track,
                      double &mean,
                      double &std);

private:
  std::vector<double> _data_gain = {239.5, 239.5, 239.5}; // Only measured of collection plane, David Caratelli
  std::vector<double> _mc_gain = {193.0, 197.0, 197.0};   // Plane 0, plane 1, plane 2

  // 23 work function (23 eV/e- in the argon)
  // 0.62 recombination factor
  double work_function = 23 / 1e6;
  double recombination_factor = 0.62;

  art::ServiceHandle<geo::Geometry> geo;
  detinfo::DetectorProperties const *detprop;  GeometryHelper geoHelper;

  //handle to tpc energy calibration provider
  const lariov::TPCEnergyCalibProvider& energyCalibProvider
    = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
};
} // namespace lee

#endif
