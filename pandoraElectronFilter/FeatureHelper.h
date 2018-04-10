////////////////////////////////////////////////////////////////////////
// Class:       FeatureHelper
// File:        FeatureHelper.h
//
////////////////////////////////////////////////////////////////////////

#ifndef FEATUREHELPER_H
#define FEATUREHELPER_H

#include "HelperBase.h"

#include "GeometryHelper.h"

#include "lardataobj/RecoBase/Hit.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace lee
{

class FeatureHelper : public HelperBase
{
public:
  explicit FeatureHelper();
  ~FeatureHelper() = default;

  /**
   * @brief      Returns if the true event is within the fiducial volume and has at least 1 electron of min_e_kinE kinetic energy. 
   * @return     bool
   */
  int true_thresholds_1eX(int &_true_nu_is_fiducial,
                          std::vector<int> &_nu_daughters_pdg, std::vector<double> &_nu_daughters_E);

  /**
   * @brief      Returns two arrays, maxangle array track and maxangle array shower. 
   * @return     void, needs two empty arrays.
   */
  void reco_maxangle(std::vector<double> &_shower_dir_x, std::vector<double> &_shower_dir_y, std::vector<double> &_shower_dir_z,
                     std::vector<double> &_track_dir_x, std::vector<double> &_track_dir_y, std::vector<double> &_track_dir_z,
                     std::vector<double> &_track_maxangle, std::vector<double> &_shower_maxangle);

  /**
   * @brief      Returns two arrays. 
   *             Add info about the hierarchy, 0 = no daughter, 1 = has showerdaughter, 2 = has trackdaughter, 3 = more than 1 daughter
   *             0 means direct neutrino daughter, 1 means shower daughter, 2 means track daughter
   * @return     void, needs two empty arrays.
   */
  void reco_hierarchy(std::vector<size_t> &_nu_object_ids, const art::Event &evt, std::string _pfp_producer,
                      std::vector<int> &object_daughter, std::vector<int> &object_is_daughter);

  /**
   * @brief      Returns the total charge inside and outside the fiducial volume of the spacepoint array that it gets. 
   * @return     float, float
   */
  void reco_spacepoint_containment(size_t pf_id, const art::Event &evt, std::string _pfp_producer,
                                   float &in_fidvol_q, float &out_fidvol_q);

  /**
   * @brief      Returns two arrays, one carring the projected distances of the spacepoints, 
   *             the other the collection plane integrals. 
   * @return     float array, actual return value of function is the maximum distance.
   */
  float reco_spacepoint_dqdx(size_t object_id, const art::Event &evt, std::string _pfp_producer,
                             std::vector<float> &sps_dqdx_distance, std::vector<float> &sps_dqdx_integral);

  /**
   * @brief      Returns two arrays, the containment ratios, 
   *             the other the dqdx from sps ratios. Both for a list of ids 
   * @return     double,double
   */
  void reco_spacepoint_ratios(std::vector<size_t> _nu_object_ids, const art::Event &evt, std::string _pfp_producer,
                              std::vector<float> &object_containment, std::vector<float> &object_dqdx_ratio);

  /**
   * @brief      Returns the PE and time of the flash used for matching. 
   * @return     double,double
   */
  void reco_flash_info(std::vector<int> &_flash_passed, std::vector<double> &_flash_PE, std::vector<double> &_flash_time,
                       double &_flash_PE_max, double &_flash_time_max);

  /**
   * @brief      Fixes the matching for the daughter particles. 
   * @return     _matched_tracks,_matched_showers
   */
  void reco_match_daughters(const art::Event &evt, std::string _pfp_producer,
                            std::vector<size_t> &_nu_shower_ids, std::vector<size_t> &_nu_track_ids,
                            std::vector<int> &_matched_showers, std::vector<int> &_matched_tracks);

private:
  double mass_e = 0.00511;           // electron mass in GeV
  double min_e_kinE = 0.02 + mass_e; // minimum kinetic energy threshold for electrons
                                     //art::ServiceHandle<geo::Geometry> geo;
                                     //detinfo::DetectorProperties const *detprop;
  GeometryHelper geoHelper;
};
} // namespace lee

#endif