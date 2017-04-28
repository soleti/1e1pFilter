////////////////////////////////////////////////////////////////////////
// Class:       ElectronSelectionAna
// Plugin Type: analyzer (art v2_06_03)
// File:        ElectronSelectionAna_module.cc
//
// Generated at Mon Apr 24 16:33:33 2017 by Corey Adams using cetskelgen
// from cetlib version v2_03_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ElectronEventSelectionAlg.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "uboone/SpaceChargeServices/SpaceChargeServiceMicroBooNE.h"
#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "TTree.h"
#include "TVector3.h"

#include <algorithm>

namespace lee {
  class ElectronSelectionAna;
}
class lee::ElectronSelectionAna : public art::EDAnalyzer {
public:
  explicit ElectronSelectionAna(fhicl::ParameterSet const & pset);
  //virtual ~ElectronSelectionAna();

  void reconfigure(fhicl::ParameterSet const & pset);
  // Plugins should not be copied or assigned.
  ElectronSelectionAna(ElectronSelectionAna const &) = delete;
  ElectronSelectionAna(ElectronSelectionAna &&) = delete;
  ElectronSelectionAna & operator = (ElectronSelectionAna const &) = delete;
  ElectronSelectionAna & operator = (ElectronSelectionAna &&) = delete;
  // Required functions.
  void analyze(art::Event const & e) override;
private:

  // functions
  void fillTree(art::Event const & e);
  void fillTruthTree(art::Event const & e);
  void fillPandoraTree(art::Event const & e);
  void fillOticalTree(art::Event const & e);

  Float_t trackEnergy(const art::Ptr<recob::Track>& track, const art::Event & e);


  // variables
  lee::ElectronEventSelectionAlg fElectronEventSelectionAlg;
  //art::ServiceHandle<spacecharge::SpaceChargeServiceMicroBooNE> SCE;
  TTree*      fTree;

  bool bool_truth   = false;
  bool bool_pandora = true;
  bool bool_optical = false;

  //Run Subrun Event
  Short_t    run;
  Short_t    subrun;
  Int_t    event;

  //Truth information
  Short_t   mcevts_truth;                               ///< number of neutrino Interactions in the spill
  std::vector<Short_t>   nuPDG_truth;                   ///< neutrino PDG code
  std::vector<Short_t>   ccnc_truth;                    ///< 0=CC 1=NC
  std::vector<Short_t>   mode_truth;                    ///< 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  std::vector<Float_t>   enu_truth;                     ///< true neutrino energy
  std::vector<Float_t>   nuvtxx_truth;                  ///< neutrino vertex x
  std::vector<Float_t>   nuvtxy_truth;                  ///< neutrino vertex y
  std::vector<Float_t>   nuvtxz_truth;                  ///< neutrino vertex z

  //PandoraNu information
  Short_t              passed;                          ///< Return true if there is at least one neutrino vertex that passed
  Short_t              nnuvtx;                          ///< Number of PandoraNu neutrino candidate vertices
  std::vector<Float_t> nuvtxx;                          ///< x coordinate
  std::vector<Float_t> nuvtxy;                          ///< y coordinate
  std::vector<Float_t> nuvtxz;                          ///< z coordinate
  std::vector<Short_t> nuvtxpdg;                        ///< PDG code assigned by PandoraNu
  std::vector<Float_t> center_of_charge_x;              ///< x Center of deposited charge
  std::vector<Float_t> center_of_charge_y;              ///< y Center of deposited charge
  std::vector<Float_t> center_of_charge_z;              ///< z Center of deposited charge

  std::vector<std::vector< TVector3 >> shwr_dir;        ///< The direction of shower for every shower connected to the passed pandoraNu candidates
  std::vector<std::vector< Float_t >> shwr_en;          ///< Shower energy
  std::vector<std::vector< Float_t >> shwr_angle;       ///< Shower opening angle

  std::vector<std::vector< TVector3 >> trck_dir;        ///< The start direction of track for every track connected to the passed pandoraNu candidates
  std::vector<std::vector< Float_t >> trck_len;         ///< Length of the track
  std::vector<std::vector< Float_t >> trck_dedxavg;     ///< Average dedx of the track


  //Optical information
  Short_t nfls;                                         ///< Number of reconstructed flashes
  std::vector<Float_t> flsTime;                         ///< Flash time (us)
  std::vector<Float_t> flsPe;                           ///< Flash total PE
  std::vector<Float_t> flsYcenter;                      ///< Flash Y center (cm)
  std::vector<Float_t> flsZcenter;                      ///< Flash Z center (cm)
  std::vector<Float_t> flsYwidth;                       ///< Flash Y width (cm)
  std::vector<Float_t> flsZwidth;                       ///< Flash Z width (cm)
};
void lee::ElectronSelectionAna::reconfigure(fhicl::ParameterSet const & pset){
  fElectronEventSelectionAlg.reconfigure(pset.get<fhicl::ParameterSet>("ElectronSelectionAlg"));
}
lee::ElectronSelectionAna::ElectronSelectionAna(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)  // ,
 // More initializers here.
{
  //initialize output tree
  std::cout << "Initializing output tree..." << std::endl;
  art::ServiceHandle<art::TFileService> tfs;
  fTree  = tfs->make<TTree>("flashtree","FlashAnalysis Tree");

  //Set branches for (Run Subrun Event)
  fTree->Branch("run",     &run,     "run/S"       );
  fTree->Branch("subrun",  &subrun,  "subrun/S"    );
  fTree->Branch("event",   &event,   "event/I"     );

  //Set branches for truth information
  fTree->Branch("mcevts_truth", &mcevts_truth,             "mcevts_truth/S"        );
  fTree->Branch("nuPDG_truth",  "std::vector<Short_t>",    &nuPDG_truth            );
  fTree->Branch("ccnc_truth",   "std::vector<Short_t>",    &ccnc_truth             );
  fTree->Branch("mode_truth",   "std::vector<Short_t>",    &mode_truth             );
  fTree->Branch("enu_truth",    "std::vector<Float_t>",    &enu_truth              );
  fTree->Branch("nuvtxx_truth", "std::vector<Float_t>",    &nuvtxx_truth           );
  fTree->Branch("nuvtxy_truth", "std::vector<Float_t>",    &nuvtxy_truth           );
  fTree->Branch("nuvtxz_truth", "std::vector<Float_t>",    &nuvtxz_truth           );

  //Set branches for PandoraNU information
  fTree->Branch("passed",                   &passed,                   "passed/S"             );
  fTree->Branch("nnuvtx",                   &nnuvtx,                   "nnuvtx/S"             );
  fTree->Branch("nuvtxx",                   "std::vector<Float_t>",    &nuvtxx                );
  fTree->Branch("nuvtxy",                   "std::vector<Float_t>",    &nuvtxy                );
  fTree->Branch("nuvtxz",                   "std::vector<Float_t>",    &nuvtxz                );
  fTree->Branch("nuvtxpdg",                 "std::vector<Short_t>",    &nuvtxpdg              );
  fTree->Branch("center_of_charge_x",       "std::vector<Float_t>",    &center_of_charge_x    );
  fTree->Branch("center_of_charge_y",       "std::vector<Float_t>",    &center_of_charge_y    );
  fTree->Branch("center_of_charge_z",       "std::vector<Float_t>",    &center_of_charge_z    );

  fTree->Branch("shwr_dir",   "std::vector<std::vector<TVector3>>",    &shwr_dir              );  
  fTree->Branch("shwr_en",    "std::vector<std::vector<Float_t>>",     &shwr_en               );  
  fTree->Branch("shwr_angle", "std::vector<std::vector<Float_t>>",     &shwr_angle            );  

  fTree->Branch("trck_dir",   "std::vector<std::vector<TVector3>>",    &trck_dir              ); 
  fTree->Branch("trck_len",   "std::vector<std::vector<Float_t>>",     &shwr_en               );  
  fTree->Branch("trck_dedxavg","std::vector<std::vector<Float_t>>",    &shwr_angle            );  


  //Set branches for optical information
  fTree->Branch("nfls",       &nfls,       "nfls/S"                 );
  fTree->Branch("flsTime",    "std::vector<Float_t>",   &flsTime    );
  fTree->Branch("flsPe",      "std::vector<Float_t>",   &flsPe      );
  fTree->Branch("flsYcenter", "std::vector<Float_t>",   &flsYcenter );
  fTree->Branch("flsZcenter", "std::vector<Float_t>",   &flsZcenter );
  fTree->Branch("flsYwidth",  "std::vector<Float_t>",   &flsYwidth  );
  fTree->Branch("flsZwidth",  "std::vector<Float_t>",   &flsZwidth  );

  this->reconfigure(pset);
}


void lee::ElectronSelectionAna::analyze(art::Event const & e)
{
  try {
  fillTree(e);
  } catch(...) {std::cerr<<"Something went wrong filling root tree"<<std::endl;}
  return;
}


void lee::ElectronSelectionAna::fillTree(art::Event const & e)
{
  // Fill run information
  run    = e.run();
  subrun = e.subRun();
  event  = e.event();

  std::cout<<"\n Begin filling variables of (run,subrun,event) \t ("<< run <<"," <<subrun <<"," <<event<< ")"<< std::endl;
  // Fill truth information
  if(!e.isRealData() && bool_truth){fillTruthTree(e);}
  // Fill PandoraNu information
  if(bool_pandora){fillPandoraTree(e);}
  // Fill optical information
  if(bool_optical){fillOticalTree(e);}

  std::cout<<"variables filled, fill tree"<<std::endl;
  fTree->Fill();
}

Float_t lee::ElectronSelectionAna::trackEnergy(const art::Ptr<recob::Track>& track, const art::Event & evt)
{
  art::InputTag pandoraNu_tag { "pandoraNu" };
  auto const& track_handle = evt.getValidHandle< std::vector< recob::Track > >( pandoraNu_tag );
  art::FindManyP<anab::Calorimetry> calo_track_ass(track_handle, evt, "pandoraNucalo");
  const std::vector<art::Ptr<anab::Calorimetry>> calos = calo_track_ass.at(track.key());

  for (size_t ical = 0; ical < calos.size(); ++ical)
  {
    if (!calos[ical]) continue;
    if (!calos[ical]->PlaneID().isValid) continue;
    int planenum = calos[ical]->PlaneID().Plane;
    if (planenum < 0 || planenum > 2) continue;
    if (planenum != 2) continue;                           // Use informartion from collection plane only

    Float_t dedxsum=0;
    Float_t dedx=0;
    Float_t counter=0;
    for (size_t iTrkHit = 0; iTrkHit < calos[ical]->dEdx().size(); ++iTrkHit)
    {
      dedx = calos[ical]->dEdx()[iTrkHit];
      if (dedx > 0 && dedx < 10)
      {
        dedxsum+=dedx;
        counter++;
      }
    }
    return dedxsum/counter;
  }
  return 0.0;
}


void lee::ElectronSelectionAna::fillTruthTree(art::Event const & e)
{
  std::cout << "Filling truth information " << std::endl;
  art::InputTag truth_tag { "generator" };
  auto const& truth_handle = e.getValidHandle< std::vector< simb::MCTruth > >( truth_tag );

  mcevts_truth=0;
  nuPDG_truth.clear();
  ccnc_truth.clear();
  mode_truth.clear();
  enu_truth.clear();
  nuvtxx_truth.clear();
  nuvtxy_truth.clear();
  nuvtxz_truth.clear();

  if (truth_handle->size() > 0) {
    for(unsigned int iList = 0; iList < truth_handle->size() ; ++iList){
      if (truth_handle->at(iList).NeutrinoSet())
      {
        simb::MCNeutrino const& neutrino = truth_handle->at(iList).GetNeutrino();
        mcevts_truth++;
        nuPDG_truth.emplace_back(neutrino.Nu().PdgCode());
        ccnc_truth.emplace_back(neutrino.CCNC());
        mode_truth.emplace_back(neutrino.Mode());
        enu_truth.emplace_back(neutrino.Nu().E());

        nuvtxx_truth.emplace_back(neutrino.Nu().Vx());
        nuvtxy_truth.emplace_back(neutrino.Nu().Vy());
        nuvtxz_truth.emplace_back(neutrino.Nu().Vz());

        geo::Point_t point;
        point.SetXYZ(neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz());

        std::cout << "True vertex \t (" << point.X() <<"," <<point.Y() <<"," <<point.Z()<< ")"<< std::endl;

        //auto const* sc = SCE->provider();
        //geo::Vector_t SCcortrue = sc->GetPosOffsets(point);
        //std::cout << "SCcor vertex \t (" << SCcortrue.X() <<"," <<SCcortrue.Y() <<"," <<SCcortrue.Z()<< ")"<< std::endl;
      }
    }
  }
}

void lee::ElectronSelectionAna::fillPandoraTree(art::Event const & e)
{
  std::cout << "Filling PandoraNu information " << std::endl;
  passed = fElectronEventSelectionAlg.eventSelected(e);

  nnuvtx=0;
  nuvtxx.clear();
  nuvtxy.clear();
  nuvtxz.clear();
  center_of_charge_x.clear();
  center_of_charge_y.clear();
  center_of_charge_z.clear();
  nuvtxpdg.clear();

  shwr_dir.clear();
  shwr_en.clear();
  shwr_angle.clear();

  trck_dir.clear();
  trck_dedxavg.clear();
  trck_len.clear();

  if(passed)
  {
    art::InputTag pandoraNu_tag { "pandoraNu" };
    auto const& pfparticle_handle = e.getValidHandle< std::vector< recob::PFParticle > >( pandoraNu_tag );
    art::FindOneP< recob::Shower > shower_per_pfpart(pfparticle_handle, e, pandoraNu_tag);
    art::FindOneP< recob::Track > track_per_pfpart(pfparticle_handle, e, pandoraNu_tag);
    const std::map<size_t,  std::vector<size_t> > & map_primpfp_shwrpfp = fElectronEventSelectionAlg.get_pfp_id_showers_from_primary();
    const std::map<size_t,  std::vector<size_t> > & map_primpfp_trckpfp = fElectronEventSelectionAlg.get_pfp_id_tracks_from_primary();

    for (auto const& pfpindex : fElectronEventSelectionAlg.get_primary_indexes()) 
    {
      if (fElectronEventSelectionAlg.get_neutrino_candidate_passed().at(pfpindex)) 
      {
        nnuvtx++;
        TVector3 neutrino_vertex = fElectronEventSelectionAlg.get_neutrino_vertex().at(pfpindex);
        nuvtxx.emplace_back(neutrino_vertex.X());
        nuvtxy.emplace_back(neutrino_vertex.Y());
        nuvtxz.emplace_back(neutrino_vertex.Z());

        TVector3 center_of_charge = fElectronEventSelectionAlg.get_center_of_charge().at(pfpindex);
        center_of_charge_x.emplace_back(center_of_charge.X());
        center_of_charge_y.emplace_back(center_of_charge.Y());
        center_of_charge_z.emplace_back(center_of_charge.Z());

        recob::PFParticle const& pfp = pfparticle_handle->at(pfpindex);
        nuvtxpdg.emplace_back(pfp.PdgCode());

        std::vector<size_t> pf_ids = map_primpfp_shwrpfp.at(pfpindex);
        std::vector<TVector3> shwr_dir_primary(pf_ids.size());
        std::vector<Float_t> shwr_angle_primary(pf_ids.size());
        std::vector<Float_t> shwr_en_primary(pf_ids.size());

        for (size_t i =0; i< pf_ids.size(); ++i) {
          auto const& shower_obj = shower_per_pfpart.at(pf_ids[i]);
          shwr_dir_primary[i]=shower_obj->Direction();
          shwr_en_primary[i]=*std::max_element(shower_obj->Energy().begin(),shower_obj->Energy().end()); 
          shwr_angle_primary[i]=shower_obj->OpenAngle();
        }
        shwr_dir.emplace_back(shwr_dir_primary);
        shwr_en.emplace_back(shwr_en_primary);
        shwr_angle.emplace_back(shwr_angle_primary);


        std::vector<size_t> pf_idt = map_primpfp_trckpfp.at(pfpindex);
        std::vector<TVector3> trck_dir_primary(pf_idt.size());
        std::vector<Float_t> trck_dedxavg_primary(pf_ids.size());
        std::vector<Float_t> trck_len_primary(pf_ids.size());

        for (size_t i =0; i< pf_idt.size(); ++i) {
          auto const& track_obj   = track_per_pfpart.at(pf_idt[i]);
          trck_dir_primary[i]     = TVector3 (track_obj->Direction().second.x(),track_obj->Direction().second.y(),track_obj->Direction().second.z());
          trck_len_primary[i]     = track_obj->Length();
          trck_dedxavg_primary[i] = trackEnergy(track_obj, e); 
        }
        trck_dir.emplace_back(trck_dir_primary);
        trck_len.emplace_back(trck_len_primary);
        trck_dedxavg.emplace_back(trck_dedxavg_primary);
      }
    }
  }
}

void lee::ElectronSelectionAna::fillOticalTree(art::Event const & e){
  std::cout << "Filling optical information " << std::endl;
  art::InputTag optical_tag{"simpleFlashBeam"};
  auto const& optical_handle = e.getValidHandle<std::vector<recob::OpFlash>>(optical_tag);
  if(optical_handle->size())
  {
    nfls=0;
    flsTime.clear();
    flsPe.clear();
    flsYcenter.clear();
    flsZcenter.clear();
    flsYwidth.clear();
    flsZwidth.clear();

    std::map<size_t, int > op_flash_indexes = fElectronEventSelectionAlg.get_op_flash_indexes();
    nfls = op_flash_indexes.size();
    for(int ifl=0; ifl< nfls; ++ifl)
    {
      recob::OpFlash const& flash = optical_handle->at(op_flash_indexes[ifl]);
      flsTime.emplace_back(flash.Time());
      flsPe.emplace_back(flash.TotalPE());
      flsYcenter.emplace_back(flash.YCenter());
      flsZcenter.emplace_back(flash.ZCenter());
      flsYwidth.emplace_back(flash.YWidth());
      flsZwidth.emplace_back(flash.ZWidth());
    }
  }
}

DEFINE_ART_MODULE(lee::ElectronSelectionAna)
