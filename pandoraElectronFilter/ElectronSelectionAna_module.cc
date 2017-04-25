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
#include "lardataobj/RecoBase/OpFlash.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

namespace lee {
  class ElectronSelectionAna;
}
class lee::ElectronSelectionAna : public art::EDAnalyzer {
public:
  explicit ElectronSelectionAna(fhicl::ParameterSet const & pset);
  virtual ~ElectronSelectionAna();

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


  // variables
  lee::ElectronEventSelectionAlg fElectronEventSelectionAlg;
  TFile*      fTFile;
  TTree*      fTree;


  //Run Subrun Event 
  Short_t    run;
  Short_t    subrun;
  Short_t    event;

  //Truth information
  std::vector<Short_t>   mcevts_truth;                  ///< number of neutrino Interactions in the spill
  std::vector<Short_t>   nuPDG_truth;                   ///< neutrino PDG code
  std::vector<Short_t>   ccnc_truth;                    ///< 0=CC 1=NC
  std::vector<Short_t>   mode_truth;                    ///< 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  std::vector<Float_t>   enu_truth;                     ///< true neutrino energy
  std::vector<Float_t>   nuvtxx_truth;                  ///< neutrino vertex x
  std::vector<Float_t>   nuvtxy_truth;                  ///< neutrino vertex y
  std::vector<Float_t>   nuvtxz_truth;                  ///< neutrino vertex z

  //PandoraNu information
  Bool_t               passed;                          ///< Return true if there is at least one neutrino vertex that passed
  Short_t              nnuvtx;                          ///< Number of PandoraNu neutrino candidate vertices
  std::vector<Float_t> nuvtxx;                          ///< x coordinate
  std::vector<Float_t> nuvtxy;                          ///< y coordinate
  std::vector<Float_t> nuvtxz;                          ///< z coordinate
  std::vector<Short_t> nuvtxpdg;                        ///< PDG code assigned by PandoraNu
  std::vector<Float_t> center_of_charge_x;              ///< x Center of deposited charge
  std::vector<Float_t> center_of_charge_y;              ///< y Center of deposited charge
  std::vector<Float_t> center_of_charge_z;              ///< z Center of deposited charge

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
  fTFile = new TFile("FlashOutput.root", "RECREATE");
  fTree  = tfs->make<TTree>("flashtree","FlashAnalysis Tree");

  //Set branches for (Run Subrun Event) 
  fTree->Branch("run",     &run,     "run/S"       );
  fTree->Branch("subrun",  &subrun,  "subrun/S"    );
  fTree->Branch("event",   &event,   "event/S"     );

  //Set branches for truth information
  fTree->Branch("mcevts_truth", &mcevts_truth,   "mcevts_truth/S"               );
  fTree->Branch("nuPDG_truth",  &nuPDG_truth,    "nuPDG_truth[mcevts_truth]/S"  );
  fTree->Branch("ccnc_truth",   &ccnc_truth,     "ccnc_truth[mcevts_truth]/S"   );
  fTree->Branch("mode_truth",   &mode_truth,     "mode_truth[mcevts_truth]/S"   );
  fTree->Branch("enu_truth",    &enu_truth,      "enu_truth[mcevts_truth]/F"    );
  fTree->Branch("nuvtxx_truth", &nuvtxx_truth,   "nuvtxx_truth[mcevts_truth]/F" );
  fTree->Branch("nuvtxy_truth", &nuvtxy_truth,   "nuvtxy_truth[mcevts_truth]/F" );
  fTree->Branch("nuvtxz_truth", &nuvtxz_truth,   "nuvtxz_truth[mcevts_truth]/F" );

  //Set branches for PandoraNU information
  fTree->Branch("passed",     &passed,     "passed/O"             );
  fTree->Branch("nnuvtx",     &nnuvtx,     "nnuvtx/S"             );
  fTree->Branch("nuvtxx",     &nuvtxx,     "nuvtxx[nnuvtx]/F"     );
  fTree->Branch("nuvtxy",     &nuvtxy,     "nuvtxy[nnuvtx]/F"     );
  fTree->Branch("nuvtxz",     &nuvtxz,     "nuvtxz[nnuvtx]/F"     );
  fTree->Branch("nuvtxpdg",   &nuvtxpdg,   "nuvtxpdg[nnuvtx]/S"   );
  fTree->Branch("center_of_charge_x",     &center_of_charge_x,     "center_of_charge_x[nnuvtx]/F"     );
  fTree->Branch("center_of_charge_y",     &center_of_charge_y,     "center_of_charge_y[nnuvtx]/F"     );
  fTree->Branch("center_of_charge_z",     &center_of_charge_z,     "center_of_charge_z[nnuvtx]/F"     );
  
  //Set branches for optical information
  fTree->Branch("nfls",       &nfls,       "nfls/S"             );
  fTree->Branch("flsTime",    &flsTime,    "flash_time[nfls]/F" );
  fTree->Branch("flsPe",      &flsPe,      "flsPe[nfls]/F"      );
  fTree->Branch("flsYcenter", &flsYcenter, "flsYcenter[nfls]/F" );
  fTree->Branch("flsZcenter", &flsZcenter, "flsZcenter[nfls]/F" );
  fTree->Branch("flsYwidth",  &flsYwidth,  "flsYwidth[nfls]/F"  );
  fTree->Branch("flsZwidth",  &flsZwidth,  "flsZwidth[nfls]/F"  );
  this->reconfigure(pset);
}
lee::ElectronSelectionAna::~ElectronSelectionAna()
{
  //write output file and tree
  std::cout << "Writing output..." << std::endl;
  fTFile->cd();
  fTree->Write("flashtree");
  fTFile->Close();
  std::cout << "Done!" << std::endl;
}
void lee::ElectronSelectionAna::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  fillTree(e);

  bool event_passed = fElectronEventSelectionAlg.eventSelected(e);
  if (event_passed){
    // Find out how many passing neutrino candidates there are:
    for (size_t i = 0; i < fElectronEventSelectionAlg.get_n_neutrino_candidates(); i ++){
      if (fElectronEventSelectionAlg.get_neutrino_candidate_passed().at(i)){
        std::cout << "Candidate " << i << " passed." << std::endl;
      }
    }
  }
  return;
}
void lee::ElectronSelectionAna::fillTree(art::Event const & e)
{

  // Fill run information 
  run    = e.run(); 
  subrun = e.subRun();
  event  = e.event();

  // Fill truth information
  //TODO

  // Fill PandoraNu information
  //TODO: PandoraNu PDG code of passed
  passed = fElectronEventSelectionAlg.eventSelected(e);
  if(passed)
  {
    // Find out how many passing neutrino candidates there are:
    for (size_t i = 0; i < fElectronEventSelectionAlg.get_n_neutrino_candidates(); i ++)
    {
      if (fElectronEventSelectionAlg.get_neutrino_candidate_passed().at(i))
      {
        nnuvtx++;

        TVector3 neutrino_vertex = fElectronEventSelectionAlg.get_neutrino_vertex().at(i);
        nuvtxx.push_back(neutrino_vertex.X());
        nuvtxy.push_back(neutrino_vertex.Y());
        nuvtxz.push_back(neutrino_vertex.Z());

        TVector3 center_of_charge = fElectronEventSelectionAlg.get_center_of_charge().at(i);
        center_of_charge_x.push_back(center_of_charge.X());
        center_of_charge_y.push_back(center_of_charge.Y());
        center_of_charge_z.push_back(center_of_charge.Z());
      }
    }
  }

  // Fill optical information
  art::InputTag optical_tag{"simpleFlashBeam"};
  auto const& optical_handle = e.getValidHandle<std::vector<recob::OpFlash>>(optical_tag);

  std::vector<int > op_flash_indexes = fElectronEventSelectionAlg.get_op_flash_indexes();
  nfls = op_flash_indexes.size();
  for(int ifl=0; ifl< nfls; ++ifl)
  {
    recob::OpFlash const& flash = optical_handle->at(ifl);
    flsTime.push_back(flash.Time());
    flsPe.push_back(flash.TotalPE());
    flsYcenter.push_back(flash.YCenter());
    flsZcenter.push_back(flash.ZCenter());
    flsYwidth.push_back(flash.YWidth());
    flsZwidth.push_back(flash.ZWidth());
  }

  fTree->Fill();
}

DEFINE_ART_MODULE(lee::ElectronSelectionAna)
