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

namespace lee {
  class ElectronSelectionAna;
}


class lee::ElectronSelectionAna : public art::EDAnalyzer {
public:
  explicit ElectronSelectionAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  void reconfigure(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  ElectronSelectionAna(ElectronSelectionAna const &) = delete;
  ElectronSelectionAna(ElectronSelectionAna &&) = delete;
  ElectronSelectionAna & operator = (ElectronSelectionAna const &) = delete;
  ElectronSelectionAna & operator = (ElectronSelectionAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.

  lee::ElectronEventSelectionAlg fElectronEventSelectionAlg;


};

void lee::ElectronSelectionAna::reconfigure(fhicl::ParameterSet const & p){
  fElectronEventSelectionAlg.reconfigure(p.get<fhicl::ParameterSet>("ElectronSelectionAlg"));
}


lee::ElectronSelectionAna::ElectronSelectionAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void lee::ElectronSelectionAna::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  
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

DEFINE_ART_MODULE(lee::ElectronSelectionAna)
