#ifndef TauAnalysis_Entanglement_EntanglementNtupleProducer_h
#define TauAnalysis_Entanglement_EntanglementNtupleProducer_h

#include "FWCore/Framework/interface/one/EDAnalyzer.h"                   // edm::one::EDAnalyzer<>
#include "FWCore/Framework/interface/Event.h"                            // edm::Event
#include "FWCore/Framework/interface/EventSetup.h"                       // edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"                  // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h"                         // edm::InputTag<>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"            // reco::GenParticle
#pragma GCC diagnostic pop
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"         // reco::GenParticleCollection

#include "TauAnalysis/Entanglement/interface/AcceptanceCuts.h"           // AcceptanceCuts
#include "TauAnalysis/Entanglement/interface/ClassicSVfitInterface.h"    // ClassicSVfitInterface
#include "TauAnalysis/Entanglement/interface/EntanglementNtuple.h"       // EntanglementNtuple
#include "TauAnalysis/Entanglement/interface/GenKinematicEventBuilder.h" // GenKinematicEventBuilder
#include "TauAnalysis/Entanglement/interface/KinematicFit.h"             // KinematicFit
#include "TauAnalysis/Entanglement/interface/StartPosFinder.h"           // StartPosFinder
#include "TauAnalysis/Entanglement/interface/StartPosTIPCompatibility.h" // StartPosTIPCompatibility
#include "TauAnalysis/Entanglement/interface/ZMF.h"                      // ZMF

#include <TTree.h>                                                       // TTree

#include <vector>                                                        // std::vector<>
#include <string>                                                        // std::string

class EntanglementNtupleProducer : public edm::one::EDAnalyzer<>
{
 public:
  explicit EntanglementNtupleProducer(const edm::ParameterSet&);
  ~EntanglementNtupleProducer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::GenParticleCollection> token_;

  GenKinematicEventBuilder* genKineEvtBuilder_woSmearing_;
  GenKinematicEventBuilder* genKineEvtBuilder_wSmearing_;

  std::vector<StartPosFinder*> startPosFinders_;
  StartPosTIPCompatibility startPosTIPCompatibility_;

  enum { kTIP, kKinFit_chi2 };
  int startPosFinder_resolveSignAmbiguity_;

  KinematicFit* kinematicFit_;

  ClassicSVfitInterface* svFit_;

  ZMF* zmf_;

  AcceptanceCuts* acceptanceCuts_;

  int collider_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights_;
  std::vector<edm::EDGetTokenT<double>> tokenWeights_;

  TTree* ntuple_pi_pi_;
  EntanglementNtuple* ntupleFiller_pi_pi_;

  TTree* ntuple_pi_rho_;
  EntanglementNtuple* ntupleFiller_pi_rho_;

  TTree* ntuple_pi_a1_;
  EntanglementNtuple* ntupleFiller_pi_a1_;

  TTree* ntuple_rho_rho_;
  EntanglementNtuple* ntupleFiller_rho_rho_;

  TTree* ntuple_rho_a1_;
  EntanglementNtuple* ntupleFiller_rho_a1_;

  TTree* ntuple_a1_a1_;
  EntanglementNtuple* ntupleFiller_a1_a1_;

  TTree* ntuple_had_had_;
  EntanglementNtuple* ntupleFiller_had_had_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_EntanglementNtupleProducer_h
