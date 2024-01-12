#ifndef TauAnalysis_Entanglement_AcceptanceCuts_h
#define TauAnalysis_Entanglement_AcceptanceCuts_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"       // edm::ParameterSet

#include "DataFormats/Candidate/interface/Candidate.h"        // reco::Candidate::LorentzVector
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle
#pragma GCC diagnostic pop

#include <vector>                                             // std::vector<>

class ParticleSelector
{
 public:
  ParticleSelector(const edm::ParameterSet& cfg);
  ~ParticleSelector();

  bool
  operator()(const reco::GenParticle& particle) const;
  bool
  operator()(const reco::Candidate::LorentzVector& p4) const;

 private:
  double minEta_;
  double maxEta_;
  double minTheta_;
  double maxTheta_;
  double minPt_;
  double maxPt_;
  double minEnergy_;
  double maxEnergy_;
};

class AcceptanceCuts
{
 public:
  AcceptanceCuts(const edm::ParameterSet& cfg);
  ~AcceptanceCuts();

  bool
  operator()(const reco::GenParticle& tauPlus, const reco::GenParticle& tauMinus) const;

 private:
  ParticleSelector chargedHadronSelector_;
  ParticleSelector photonSelector_;
  ParticleSelector tauJetSelector_;
};

#endif // TauAnalysis_Entanglement_AcceptanceCuts_h
