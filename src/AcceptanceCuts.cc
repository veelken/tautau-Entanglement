#include "TauAnalysis/Entanglement/interface/AcceptanceCuts.h"

#include "TauAnalysis/Entanglement/interface/comp_visP4.h"            // comp_visP4()
#include "TauAnalysis/Entanglement/interface/findDecayProducts.h"     // findDecayProducts()
#include "TauAnalysis/Entanglement/interface/get_particles_of_type.h" // get_chargedHadrons(), get_photons()

namespace
{
  double
  readParameter(const edm::ParameterSet& cfg, const std::string& parameterName, double defaultValue)
  {
    if ( cfg.exists(parameterName) )
    {
      return cfg.getParameter<double>(parameterName);
    }
    else
    {
      return defaultValue;
    }
  }
}

ParticleSelector::ParticleSelector(const edm::ParameterSet& cfg)
  : minEta_(readParameter(cfg, "minEta", -1.e+6))
  , maxEta_(readParameter(cfg, "maxEta", +1.e+6))
  , minTheta_(readParameter(cfg, "minTheta", -1.e+6))
  , maxTheta_(readParameter(cfg, "maxTheta", +1.e+6))
  , minPt_(readParameter(cfg, "minPt", -1.e+6))
  , maxPt_(readParameter(cfg, "maxPt", +1.e+6))
  , minEnergy_(readParameter(cfg, "minEnergy", -1.e+6))
  , maxEnergy_(readParameter(cfg, "maxEnergy", +1.e+6))
{}

ParticleSelector::~ParticleSelector()
{}

bool
ParticleSelector::operator()(const reco::GenParticle& particle) const
{
  return this->operator()(particle.p4());
}

bool
ParticleSelector::operator()(const reco::Candidate::LorentzVector& p4) const
{
  if ( p4.eta()    > minEta_    && p4.eta()    < maxEta_    && 
       p4.theta()  > minTheta_  && p4.theta()  < maxTheta_  && 
       p4.pt()     > minPt_     && p4.pt()     < maxPt_     && 
       p4.energy() > minEnergy_ && p4.energy() < maxEnergy_ )
  {
    return true;
  }
  else
  {
    return false;
  }
}

AcceptanceCuts::AcceptanceCuts(const edm::ParameterSet& cfg)
  : chargedHadronSelector_(cfg.getParameter<edm::ParameterSet>("chargedHadrons"))
  , photonSelector_(cfg.getParameter<edm::ParameterSet>("photons"))
  , tauJetSelector_(cfg.getParameter<edm::ParameterSet>("tauJets"))
{}

AcceptanceCuts::~AcceptanceCuts()
{}

namespace
{
  bool
  passesAcceptanceCuts(const std::vector<const reco::GenParticle*>& particles, const ParticleSelector& selector)
  {
    for ( const reco::GenParticle* particle : particles )
    {
      if ( !selector(*particle) ) 
      {
        // CV: at least one particle given as function argument fails selection
        return false;
      }
    }
    // CV: all particles given as function argument pass selection
    return true;
  }
}

bool
AcceptanceCuts::operator()(const reco::GenParticle& tauPlus, const reco::GenParticle& tauMinus) const
{
  std::vector<const reco::GenParticle*> tauPlus_daughters;
  findDecayProducts(&tauPlus, tauPlus_daughters, false);
  reco::Candidate::LorentzVector visTauPlusP4 = comp_visP4(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_ch = get_chargedHadrons(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_y = get_photons(tauPlus_daughters);

  std::vector<const reco::GenParticle*> tauMinus_daughters;
  findDecayProducts(&tauMinus, tauMinus_daughters, false);
  reco::Candidate::LorentzVector visTauMinusP4 = comp_visP4(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_ch = get_chargedHadrons(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_y = get_photons(tauMinus_daughters);

  if ( tauJetSelector_(visTauPlusP4)                             &&
       passesAcceptanceCuts(tauPlus_ch,  chargedHadronSelector_) &&
       passesAcceptanceCuts(tauPlus_y,   photonSelector_)        &&
       tauJetSelector_(visTauMinusP4)                            &&
       passesAcceptanceCuts(tauMinus_ch, chargedHadronSelector_) &&
       passesAcceptanceCuts(tauMinus_y,  photonSelector_)        )
  {
    return true;
  }
  else
  {
    return false;
  }
}
