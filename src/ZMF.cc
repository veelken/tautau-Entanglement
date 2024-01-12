#include "TauAnalysis/Entanglement/interface/ZMF.h"

#include <iostream> // std::cout
#include <string>   // std::string

ZMF::ZMF(const edm::ParameterSet& cfg)
  : polarimetricVector_(cfg)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{}

ZMF::~ZMF()
{}

KinematicEvent
ZMF::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<ZMF::operator()>:\n"; 
  }

  KinematicEvent kineEvt_zmf = kineEvt;
  kineEvt_zmf.recoilP4_ = kineEvt.visTauPlusP4_ + kineEvt.visTauMinusP4_;
  kineEvt_zmf.recoilCov_ = kineEvt.visTauPlusCov_ + kineEvt.visTauMinusCov_;
  kineEvt_zmf.nuTauPlusP4_ = reco::Candidate::LorentzVector(0.,0.,0.,0.);
  kineEvt_zmf.nuTauPlusP4_isValid_ = true;
  kineEvt_zmf.nuTauPlusCov_ = math::Matrix3x3();
  kineEvt_zmf.tauPlusP4_ = kineEvt.visTauPlusP4_;
  kineEvt_zmf.tauPlusP4_isValid_ = true;
  kineEvt_zmf.nuTauMinusP4_ = reco::Candidate::LorentzVector(0.,0.,0.,0.);
  kineEvt_zmf.nuTauMinusP4_isValid_ = true;       
  kineEvt_zmf.nuTauMinusCov_ = math::Matrix3x3();
  kineEvt_zmf.tauMinusP4_ = kineEvt.visTauMinusP4_;
  kineEvt_zmf.tauMinusP4_isValid_ = true;
  reco::Candidate::Vector hPlus = polarimetricVector_(kineEvt_zmf, pol::kTauPlus);
  kineEvt_zmf.hPlus_ = hPlus;
  kineEvt_zmf.hPlus_isValid_ = true;
  reco::Candidate::Vector hMinus = polarimetricVector_(kineEvt_zmf, pol::kTauMinus);
  kineEvt_zmf.hMinus_ = hMinus;
  kineEvt_zmf.hMinus_isValid_ = true;

  return kineEvt_zmf;
}

