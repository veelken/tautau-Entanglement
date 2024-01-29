#include "TauAnalysis/Entanglement/interface/comp_normalDecayPlane.h"

#include "TauAnalysis/Entanglement/interface/constants.h"                  // mTau
#include "TauAnalysis/Entanglement/interface/getP4_rf.h"                   // getP4_rf()
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoBase.h" // pol::kTauPlus, pol::kTauMinus
#include "TauAnalysis/Entanglement/interface/square.h"                     // square()

#include <Math/Boost.h>                                                    // Boost

reco::Candidate::Vector
comp_normalDecayPlane(const KinematicEvent& evt, int tau)
{
  reco::Candidate::LorentzVector tauP4;
  reco::Candidate::LorentzVector visP4;
  if ( tau == pol::kTauPlus )
  {
    tauP4 = evt.tauPlusP4();
    visP4 = evt.visTauPlusP4();
  }
  else if ( tau == pol::kTauMinus )
  {
    tauP4 = evt.tauMinusP4();
    visP4 = evt.visTauMinusP4();
  }
  else assert(0);
  ROOT::Math::Boost boost_trf = ROOT::Math::Boost(tauP4.BoostToCM());
  const double sf = 1.01;
  double tauPx = sf*tauP4.px();
  double tauPy = sf*tauP4.py();
  double tauPz = sf*tauP4.pz();
  double tauE = std::sqrt(square(tauPx) + square(tauPy) + square(tauPz) + square(mTau));
  reco::Candidate::LorentzVector tauP4_dash(tauPx, tauPy, tauPz, tauE);
  reco::Candidate::LorentzVector tauP4_trf = getP4_rf(tauP4_dash, boost_trf);
  reco::Candidate::LorentzVector visP4_trf = getP4_rf(visP4, boost_trf);
  reco::Candidate::Vector normalDecayPlane = tauP4_trf.Vect().Cross(visP4_trf.Vect());
  return normalDecayPlane.unit();
}
