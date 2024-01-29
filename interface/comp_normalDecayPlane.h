#ifndef TauAnalysis_Entanglement_comp_normalDecayPlane_h
#define TauAnalysis_Entanglement_comp_normalDecayPlane_h

#include "DataFormats/Candidate/interface/Candidate.h"         // Candidate::Vector

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent

reco::Candidate::Vector
comp_normalDecayPlane(const KinematicEvent& evt, int tau);

#endif // TauAnalysis_Entanglement_comp_normalDecayPlane_h
