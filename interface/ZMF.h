#ifndef TauAnalysis_Entanglement_ZMF_h
#define TauAnalysis_Entanglement_ZMF_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"            // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"     // KinematicEvent
#include "TauAnalysis/Entanglement/interface/PolarimetricVector.h" // PolarimetricVector

class ZMF
{
 public:
  ZMF(const edm::ParameterSet& cfg);
  ~ZMF();

  KinematicEvent
  operator()(const KinematicEvent& evt);

 private:
  PolarimetricVector polarimetricVector_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_ZMF_h
