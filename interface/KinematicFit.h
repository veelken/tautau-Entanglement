#ifndef TauAnalysis_Entanglement_KinematicFit_h
#define TauAnalysis_Entanglement_KinematicFit_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"            // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"     // KinematicEvent
#include "TauAnalysis/Entanglement/interface/PolarimetricVector.h" // PolarimetricVector

class KinematicFit
{
 public:
  KinematicFit(const edm::ParameterSet& cfg);
  ~KinematicFit();

  KinematicEvent
  operator()(const KinematicEvent& evt);

 private:
  PolarimetricVector polarimetricVector_;

  int collider_;

  bool skip_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_KinematicFit_h
