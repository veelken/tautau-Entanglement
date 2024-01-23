#ifndef TauAnalysis_Entanglement_ClassicSVfitInterface_h
#define TauAnalysis_Entanglement_ClassicSVfitInterface_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"                   // edm::ParameterSet
#include "DataFormats/Candidate/interface/Candidate.h"                    // reco::Candidate::LorentzVector

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"              // ClassicSVfit
#include "TauAnalysis/ClassicSVfit/interface/HistogramAdapterDiTauSpin.h" // HistogramAdapterDiTauSpin

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"            // KinematicEvent
#include "TauAnalysis/Entanglement/interface/PolarimetricVector.h"        // PolarimetricVector
#include "TauAnalysis/Entanglement/interface/Resolutions.h"               // Resolutions

class ClassicSVfitInterface
{
 public:
  ClassicSVfitInterface(const edm::ParameterSet& cfg);
  ~ClassicSVfitInterface();

  KinematicEvent
  operator()(const KinematicEvent& evt);

 private:
  Resolutions* resolutions_;
  int collider_;

  bool applyHiggsMassConstraint_;

  ClassicSVfit* svFitAlgo_;
  classic_svFit::HistogramAdapterDiTauSpin* histogramAdapter_;

  bool skip_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_ClassicSVfitInterface_h
