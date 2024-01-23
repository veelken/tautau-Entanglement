#include "TauAnalysis/Entanglement/interface/ClassicSVfitInterface.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#include "DataFormats/TauReco/interface/PFTau.h"                          // reco::PFTau::hadronicDecayMode
#pragma GCC diagnostic pop

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"              // ClassicSVfit
#include "TauAnalysis/ClassicSVfit/interface/HistogramAdapterDiTauSpin.h" // HistogramAdapterDiTauSpin
#include "TauAnalysis/ClassicSVfit/interface/MeasuredEvent.h"             // MeasuredEvent
#include "TauAnalysis/ClassicSVfit/interface/MeasuredMEt.h"               // MeasuredMEt
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"         // MeasuredTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"     // HistogramAdapterDiTau

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"                 // kLHC, kSuperKEKB, mHiggs
#include "TauAnalysis/Entanglement/interface/convert_to_TMatrixD.h"       // convert_to_TMatrixD()
#include "TauAnalysis/Entanglement/interface/get_decayMode.h"             // is3Prong()

#include <iostream>                                                       // std::cout
#include <string>                                                         // std::string

using namespace classic_svFit;

ClassicSVfitInterface::ClassicSVfitInterface(const edm::ParameterSet& cfg)
  : resolutions_(nullptr)
  , svFitAlgo_(nullptr)
  , histogramAdapter_(nullptr)
  , applyHiggsMassConstraint_(cfg.getParameter<bool>("applyHiggsMassConstraint"))
  , skip_(cfg.getParameter<bool>("skip"))
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  edm::ParameterSet cfg_resolutions = cfg.getParameter<edm::ParameterSet>("resolutions");
  resolutions_ = new Resolutions(cfg_resolutions);

  std::string collider = cfg.getParameter<std::string>("collider");
  if      ( collider == "LHC"       ) collider_ = kLHC;
  else if ( collider == "SuperKEKB" ) collider_ = kSuperKEKB;
  else throw cmsException("KinematicFit", __LINE__)
    << "Invalid Configuration parameter 'collider' = " << collider << " !!\n";

  svFitAlgo_ = new ClassicSVfit();
  svFitAlgo_->enableTauFlightLength();
  //svFitAlgo_->disableTauFlightLength();  
  if ( applyHiggsMassConstraint_ )
  {
    std::cout << "Enabling SVfit mass-constraint (mH = " << mHiggs << " GeV).\n";
    svFitAlgo_->enableDiTauMassConstraint(mHiggs);
  }
  else
  {
    std::cout << "Disabling SVfit mass-constraint.\n";
    svFitAlgo_->disableDiTauMassConstraint();
    if ( collider_ == kLHC )
    {
      svFitAlgo_->enableLogM(6.);
    }
    else
    {
      svFitAlgo_->disableLogM();
    }
  }
  //svFitAlgo_->setMaxObjFunctionCalls(100000); // CV: default is 100000 evaluations of integrand per event
//svFitAlgo_->setMaxObjFunctionCalls(1000);
  histogramAdapter_ = new HistogramAdapterDiTauSpin();
  svFitAlgo_->setHistogramAdapter(histogramAdapter_);
//svFitAlgo_->setVerbosity(2);
}

ClassicSVfitInterface::~ClassicSVfitInterface()
{
  delete resolutions_;
  delete svFitAlgo_;
  // CV: histogramAdapter deleted in destructor of ClassicSVfit class !!
  //delete histogramAdapter_;
}

namespace
{
  std::vector<MeasuredHadTauDecayProduct>
  buildMeasuredHadTauDecayProduct(const std::vector<KinematicParticle>& decayProducts)
  {
    std::vector<MeasuredHadTauDecayProduct> measuredHadTauDecayProducts;
    for ( const KinematicParticle& decayProduct : decayProducts )
    {
      if ( std::abs(decayProduct.pdgId()) == 11 || std::abs(decayProduct.pdgId()) == 13 ) continue;
      if ( decayProduct.charge() != 0 || decayProduct.pdgId() == 111 )
      {
        const reco::Candidate::LorentzVector& decayProductP4 = decayProduct.p4();
        MeasuredHadTauDecayProduct measuredHadTauDecayProduct(
          decayProduct.charge(),
          decayProductP4.pt(), decayProductP4.eta(), decayProductP4.phi(), decayProductP4.mass());
        measuredHadTauDecayProducts.push_back(measuredHadTauDecayProduct);
      }
    }
    return measuredHadTauDecayProducts;
  }

  MeasuredTauLepton
  buildMeasuredTauLepton(int charge, const reco::Candidate::LorentzVector& visTauP4, int decayMode,
                         const reco::Candidate::Point& decayVertex, const math::Matrix3x3& decayVertexCov, 
                         const Resolutions& resolutions,
                         const std::vector<MeasuredHadTauDecayProduct>& measuredHadTauDecayProducts)
  {
    MeasuredTauLepton measuredTauLepton;
    if ( is3Prong(decayMode) )
    {
      measuredTauLepton = MeasuredTauLepton(
        MeasuredTauLepton::kTauToHadDecay, 
        charge, visTauP4.pt(), visTauP4.eta(), visTauP4.phi(), visTauP4.mass(),
        decayVertex, convert_to_TMatrixD(decayVertexCov),
        decayMode, &measuredHadTauDecayProducts);
    }
    else
    {
      measuredTauLepton = MeasuredTauLepton(
        MeasuredTauLepton::kTauToHadDecay, 
        charge, visTauP4.pt(), visTauP4.eta(), visTauP4.phi(), visTauP4.mass(),
        decayVertex, 1.e+3, resolutions.tipResolution_perp(),
        decayMode, &measuredHadTauDecayProducts);
    }
    return measuredTauLepton;
  }

  reco::Candidate::Vector
  buildPolarimeterVector(double h_r, double h_n, double h_k)
  {
    reco::Candidate::Vector h(h_r, h_n, h_k);
    return h.unit();
  }
}

KinematicEvent
ClassicSVfitInterface::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<ClassicSVfitInterface::operator()>:\n"; 
  }

  if ( skip_ )
  {
    return kineEvt;
  }

  KinematicEvent kineEvt_svFit = kineEvt;

  reco::Candidate::Point tauPlusDecayVertex;
  if ( is3Prong(kineEvt.tauPlus_decayMode()) )
  {
    assert(kineEvt.svTauPlus_isValid());
    tauPlusDecayVertex = kineEvt.svTauPlus();
  }
  else
  {
    tauPlusDecayVertex = kineEvt.tipPCATauPlus();
  }
  std::vector<MeasuredHadTauDecayProduct> measuredTauPlusDecayProducts = buildMeasuredHadTauDecayProduct(
    kineEvt.daughtersTauPlus());
  MeasuredTauLepton measuredTauPlus = buildMeasuredTauLepton(
    +1, kineEvt.visTauPlusP4(), kineEvt.tauPlus_decayMode(),
    tauPlusDecayVertex, kineEvt.svTauPlusCov(),
    *resolutions_,
    measuredTauPlusDecayProducts);

  reco::Candidate::Point tauMinusDecayVertex;
  if ( is3Prong(kineEvt.tauMinus_decayMode()) )
  {
    assert(kineEvt.svTauMinus_isValid());
    tauMinusDecayVertex = kineEvt.svTauMinus();
  }
  else
  {
    tauMinusDecayVertex = kineEvt.tipPCATauMinus();
  }
  std::vector<MeasuredHadTauDecayProduct> measuredTauMinusDecayProducts = buildMeasuredHadTauDecayProduct(
    kineEvt.daughtersTauMinus());
  MeasuredTauLepton measuredTauMinus = buildMeasuredTauLepton(
    -1, kineEvt.visTauMinusP4(), kineEvt.tauMinus_decayMode(),
    tauMinusDecayVertex, kineEvt.svTauMinusCov(),
    *resolutions_,
    measuredTauMinusDecayProducts);

  std::vector<MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(measuredTauPlus);
  measuredTauLeptons.push_back(measuredTauMinus);

  MeasuredMEt measuredMEt;
  if ( collider_ == kLHC )
  {
    reco::Candidate::LorentzVector MEt = kineEvt.recoilP4() - (kineEvt.visTauPlusP4() + kineEvt.visTauMinusP4());
    TMatrixD covMEt = convert_to_TMatrixD(kineEvt.recoilCov().Sub<math::Matrix2x2>(0,0));
    measuredMEt = MeasuredMEt(MEt.px(), MEt.py(), covMEt);
  }
  else if ( collider_ == kSuperKEKB )
  {
    reco::Candidate::LorentzVector MEt = kineEvt.recoilP4() - (kineEvt.visTauPlusP4() + kineEvt.visTauMinusP4());
    TMatrixD covMEt = convert_to_TMatrixD(kineEvt.recoilCov());
    measuredMEt = MeasuredMEt(MEt.px(), MEt.py(), MEt.pz(), MEt.energy(), covMEt);
  } else assert(0);

  MeasuredEvent measuredEvent(measuredTauLeptons, { measuredMEt }, kineEvt.pv(), convert_to_TMatrixD(kineEvt.pvCov()));

  svFitAlgo_->setStartPosition(kineEvt.tauPlusP4(), kineEvt.tauMinusP4());
  svFitAlgo_->integrate(measuredEvent);
  bool isValidSolution = svFitAlgo_->isValidSolution();

  if ( isValidSolution  )
  {
    const HistogramAdapterDiTauSpin* fittedDiTau = dynamic_cast<const HistogramAdapterDiTauSpin*>(svFitAlgo_->getHistogramAdapter());
    assert(fittedDiTau);
    const HistogramAdapterTau* fittedTau1 = fittedDiTau->tau1();
    const HistogramAdapterTau* fittedTau2 = fittedDiTau->tau2();

    reco::Candidate::LorentzVector tauPlusP4, tauMinusP4;
    if ( fittedTau1->getCharge() == +1 && fittedTau2->getCharge() == -1 )
    {
      tauPlusP4 = fittedTau1->getP4();
      tauMinusP4 = fittedTau2->getP4();
    }
    else if ( fittedTau1->getCharge() == -1 && fittedTau2->getCharge() == +1 )
    {
      tauPlusP4 = fittedTau2->getP4();
      tauMinusP4 = fittedTau1->getP4();
    } else assert(0);

    kineEvt_svFit.nuTauPlusP4_ = tauPlusP4 - kineEvt.visTauPlusP4();
    kineEvt_svFit.nuTauPlusP4_isValid_ = true;
    kineEvt_svFit.tauPlusP4_ = tauPlusP4;
    kineEvt_svFit.tauPlusP4_isValid_ = true;
    kineEvt_svFit.hPlus_ = buildPolarimeterVector(fittedDiTau->getBp_r(), fittedDiTau->getBp_n(), fittedDiTau->getBp_k());
    kineEvt_svFit.hPlus_isValid_ = true;
    kineEvt_svFit.nuTauMinusP4_ = tauMinusP4 - kineEvt.visTauMinusP4();
    kineEvt_svFit.nuTauMinusP4_isValid_ = true;       
    kineEvt_svFit.tauMinusP4_ = tauMinusP4;
    kineEvt_svFit.tauMinusP4_isValid_ = true;
    kineEvt_svFit.hMinus_ = buildPolarimeterVector(fittedDiTau->getBm_r(), fittedDiTau->getBm_n(), fittedDiTau->getBm_k());
    kineEvt_svFit.hMinus_isValid_ = true;
    kineEvt_svFit.svFit_isValid_ = true;
  } 
  else
  {
    kineEvt_svFit.svFit_isValid_ = false;
  }

  return kineEvt_svFit;
}

