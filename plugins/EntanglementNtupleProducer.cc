#include "TauAnalysis/Entanglement/plugins/EntanglementNtupleProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"                 // edm::Service<>
#include "CommonTools/UtilAlgos/interface/TFileService.h"             // TFileService

#include "DataFormats/Common/interface/Handle.h"                      // edm::Handle<>
#include "DataFormats/Math/interface/angle.h"                         // angle()

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#include "DataFormats/TauReco/interface/PFTau.h"                      // reco::PFTau::hadronicDecayMode
#pragma GCC diagnostic pop

#include "TauAnalysis/Entanglement/interface/cmsException.h"          // cmsException
#include "TauAnalysis/Entanglement/interface/comp_visP4.h"            // comp_visP4()
#include "TauAnalysis/Entanglement/interface/constants.h"             // kLHC, kSuperKEKB
#include "TauAnalysis/Entanglement/interface/findDecayProducts.h"     // findDecayProducts()
#include "TauAnalysis/Entanglement/interface/findLastTau.h"           // findLastTau()
#include "TauAnalysis/Entanglement/interface/get_decayMode.h"         // get_decayMode(), isHadTauDecay()
#include "TauAnalysis/Entanglement/interface/get_particles_of_type.h" // get_chargedKaons(), get_neutralKaons(), get_photons()

#include <TString.h>                                                  // Form()

#include <iostream>                                                   // std::cout

EntanglementNtupleProducer::EntanglementNtupleProducer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , genKineEvtBuilder_woSmearing_(nullptr)
  , genKineEvtBuilder_wSmearing_(nullptr)
  , startPosTIPCompatibility_(cfg)
  , kinematicFit_(nullptr)
  , svFit_(nullptr)
  , acceptanceCuts_(nullptr)
  , ntuple_pi_pi_(nullptr)
  , ntupleFiller_pi_pi_(nullptr)
  , ntuple_pi_rho_(nullptr)
  , ntupleFiller_pi_rho_(nullptr)
  , ntuple_pi_a1_(nullptr)
  , ntupleFiller_pi_a1_(nullptr)
  , ntuple_rho_rho_(nullptr)
  , ntupleFiller_rho_rho_(nullptr)
  , ntuple_rho_a1_(nullptr)
  , ntupleFiller_rho_a1_(nullptr)
  , ntuple_a1_a1_(nullptr)
  , ntupleFiller_a1_a1_(nullptr)
  , ntuple_had_had_(nullptr)
  , ntupleFiller_had_had_(nullptr)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenParticleCollection>(src_);

  edm::ParameterSet cfg_resolutions = cfg.getParameter<edm::ParameterSet>("resolutions");

  bool applySmearing = cfg.getParameter<bool>("applySmearing");

  edm::ParameterSet cfg_woSmearing = cfg;
  cfg_woSmearing.addParameter<edm::ParameterSet>("resolutions", cfg_resolutions);
  cfg_woSmearing.addParameter<bool>("applySmearing", false);
  genKineEvtBuilder_woSmearing_ = new GenKinematicEventBuilder(cfg_woSmearing);

  edm::ParameterSet cfg_wSmearing = cfg;
  cfg_wSmearing.addParameter<edm::ParameterSet>("resolutions", cfg_resolutions);
  cfg_wSmearing.addParameter<bool>("applySmearing", applySmearing);
  genKineEvtBuilder_wSmearing_ = new GenKinematicEventBuilder(cfg_wSmearing);
  
  std::string hAxis = cfg.getParameter<std::string>("hAxis");

  std::string collider = cfg.getParameter<std::string>("collider");
  if      ( collider == "LHC"       ) collider_ = kLHC;
  else if ( collider == "SuperKEKB" ) collider_ = kSuperKEKB;
  else throw cmsException("EntanglementNtupleProducer", __LINE__)
    << "Invalid Configuration parameter 'collider' = " << collider << " !!\n";

  edm::ParameterSet cfg_startPosFinder = cfg.getParameter<edm::ParameterSet>("startPosFinder");
  cfg_startPosFinder.addParameter<edm::ParameterSet>("resolutions", cfg_resolutions);
  cfg_startPosFinder.addParameter<std::string>("hAxis", hAxis);
  cfg_startPosFinder.addParameter<std::string>("collider", collider);
  cfg_startPosFinder.addUntrackedParameter<int>("verbosity", verbosity_);
  cfg_startPosFinder.addUntrackedParameter<bool>("cartesian", cartesian_);
  std::vector<int> startPos_algos = cfg_startPosFinder.getParameter<std::vector<int>>("algos");
  for ( int startPos_algo : startPos_algos )
  {
    edm::ParameterSet cfg_startPosFinder_algo = cfg_startPosFinder;
    cfg_startPosFinder_algo.addParameter<int>("algo", startPos_algo);
    startPosFinders_.push_back(new StartPosFinder(cfg_startPosFinder_algo));
  }
  std::string resolveSignAmbiguity = cfg_startPosFinder.getParameter<std::string>("resolveSignAmbiguity");
  if      ( resolveSignAmbiguity == "TIP"         ) startPosFinder_resolveSignAmbiguity_ = kTIP;
  else if ( resolveSignAmbiguity == "kinFit_chi2" ) startPosFinder_resolveSignAmbiguity_ = kKinFit_chi2;
  else throw cmsException("StartPosFinder::StartPosFinder", __LINE__) 
    << "Invalid Configuration parameter 'resolveSignAmbiguity' = " << resolveSignAmbiguity << " !!\n";

  edm::ParameterSet cfg_kinematicFit = cfg.getParameter<edm::ParameterSet>("kinematicFit");
  cfg_kinematicFit.addParameter<std::string>("hAxis", hAxis);
  cfg_kinematicFit.addParameter<std::string>("collider", collider);
  cfg_kinematicFit.addUntrackedParameter<int>("verbosity", verbosity_);
  cfg_kinematicFit.addUntrackedParameter<bool>("cartesian", cartesian_);
  kinematicFit_ = new KinematicFit(cfg_kinematicFit);

  edm::ParameterSet cfg_svFit = cfg.getParameter<edm::ParameterSet>("svFit");
  cfg_svFit.addParameter<edm::ParameterSet>("resolutions", cfg_resolutions);
  cfg_svFit.addParameter<std::string>("collider", collider);
  cfg_svFit.addUntrackedParameter<int>("verbosity", verbosity_);
  cfg_svFit.addUntrackedParameter<bool>("cartesian", cartesian_);
  svFit_ = new ClassicSVfitInterface(cfg_svFit);

  edm::ParameterSet cfg_zmf;
  cfg_zmf.addUntrackedParameter<int>("verbosity", verbosity_);
  cfg_zmf.addUntrackedParameter<bool>("cartesian", cartesian_);
  zmf_ = new ZMF(cfg_zmf);

  acceptanceCuts_ = new AcceptanceCuts(cfg.getParameter<edm::ParameterSet>("acceptanceCuts"));

  srcWeights_ = cfg.getParameter<vInputTag>("srcEvtWeights");
  for ( const edm::InputTag& srcWeight : srcWeights_ )
  {
    tokenWeights_.push_back(consumes<double>(srcWeight));
  }
}

EntanglementNtupleProducer::~EntanglementNtupleProducer()
{
  delete genKineEvtBuilder_woSmearing_;
  delete genKineEvtBuilder_wSmearing_;

  for ( StartPosFinder* startPosFinder : startPosFinders_ )
  {
    delete startPosFinder;
  }

  delete kinematicFit_;

  delete svFit_;

  delete zmf_;

  delete acceptanceCuts_;

  // CV: don't delete TTree objects, as these are handled by TFileService

  delete ntupleFiller_pi_pi_;
  delete ntupleFiller_pi_rho_;
  delete ntupleFiller_pi_a1_;
  delete ntupleFiller_rho_rho_;
  delete ntupleFiller_rho_a1_;
  delete ntupleFiller_a1_a1_;
  delete ntupleFiller_had_had_;
}

void EntanglementNtupleProducer::beginJob()
{
  edm::Service<TFileService> fs;

  ntuple_pi_pi_ = fs->make<TTree>("pi_pi", "pi_pi");
  ntupleFiller_pi_pi_ = new EntanglementNtuple(ntuple_pi_pi_);
  ntuple_pi_rho_ = fs->make<TTree>("pi_rho", "pi_rho");
  ntupleFiller_pi_rho_ = new EntanglementNtuple(ntuple_pi_rho_);
  ntuple_pi_a1_ = fs->make<TTree>("pi_a1", "pi_a1");
  ntupleFiller_pi_a1_ = new EntanglementNtuple(ntuple_pi_a1_);
  ntuple_rho_rho_ = fs->make<TTree>("rho_rho", "rho_rho");
  ntupleFiller_rho_rho_ = new EntanglementNtuple(ntuple_rho_rho_);
  ntuple_rho_a1_ = fs->make<TTree>("rho_a1", "rho_a1");
  ntupleFiller_rho_a1_ = new EntanglementNtuple(ntuple_rho_a1_);
  ntuple_a1_a1_ = fs->make<TTree>("a1_a1", "a1_a1");
  ntupleFiller_a1_a1_ = new EntanglementNtuple(ntuple_a1_a1_);
  ntuple_had_had_ = fs->make<TTree>("had_had", "had_had");
  ntupleFiller_had_had_ = new EntanglementNtuple(ntuple_had_had_);
}

namespace
{
  bool
  isHigherTIPCompatibility(const KinematicEvent& kineEvt1, const KinematicEvent& kineEvt2)
  {
    return kineEvt1.startPosTIPCompatibility() > kineEvt2.startPosTIPCompatibility();
  }
}

void EntanglementNtupleProducer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<EntanglementNtupleProducer::analyze>:\n";
    std::cout << " moduleLabel = " << moduleLabel_ << "\n";
  }

  double evtWeight = 1.;
  for ( const edm::EDGetTokenT<double>& tokenWeight : tokenWeights_ )
  {
    edm::Handle<double> weight;
    evt.getByToken(tokenWeight, weight);
    evtWeight *= (*weight);
  }
  if ( verbosity_ >= 1 )
  {
    std::cout << "evtWeight = " << evtWeight << "\n";
  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(token_, genParticles);

  KinematicEvent kineEvt_gen = (*genKineEvtBuilder_woSmearing_)(*genParticles);
  if ( verbosity_ >= 1 )
  {
    printKinematicEvent("kineEvt_gen", kineEvt_gen, verbosity_, cartesian_);
  }

  KinematicEvent kineEvt_gen_smeared = (*genKineEvtBuilder_wSmearing_)(*genParticles);
  if ( verbosity_ >= 1 )
  {
    printKinematicEvent("kineEvt_gen_smeared", kineEvt_gen_smeared, verbosity_, cartesian_);
  }

  std::vector<KinematicEvent> kineEvts_startPos;
  for ( StartPosFinder* startPosFinder : startPosFinders_ )
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "executing startPosFinder algo #" << startPosFinder->get_algo() << "...";
    }

    std::vector<KinematicEvent> startPosFinder_gen_solutions  = (*startPosFinder)(kineEvt_gen);
    size_t numSolutions_gen = startPosFinder_gen_solutions.size();
    for ( size_t idxSolution = 0; idxSolution < numSolutions_gen; ++idxSolution )
    {
      KinematicEvent& kineEvt_startPos = startPosFinder_gen_solutions.at(idxSolution);
      double tauPlus_residual = angle(kineEvt_startPos.tauPlusP4().Vect(), kineEvt_gen.tauPlusP4().Vect());
      double tauMinus_residual = angle(kineEvt_startPos.tauMinusP4().Vect(), kineEvt_gen.tauMinusP4().Vect());
      // CV: large negative values of tipCompatibility indicate large residuals,
      //     while a tipCompatibility of zero indicates a perfect match
      kineEvt_startPos.startPosTIPCompatibility_ = -(pow(tauPlus_residual, 2) + pow(tauMinus_residual, 2));
    }
    std::sort(startPosFinder_gen_solutions.begin(), startPosFinder_gen_solutions.end(), isHigherTIPCompatibility);
    int gen_startPosSign = startPosFinder_gen_solutions.front().startPosSign();
    if ( verbosity_ >= 3 )
    {
      std::cout << "gen_startPosSign = " << gen_startPosSign << "\n";
    }

    std::vector<KinematicEvent> startPosFinder_gen_smeared_solutions  = (*startPosFinder)(kineEvt_gen_smeared);
    size_t numSolutions_gen_smeared = startPosFinder_gen_smeared_solutions.size();
    if ( verbosity_ >= 1 )
    {
      std::cout << "#solutions = " << numSolutions_gen_smeared << "\n";
    }
    for ( size_t idxSolution = 0; idxSolution < numSolutions_gen_smeared; ++idxSolution )
    {
      KinematicEvent& kineEvt_startPos = startPosFinder_gen_smeared_solutions.at(idxSolution);
      std::string label = Form("kineEvt_startPos (algo #%i, solution #%lu)", startPosFinder->get_algo(), idxSolution);
      kineEvt_startPos.label_ = label;
      kineEvt_startPos.startPosTIPCompatibility_ = startPosTIPCompatibility_(kineEvt_startPos);
      if ( verbosity_ >= 3 )
      {
        std::cout << "solution #" << idxSolution << ": tipCompatibility = " << kineEvt_startPos.startPosTIPCompatibility_ << "\n";
      }
      kineEvt_startPos.startPosSign_isCorrect_ = ( kineEvt_startPos.startPosSign_ == gen_startPosSign ) ? true : false;
      kineEvts_startPos.push_back(kineEvt_startPos);
    }
  }

  // CV: sort solutions by decreasing compatibility with transverse impact parameters;
  //     the transverse impact parameter-compatibility of solutions is computed as described in the paper arXiv:hep-ph/9307269 
  std::sort(kineEvts_startPos.begin(), kineEvts_startPos.end(), isHigherTIPCompatibility);

  // CV: choose solution with best transverse impact parameter-compatibility
  //    (and run kinematic fit only for this solution)
  if ( startPosFinder_resolveSignAmbiguity_ == kTIP )
  {
    if ( kineEvts_startPos.size() >= 1 ) kineEvts_startPos = { kineEvts_startPos.front() };
  }

  KinematicEvent kineEvt_startPos_bestfit;
  KinematicEvent kineEvt_kinFit_bestfit;
  double kinFitChi2_bestfit = -1.;
  int kinFitStatus_bestfit = -1;
  bool isFirst = true;
  for ( const KinematicEvent& kineEvt_startPos : kineEvts_startPos )
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "executing KinematicFit for " << kineEvt_startPos.label() << "...\n";
      printKinematicEvent("kineEvt_startPos", kineEvt_startPos, verbosity_, cartesian_);
    }

    KinematicEvent kineEvt_kinFit = (*kinematicFit_)(kineEvt_startPos);
    if ( verbosity_ >= 1 )
    {
      printKinematicEvent("kineEvt_kinFit", kineEvt_kinFit, verbosity_, cartesian_);
    }
 
    bool isBetterFit = kineEvt_kinFit.kinFitStatus() >  kinFitStatus_bestfit || 
                      (kineEvt_kinFit.kinFitStatus() == kinFitStatus_bestfit && kineEvt_kinFit.kinFitChi2() < kinFitChi2_bestfit);
    if ( verbosity_ >= 1 )
    {
      std::cout << "isBetterFit = " << isBetterFit << "\n";
    }
    if ( isFirst || isBetterFit )
    {
      kineEvt_startPos_bestfit = kineEvt_startPos;
      kineEvt_kinFit_bestfit = kineEvt_kinFit;
      kinFitChi2_bestfit = kineEvt_kinFit.kinFitChi2();
      kinFitStatus_bestfit = kineEvt_kinFit.kinFitStatus();
      isFirst = false;
    }
  }
  if ( !(kinFitStatus_bestfit == 1 || (kinFitStatus_bestfit == 0 && kinFitChi2_bestfit < 1.e+2)) && verbosity_ >= 0 )
  {
    std::cerr << "WARNING: KinematicFit failed to converge !!" << std::endl;
  }

  if ( verbosity_ >= 1 )
  {
    std::cout << "best fit:\n";
    printKinematicEvent("kineEvt_startPos", kineEvt_startPos_bestfit, verbosity_, cartesian_);
    printKinematicEvent("kineEvt_kinFit", kineEvt_kinFit_bestfit, verbosity_, cartesian_);
  }

  KinematicEvent kineEvt_svFit;
  if ( kineEvts_startPos.size() >= 1 )
  {
    const KinematicEvent kineEvt_startPos = kineEvts_startPos[0];
    if ( isHadTauDecay(kineEvt_startPos.tauPlus_decayMode()) && isHadTauDecay(kineEvt_startPos.tauMinus_decayMode()) )
    {
      kineEvt_svFit = (*svFit_)(kineEvt_startPos);
      if ( verbosity_ >= 1 )
      {
        printKinematicEvent("kineEvt_svFit", kineEvt_svFit, verbosity_, cartesian_);
      }
    }
  }

  KinematicEvent kineEvt_zmf = (*zmf_)(kineEvt_gen_smeared);
  if ( verbosity_ >= 1 )
  {
    printKinematicEvent("kineEvt_zmf", kineEvt_zmf, verbosity_, cartesian_);
  }

  const reco::GenParticle* tauPlus  = nullptr;
  const reco::GenParticle* tauMinus = nullptr;
  double tauPlus_enLoss = 0.;
  double tauMinus_enLoss = 0.;
  for ( const reco::GenParticle& genParticle : *genParticles )
  {
    if ( genParticle.pdgId() == -15 && !tauPlus  )
    {
      const double tauPlus_enBefore = genParticle.energy();
      tauPlus  = findLastTau(&genParticle);
      const double tauPlus_enAfter = tauPlus->energy();
      tauPlus_enLoss = std::max(tauPlus_enBefore - tauPlus_enAfter, 0.);
    }
    if ( genParticle.pdgId() == +15 && !tauMinus )
    {
      const double tauMinus_enBefore = genParticle.energy();
      tauMinus = findLastTau(&genParticle);
      const double tauMinus_enAfter = tauMinus->energy();
      tauMinus_enLoss = std::max(tauMinus_enBefore - tauMinus_enAfter, 0.);
    }
  }
  if ( !(tauPlus && tauMinus) ) 
  {
    std::cerr << "WARNING: Failed to find tau+ tau- pair --> skipping the event !!\n";
    return;
  }

  bool passesAcceptanceCuts = (*acceptanceCuts_)(*tauPlus, *tauMinus);

  std::vector<const reco::GenParticle*> tauPlus_daughters;
  findDecayProducts(tauPlus, tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_ch = get_chargedHadrons(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_pi0 = get_neutralPions(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_nu = get_neutrinos(tauPlus_daughters);
  int tauPlus_decaymode = get_decayMode(tauPlus_ch, tauPlus_pi0, tauPlus_nu);
  int tauPlus_nChargedKaons = get_chargedKaons(tauPlus_daughters).size();
  int tauPlus_nNeutralKaons = get_neutralKaons(tauPlus_daughters).size();
  std::vector<const reco::GenParticle*> tauPlus_y = get_photons(tauPlus_daughters);
  int tauPlus_nPhotons = tauPlus_y.size();
  double tauPlus_sumPhotonEn = comp_visP4(tauPlus_y).energy() + tauPlus_enLoss;

  std::vector<const reco::GenParticle*> tauMinus_daughters;
  findDecayProducts(tauMinus, tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_ch = get_chargedHadrons(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_pi0 = get_neutralPions(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_nu = get_neutrinos(tauMinus_daughters);
  int tauMinus_decaymode = get_decayMode(tauMinus_ch, tauMinus_pi0, tauMinus_nu);
  int tauMinus_nChargedKaons = get_chargedKaons(tauMinus_daughters).size();
  int tauMinus_nNeutralKaons = get_neutralKaons(tauMinus_daughters).size();
  std::vector<const reco::GenParticle*> tauMinus_y = get_photons(tauMinus_daughters);
  int tauMinus_nPhotons = tauMinus_y.size();
  double tauMinus_sumPhotonEn = comp_visP4(tauMinus_y).energy() + tauMinus_enLoss;

  EntanglementNtuple* ntupleFiller = nullptr;
  if ( tauPlus_decaymode == reco::PFTau::kOneProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero )
  {
    ntupleFiller = ntupleFiller_pi_pi_;
  }
  else if ( (tauPlus_decaymode == reco::PFTau::kOneProng0PiZero   && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero  ) ||
            (tauPlus_decaymode == reco::PFTau::kOneProng1PiZero   && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero  ) )
  {
    ntupleFiller = ntupleFiller_pi_rho_;
  }
  else if ( (tauPlus_decaymode == reco::PFTau::kOneProng0PiZero   && tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero) ||
            (tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng0PiZero  ) )
  {
    ntupleFiller = ntupleFiller_pi_a1_;
  }
  else if (  tauPlus_decaymode == reco::PFTau::kOneProng1PiZero   && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero  )
  {
    ntupleFiller = ntupleFiller_rho_rho_;
  }
  else if ( (tauPlus_decaymode == reco::PFTau::kOneProng1PiZero   && tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero) ||
            (tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero && tauMinus_decaymode == reco::PFTau::kOneProng1PiZero  ) )
  {
    ntupleFiller = ntupleFiller_rho_a1_;
  }
  else if ( tauPlus_decaymode == reco::PFTau::kThreeProng0PiZero && tauMinus_decaymode == reco::PFTau::kThreeProng0PiZero )
  {
    ntupleFiller = ntupleFiller_a1_a1_;
  }
  
  if ( ntupleFiller )
  {
    ntupleFiller->fillBranches(
      evt, 
      &kineEvt_gen,
      tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons, tauPlus_sumPhotonEn,
      tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons, tauMinus_sumPhotonEn,
      &kineEvt_gen_smeared, &kineEvt_startPos_bestfit, &kineEvt_kinFit_bestfit, &kineEvt_svFit, &kineEvt_zmf,
      passesAcceptanceCuts,
      evtWeight);
    ntupleFiller_had_had_->fillBranches(
      evt,
      &kineEvt_gen,
      tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons, tauPlus_sumPhotonEn,
      tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons, tauMinus_sumPhotonEn,
      &kineEvt_gen_smeared, &kineEvt_startPos_bestfit, &kineEvt_kinFit_bestfit, &kineEvt_svFit, &kineEvt_zmf,
      passesAcceptanceCuts,
      evtWeight);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(EntanglementNtupleProducer);
