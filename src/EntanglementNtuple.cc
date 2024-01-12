#include "TauAnalysis/Entanglement/interface/EntanglementNtuple.h"

#include "TauAnalysis/Entanglement/interface/comp_mT.h" // comp_mT()

EntanglementNtuple::EntanglementNtuple(TTree* ntuple)
  : ntuple_(ntuple)
  , run_(0)
  , lumi_(0)
  , event_(0)
  , entry_(1)
  , branches_KinematicEvent_gen_("gen")
  , branches_KinematicEvent_gen_smeared_("gen_smeared")
  , branches_KinematicEvent_startPos_("startPos")
  , branches_KinematicEvent_kinFit_("kinFit")
  , branches_KinematicEvent_svFit_("svFit")
  , branches_KinematicEvent_zmf_("zmf")
  , passesAcceptanceCuts_(1)
  , evtWeight_(1.)
{
  ntuple_->Branch("run", &run_, "run/i");
  ntuple_->Branch("lumi", &lumi_, "lumi/i");
  ntuple_->Branch("event", &event_, "event/l");

  ntuple_->Branch("entry", &entry_, "entry/l");

  branches_KinematicEvent_gen_.initBranches(ntuple);
  createBranchI(ntuple_, "gen", "tauPlus_nChargedKaons", &tauPlus_nChargedKaons_gen_); 
  createBranchI(ntuple_, "gen", "tauPlus_nNeutralKaons", &tauPlus_nNeutralKaons_gen_);
  createBranchI(ntuple_, "gen", "tauPlus_nPhotons", &tauPlus_nPhotons_gen_); 
  createBranchF(ntuple_, "gen", "tauPlus_sumPhotonEn", &tauPlus_sumPhotonEn_gen_);
  createBranchF(ntuple_, "gen", "tauPlus_mT", &tauPlus_mT_gen_);
  createBranchI(ntuple_, "gen", "tauMinus_nChargedKaons", &tauMinus_nChargedKaons_gen_); 
  createBranchI(ntuple_, "gen", "tauMinus_nNeutralKaons", &tauMinus_nNeutralKaons_gen_);
  createBranchI(ntuple_, "gen", "tauMinus_nPhotons", &tauMinus_nPhotons_gen_); 
  createBranchF(ntuple_, "gen", "tauMinus_sumPhotonEn", &tauMinus_sumPhotonEn_gen_);
  createBranchF(ntuple_, "gen", "tauMinus_mT", &tauMinus_mT_gen_);

  branches_KinematicEvent_gen_smeared_.initBranches(ntuple);

  branches_KinematicEvent_startPos_.initBranches(ntuple);
  createBranchI(ntuple_, "startPos", "isCorrectSign", &startPos_isCorrectSign_);

  branches_KinematicEvent_kinFit_.initBranches(ntuple);
  createBranchI(ntuple_, "kinFit", "status", &kinFit_status_);
  createBranchF(ntuple_, "kinFit", "chi2", &kinFit_chi2_);
  for ( int idxRow = 0; idxRow < kinFit::numParameters; ++idxRow )
  {
    for ( int idxColumn = 0; idxColumn < kinFit::numParameters; ++idxColumn )
    {
      std::string branchName = Form("cov_r%02ic%02i", idxRow, idxColumn);
      createBranchF(ntuple_, "kinFit", branchName.c_str(), &kinFit_cov_[idxRow][idxColumn]);
    }
  }

  branches_KinematicEvent_svFit_.initBranches(ntuple);
  createBranchI(ntuple_, "svFit", "status", &svFit_status_);

  branches_KinematicEvent_zmf_.initBranches(ntuple);

  ntuple_->Branch("passesAcceptanceCuts", &passesAcceptanceCuts_, "passesAcceptanceCuts/O");

  ntuple_->Branch("evtWeight", &evtWeight_, "evtWeight/F");
}

EntanglementNtuple::~EntanglementNtuple()
{}

void
EntanglementNtuple::fillBranches(const edm::Event& evt,
                                 const KinematicEvent* kineEvt_gen, 
                                 Int_t tauPlus_nChargedKaons, Int_t tauPlus_nNeutralKaons, Int_t tauPlus_nPhotons, double tauPlus_sumPhotonEn,
                                 Int_t tauMinus_nChargedKaons, Int_t tauMinus_nNeutralKaons, Int_t tauMinus_nPhotons, double tauMinus_sumPhotonEn,
                                 const KinematicEvent* kineEvt_gen_smeared,
                                 const KinematicEvent* kineEvt_startPos,
                                 const KinematicEvent* kineEvt_kinFit,
                                 const KinematicEvent* kineEvt_svFit,
                                 const KinematicEvent* kineEvt_zmf,
                                 bool passesAcceptanceCuts,
                                 double evtWeight)
{
  run_ = evt.id().run();
  lumi_ = evt.id().luminosityBlock();
  event_ = evt.id().event();

  assert(kineEvt_gen);
  branches_KinematicEvent_gen_.fillBranches(*kineEvt_gen);
  tauPlus_nChargedKaons_gen_  = tauPlus_nChargedKaons;
  tauPlus_nNeutralKaons_gen_  = tauPlus_nNeutralKaons;
  tauPlus_nPhotons_gen_       = tauPlus_nPhotons;
  tauPlus_sumPhotonEn_gen_    = tauPlus_sumPhotonEn;
  const reco::Candidate::LorentzVector& visTauPlusP4_gen = kineEvt_gen->visTauPlusP4();
  const reco::Candidate::LorentzVector& nuTauPlusP4_gen = kineEvt_gen->nuTauPlusP4();
  tauPlus_mT_gen_             = comp_mT(
    visTauPlusP4_gen.mass(), visTauPlusP4_gen.px(), visTauPlusP4_gen.py(), 
    0., nuTauPlusP4_gen.px(), nuTauPlusP4_gen.py());
  tauMinus_nChargedKaons_gen_ = tauMinus_nChargedKaons;
  tauMinus_nNeutralKaons_gen_ = tauMinus_nNeutralKaons;
  tauMinus_nPhotons_gen_      = tauMinus_nPhotons;
  tauMinus_sumPhotonEn_gen_   = tauMinus_sumPhotonEn;
  const reco::Candidate::LorentzVector& visTauMinusP4_gen = kineEvt_gen->visTauMinusP4();
  const reco::Candidate::LorentzVector& nuTauMinusP4_gen = kineEvt_gen->nuTauMinusP4();
  tauMinus_mT_gen_            = comp_mT(
    visTauMinusP4_gen.mass(), visTauMinusP4_gen.px(), visTauMinusP4_gen.py(), 
    0., nuTauMinusP4_gen.px(), nuTauMinusP4_gen.py());

  if ( kineEvt_gen_smeared )
  {
    branches_KinematicEvent_gen_smeared_.fillBranches(*kineEvt_gen_smeared);
  }

  if ( kineEvt_startPos )
  {
    branches_KinematicEvent_startPos_.fillBranches(*kineEvt_startPos);
    startPos_isCorrectSign_ = kineEvt_startPos->startPosSign_isCorrect();
  }

  if ( kineEvt_kinFit )
  {
    branches_KinematicEvent_kinFit_.fillBranches(*kineEvt_kinFit);
    if ( kineEvt_kinFit->kinFit_isValid() )
    {
      kinFit_status_ = kineEvt_kinFit->kinFitStatus();
      kinFit_chi2_ = kineEvt_kinFit->kinFitChi2();
      const math::MatrixPxP& kinFitCov = kineEvt_kinFit->kinFitCov();
      for ( int idxRow = 0; idxRow < kinFit::numParameters; ++idxRow )
      {
        for ( int idxColumn = 0; idxColumn < kinFit::numParameters; ++idxColumn )
        {
          kinFit_cov_[idxRow][idxColumn] = kinFitCov[idxRow][idxColumn];
        }
      }
    }
    else
    {
      kinFit_status_ = -1;
      kinFit_chi2_ = -1.;
      for ( int idxRow = 0; idxRow < kinFit::numParameters; ++idxRow )
      {
        for ( int idxColumn = 0; idxColumn < kinFit::numParameters; ++idxColumn )
        {
          kinFit_cov_[idxRow][idxColumn] = 0.;
        }
      }
    } 
  }

  if ( kineEvt_svFit )
  {
    branches_KinematicEvent_svFit_.fillBranches(*kineEvt_svFit);
    if ( kineEvt_svFit->svFit_isValid() )
    {
      svFit_status_ = +1;
    }
    else
    {
      svFit_status_ = -1;
    }
  }

  if ( kineEvt_zmf )
  {
    branches_KinematicEvent_zmf_.fillBranches(*kineEvt_zmf);
  }

  passesAcceptanceCuts_ = passesAcceptanceCuts;

  evtWeight_ = evtWeight;

  ntuple_->Fill();

  ++entry_;
}
