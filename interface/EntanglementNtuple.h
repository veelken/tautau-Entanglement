#ifndef TauAnalysis_Entanglement_EntanglementNtuple_h
#define TauAnalysis_Entanglement_EntanglementNtuple_h

#include "FWCore/Framework/interface/Event.h"                     // edm::Event

#include "DataFormats/Candidate/interface/Candidate.h"            // reco::Candidate::LorentzVector

#include "TauAnalysis/Entanglement/interface/getP4_rf.h"          // getP4_rf()
#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"    // KinematicEvent
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // numMeasurements, numParameters

#include <Math/Boost.h>                                           // ROOT::Math::Boost
#include <TString.h>                                              // Form()
#include <TTree.h>                                                // TTree
#include <Rtypes.h>                                               // Float_t, Int_t

#include <vector>                                                 // std::vector<>
#include <string>                                                 // std::string

namespace
{
  template<typename T>
  void
  createBranchT(TTree* ntuple, const std::string& label, const std::string& branchName, T* address, const std::string& branchType)
  {
    std::string branchName_full = label;
    if ( branchName_full != "" ) branchName_full.append("_");
    branchName_full.append(branchName);
    std::string leaflist = Form("%s/%s", branchName_full.c_str(), branchType.c_str());
    ntuple->Branch(branchName_full.c_str(), address, leaflist.c_str());
  }

  void
  createBranchF(TTree* ntuple, const std::string& label, const std::string& branchName, Float_t* address)
  {
    createBranchT(ntuple, label, branchName, address, "F");
  }

  void
  createBranchI(TTree* ntuple, const std::string& label, const std::string& branchName, Int_t* address)
  {
    createBranchT(ntuple, label, branchName, address, "I");
  }

  double
  comp_z(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visTauP4, const ROOT::Math::Boost& boost_ttrf)
  {
    reco::Candidate::LorentzVector tauP4_ttrf = getP4_rf(tauP4, boost_ttrf);
    reco::Candidate::LorentzVector visTauP4_ttrf = getP4_rf(visTauP4, boost_ttrf);
    double z = visTauP4_ttrf.energy()/tauP4_ttrf.energy();
    return z;
  }
}

class EntanglementNtuple
{
 public:
  explicit EntanglementNtuple(TTree* ntuple);
  ~EntanglementNtuple();
    
  void
  fillBranches(const edm::Event& evt,
               const KinematicEvent* kineEvt_gen,
               Int_t tauPlus_nChargedKaons, Int_t tauPlus_nNeutralKaons, Int_t tauPlus_nPhotons, double tauPlus_sumPhotonEn,
               Int_t tauMinus_nChargedKaons, Int_t tauMinus_nNeutralKaons, Int_t tauMinus_nPhotons, double tauMinus_sumPhotonEn,
               const KinematicEvent* kineEvt_gen_smeared, 
               const KinematicEvent* kineEvt_startPos,
               const KinematicEvent* kineEvt_kinFit,
               const KinematicEvent* kineEvt_svFit,
               const KinematicEvent* kineEvt_zmf,
               bool passesAcceptanceCuts,
               double evtWeight);

 private:

  TTree* ntuple_;

  UInt_t run_;                       // run number
  UInt_t lumi_;                      // luminosity-section number
  ULong64_t event_;                  // event number  

  ULong64_t entry_;

  class branchType_KinematicEvent
  {
   public:
    branchType_KinematicEvent(const std::string& label)
      : label_(label)
      , pv_x_(0.)
      , pv_y_(0.)
      , pv_z_(0.)
      , hPlus_r_(0.)
      , hPlus_n_(0.)
      , hPlus_k_(0.)
      , tauPlus_pt_(0.)
      , tauPlus_eta_(0.)
      , tauPlus_phi_(0.)
      , tauPlus_mass_(0.)
      , tauPlus_decayMode_(-1)
      , tauPlus_tip_(0.)
      , svPlus_x_(0.)
      , svPlus_y_(0.)
      , svPlus_z_(0.)
      , visPlus_pt_(0.)
      , visPlus_eta_(0.)
      , visPlus_phi_(0.)
      , visPlus_mass_(0.)
      , nuPlus_pt_(0.)
      , nuPlus_eta_(0.)
      , nuPlus_phi_(0.)
      , nuPlus_mass_(0.)
      , zPlus_(0.)
      , hMinus_r_(0.)
      , hMinus_n_(0.)
      , hMinus_k_(0.)
      , tauMinus_pt_(0.)
      , tauMinus_eta_(0.)
      , tauMinus_phi_(0.)
      , tauMinus_mass_(0.)
      , tauMinus_decayMode_(-1)
      , tauMinus_tip_(0.)
      , svMinus_x_(0.)
      , svMinus_y_(0.)
      , svMinus_z_(0.)
      , visMinus_pt_(0.)
      , visMinus_eta_(0.)
      , visMinus_phi_(0.)
      , visMinus_mass_(0.)
      , nuMinus_pt_(0.)
      , nuMinus_eta_(0.)
      , nuMinus_phi_(0.)
      , nuMinus_mass_(0.)
      , zMinus_(0.)
      , mTauTau_(0.)
      , mVis_(0.)
      , cosThetaStar_(0.)
    {}
    ~branchType_KinematicEvent()
    {}

    void
    initBranches(TTree* ntuple)
    {
      createBranchF(ntuple, label_, "pv_x", &pv_x_);
      createBranchF(ntuple, label_, "pv_y", &pv_y_);
      createBranchF(ntuple, label_, "pv_z", &pv_z_);

      createBranchF(ntuple, label_, "recoil_pt", &recoil_pt_);
      createBranchF(ntuple, label_, "recoil_eta", &recoil_eta_);
      createBranchF(ntuple, label_, "recoil_phi", &recoil_phi_);
      createBranchF(ntuple, label_, "recoil_mass", &recoil_mass_);

      createBranchF(ntuple, label_, "hPlus_r", &hPlus_r_);
      createBranchF(ntuple, label_, "hPlus_n", &hPlus_n_);
      createBranchF(ntuple, label_, "hPlus_k", &hPlus_k_);
      createBranchF(ntuple, label_, "tauPlus_pt", &tauPlus_pt_);
      createBranchF(ntuple, label_, "tauPlus_eta", &tauPlus_eta_);
      createBranchF(ntuple, label_, "tauPlus_phi", &tauPlus_phi_);
      createBranchF(ntuple, label_, "tauPlus_mass", &tauPlus_mass_);
      createBranchI(ntuple, label_, "tauPlus_decayMode", &tauPlus_decayMode_);
      createBranchF(ntuple, label_, "tauPlus_tip", &tauPlus_tip_);
      createBranchF(ntuple, label_, "svPlus_x", &svPlus_x_);
      createBranchF(ntuple, label_, "svPlus_y", &svPlus_y_);
      createBranchF(ntuple, label_, "svPlus_z", &svPlus_z_);
      createBranchF(ntuple, label_, "visPlus_pt", &visPlus_pt_);
      createBranchF(ntuple, label_, "visPlus_eta", &visPlus_eta_);
      createBranchF(ntuple, label_, "visPlus_phi", &visPlus_phi_);
      createBranchF(ntuple, label_, "visPlus_mass", &visPlus_mass_);
      createBranchF(ntuple, label_, "nuPlus_pt", &nuPlus_pt_);
      createBranchF(ntuple, label_, "nuPlus_eta", &nuPlus_eta_);
      createBranchF(ntuple, label_, "nuPlus_phi", &nuPlus_phi_);
      createBranchF(ntuple, label_, "nuPlus_mass", &nuPlus_mass_);
      createBranchF(ntuple, label_, "zPlus", &zPlus_);

      createBranchF(ntuple, label_, "hMinus_r", &hMinus_r_);
      createBranchF(ntuple, label_, "hMinus_n", &hMinus_n_);
      createBranchF(ntuple, label_, "hMinus_k", &hMinus_k_);
      createBranchF(ntuple, label_, "tauMinus_pt", &tauMinus_pt_);
      createBranchF(ntuple, label_, "tauMinus_eta", &tauMinus_eta_);
      createBranchF(ntuple, label_, "tauMinus_phi", &tauMinus_phi_);
      createBranchF(ntuple, label_, "tauMinus_mass", &tauMinus_mass_);
      createBranchI(ntuple, label_, "tauMinus_decayMode", &tauMinus_decayMode_);
      createBranchF(ntuple, label_, "tauMinus_tip", &tauMinus_tip_);
      createBranchF(ntuple, label_, "svMinus_x", &svMinus_x_);
      createBranchF(ntuple, label_, "svMinus_y", &svMinus_y_);
      createBranchF(ntuple, label_, "svMinus_z", &svMinus_z_);
      createBranchF(ntuple, label_, "visMinus_pt", &visMinus_pt_);
      createBranchF(ntuple, label_, "visMinus_eta", &visMinus_eta_);
      createBranchF(ntuple, label_, "visMinus_phi", &visMinus_phi_);
      createBranchF(ntuple, label_, "visMinus_mass", &visMinus_mass_);
      createBranchF(ntuple, label_, "nuMinus_pt", &nuMinus_pt_);
      createBranchF(ntuple, label_, "nuMinus_eta", &nuMinus_eta_);
      createBranchF(ntuple, label_, "nuMinus_phi", &nuMinus_phi_);
      createBranchF(ntuple, label_, "nuMinus_mass", &nuMinus_mass_);
      createBranchF(ntuple, label_, "zMinus", &zMinus_);

      createBranchF(ntuple, label_, "mTauTau", &mTauTau_);
      createBranchF(ntuple, label_, "mVis", &mVis_);
      createBranchF(ntuple, label_, "cosThetaStar", &cosThetaStar_);
    }

    void
    fillBranches(const KinematicEvent& kineEvt)
    {
      const reco::Candidate::Point& pv = kineEvt.pv();
      pv_x_               = pv.x();
      pv_y_               = pv.y();
      pv_z_               = pv.z();

      const reco::Candidate::LorentzVector& recoilP4 = kineEvt.recoilP4();
      recoil_pt_          = recoilP4.pt();
      recoil_eta_         = recoilP4.eta();
      recoil_phi_         = recoilP4.phi();
      recoil_mass_        = recoilP4.mass();
      ROOT::Math::Boost boost_ttrf = ROOT::Math::Boost(recoilP4.BoostToCM());

      if ( kineEvt.hPlus_isValid() )
      {
        reco::Candidate::Vector hPlus = kineEvt.hPlus();
        hPlus_r_          = hPlus.x();
        hPlus_n_          = hPlus.y();
        hPlus_k_          = hPlus.z();
      }
      else
      {
        hPlus_r_          = 0.;
        hPlus_n_          = 0.;
        hPlus_k_          = 0.;
      }
      reco::Candidate::LorentzVector tauPlusP4;
      if ( kineEvt.tauPlusP4_isValid() )
      {
        tauPlusP4 = kineEvt.tauPlusP4();
        tauPlus_pt_       = tauPlusP4.pt();
        tauPlus_eta_      = tauPlusP4.eta();
        tauPlus_phi_      = tauPlusP4.phi();
        tauPlus_mass_     = tauPlusP4.mass();
      }
      else
      {
        tauPlus_pt_       = 0.;
        tauPlus_eta_      = 0.;
        tauPlus_phi_      = 0.;
        tauPlus_mass_     = 0.;
      }
      tauPlus_decayMode_  = kineEvt.tauPlus_decayMode();
      tauPlus_tip_        = std::sqrt((kineEvt.tipPCATauPlus() - kineEvt.pv()).mag2());
      if ( kineEvt.svTauPlus_isValid() )
      {
        const reco::Candidate::Point& svPlus = kineEvt.svTauPlus();
        svPlus_x_         = svPlus.x();
        svPlus_y_         = svPlus.y();
        svPlus_z_         = svPlus.z();
      }
      else
      {
        svPlus_x_         = 0.;
        svPlus_y_         = 0.;
        svPlus_z_         = 0.;
      }
      const reco::Candidate::LorentzVector& visPlusP4 = kineEvt.visTauPlusP4();
      visPlus_pt_         = visPlusP4.pt();
      visPlus_eta_        = visPlusP4.eta();
      visPlus_phi_        = visPlusP4.phi();
      visPlus_mass_       = visPlusP4.mass();
      if ( kineEvt.nuTauPlusP4_isValid() )
      {
        const reco::Candidate::LorentzVector& nuPlusP4 = kineEvt.nuTauPlusP4();
        nuPlus_pt_        = nuPlusP4.pt();
        nuPlus_eta_       = nuPlusP4.eta();
        nuPlus_phi_       = nuPlusP4.phi();
        nuPlus_mass_      = nuPlusP4.mass();
      }
      else
      {        
        nuPlus_pt_        = 0.;
        nuPlus_eta_       = 0.;
        nuPlus_phi_       = 0.;
        nuPlus_mass_      = 0.;
      }
      if ( kineEvt.tauPlusP4_isValid() )
      {
        zPlus_            = comp_z(tauPlusP4, visPlusP4, boost_ttrf);
      }

      if ( kineEvt.hMinus_isValid() )
      {
        reco::Candidate::Vector hMinus = kineEvt.hMinus();
        hMinus_r_         = hMinus.x();
        hMinus_n_         = hMinus.y();
        hMinus_k_         = hMinus.z();
      }
      else
      {
        hMinus_r_         = 0.;
        hMinus_n_         = 0.;
        hMinus_k_         = 0.;
      }
      reco::Candidate::LorentzVector tauMinusP4;
      if ( kineEvt.tauMinusP4_isValid() )
      {
        tauMinusP4 = kineEvt.tauMinusP4();
        tauMinus_pt_      = tauMinusP4.pt();
        tauMinus_eta_     = tauMinusP4.eta();
        tauMinus_phi_     = tauMinusP4.phi();
        tauMinus_mass_    = tauMinusP4.mass();
      }
      else
      {
        tauMinus_pt_      = 0.;
        tauMinus_eta_     = 0.;
        tauMinus_phi_     = 0.;
        tauMinus_mass_    = 0.;
      }
      tauMinus_decayMode_ = kineEvt.tauMinus_decayMode();
      tauMinus_tip_       = std::sqrt((kineEvt.tipPCATauMinus() - kineEvt.pv()).mag2());
      if ( kineEvt.svTauMinus_isValid() )
      {
        const reco::Candidate::Point& svMinus = kineEvt.svTauMinus();
        svMinus_x_        = svMinus.x();
        svMinus_y_        = svMinus.y();
        svMinus_z_        = svMinus.z();
      }
      else
      {
        svMinus_x_        = 0.;
        svMinus_y_        = 0.;
        svMinus_z_        = 0.;
      }
      const reco::Candidate::LorentzVector& visMinusP4 = kineEvt.visTauMinusP4();
      visMinus_pt_        = visMinusP4.pt();
      visMinus_eta_       = visMinusP4.eta();
      visMinus_phi_       = visMinusP4.phi();
      visMinus_mass_      = visMinusP4.mass();
      if ( kineEvt.nuTauMinusP4_isValid() )
      {
        const reco::Candidate::LorentzVector& nuMinusP4 = kineEvt.nuTauMinusP4();
        nuMinus_pt_       = nuMinusP4.pt();
        nuMinus_eta_      = nuMinusP4.eta();
        nuMinus_phi_      = nuMinusP4.phi();
        nuMinus_mass_     = nuMinusP4.mass();
      }
      else
      {        
        nuMinus_pt_       = 0.;
        nuMinus_eta_      = 0.;
        nuMinus_phi_      = 0.;
        nuMinus_mass_     = 0.;
      }
      if ( kineEvt.tauMinusP4_isValid() )
      {
        zMinus_           = comp_z(tauMinusP4, visMinusP4, boost_ttrf);
      }

      if ( kineEvt.tauPlusP4_isValid() && kineEvt.tauMinusP4_isValid() )
      {
        mTauTau_          = (tauPlusP4 + tauMinusP4).mass();
        reco::Candidate::LorentzVector tauMinusP4_ttrf = getP4_rf(tauMinusP4, boost_ttrf);
        cosThetaStar_     = cos(tauMinusP4_ttrf.theta());
      }
      else
      {
        mTauTau_          = 0.;
        cosThetaStar_     = 0.;
      }
      mVis_               = (visPlusP4 + visMinusP4).mass();
    }

    std::string label_;

    Float_t pv_x_;                   // x coordinate of primary event vertex (PV)
    Float_t pv_y_;                   // y coordinate of primary event vertex (PV)
    Float_t pv_z_;                   // z coordinate of primary event vertex (PV)

    Float_t recoil_pt_;              // transverse momentum (in laboratory frame) of tau pair, reconstructed via recoil
    Float_t recoil_eta_;             // pseudo-rapidity (in laboratory frame) of tau pair, reconstructed via recoil
    Float_t recoil_phi_;             // azimuthal angle (in laboratory frame) of tau pair, reconstructed via recoil
    Float_t recoil_mass_;            // mass component of tau pair, reconstructed via recoil

    Float_t hPlus_r_;                // r component (in helicity frame) of polarimetric vector of tau+
    Float_t hPlus_n_;                // n component (in helicity frame) of polarimetric vector of tau+
    Float_t hPlus_k_;                // k component (in helicity frame) of polarimetric vector of tau+
    Float_t tauPlus_pt_;             // transverse momentum (in laboratory frame) of tau+
    Float_t tauPlus_eta_;            // pseudo-rapidity (in laboratory frame) of tau+
    Float_t tauPlus_phi_;            // azimuthal angle (in laboratory frame) of tau+
    Float_t tauPlus_mass_;           // mass component of tau+ four-vector
    Int_t tauPlus_decayMode_;        // tau+ decay mode
    Float_t tauPlus_tip_;            // transverse impact parameter of "leading" (highest pT) track of tau+
    Float_t svPlus_x_;               // x coordinate of tau+ decay vertex (SV+)
    Float_t svPlus_y_;               // y coordinate of tau+ decay vertex (SV+)
    Float_t svPlus_z_;               // z coordinate of tau+ decay vertex (SV+)
    Float_t visPlus_pt_;             // transverse momentum (in laboratory frame) of visible decay products of tau+
    Float_t visPlus_eta_;            // pseudo-rapidity (in laboratory frame) of visible decay products of tau+
    Float_t visPlus_phi_;            // azimuthal angle (in laboratory frame) of visible decay products of tau+
    Float_t visPlus_mass_;           // mass of visible decay products of tau+
    Float_t nuPlus_pt_;              // transverse momentum (in laboratory frame) of neutrino produced in tau+ decay
    Float_t nuPlus_eta_;             // pseudo-rapidity (in laboratory frame) of neutrino produced in tau+ decay
    Float_t nuPlus_phi_;             // azimuthal angle (in laboratory frame) of neutrino produced in tau+ decay
    Float_t nuPlus_mass_;            // mass of neutrino produced in tau+ decay
    Float_t zPlus_;                  // fraction of tau+ energy (in tau-pair restframe) carried by visible decay products of tau+

    Float_t hMinus_r_;               // r component (in helicity frame) of polarimetric vector of tau-
    Float_t hMinus_n_;               // n component (in helicity frame) of polarimetric vector of tau-
    Float_t hMinus_k_;               // k component (in helicity frame) of polarimetric vector of tau-
    Float_t tauMinus_pt_;            // transverse momentum (in laboratory frame) of tau-
    Float_t tauMinus_eta_;           // pseudo-rapidity (in laboratory frame) of tau-
    Float_t tauMinus_phi_;           // azimuthal angle (in laboratory frame) of tau-
    Float_t tauMinus_mass_;          // mass component of tau- four-vector
    Int_t tauMinus_decayMode_;       // tau- decay mode
    Float_t tauMinus_tip_;           // transverse impact parameter of "leading" (highest pT) track of tau-
    Float_t svMinus_x_;              // x coordinate of tau- decay vertex (SV-)
    Float_t svMinus_y_;              // y coordinate of tau- decay vertex (SV-)
    Float_t svMinus_z_;              // z coordinate of tau- decay vertex (SV-)
    Float_t visMinus_pt_;            // transverse momentum (in laboratory frame) of visible decay products of tau-
    Float_t visMinus_eta_;           // pseudo-rapidity (in laboratory frame) of visible decay products of tau-
    Float_t visMinus_phi_;           // azimuthal angle (in laboratory frame) of visible decay products of tau-
    Float_t visMinus_mass_;          // mass of visible decay products of tau-
    Float_t nuMinus_pt_;             // transverse momentum (in laboratory frame) of neutrino produced in tau- decay
    Float_t nuMinus_eta_;            // pseudo-rapidity (in laboratory frame) of neutrino produced in tau- decay
    Float_t nuMinus_phi_;            // azimuthal angle (in laboratory frame) of neutrino produced in tau- decay
    Float_t nuMinus_mass_;           // mass of neutrino produced in tau- decay
    Float_t zMinus_;                 // fraction of tau- energy (in tau-pair restframe) carried by visible decay products of tau-

    Float_t mTauTau_;                // mass of tau pair
    Float_t mVis_;                   // mass of visible decay products of tau pair
    Float_t cosThetaStar_;           // polar angle of tau- in tau-pair restframe
  };

  branchType_KinematicEvent branches_KinematicEvent_gen_;
  Int_t tauPlus_nChargedKaons_gen_;  // number of charged Kaons produced in decay of tau-
  Int_t tauPlus_nNeutralKaons_gen_;  // number of neutral Kaons produced in decay of tau-
  Int_t tauPlus_nPhotons_gen_;       // number of photons radiated from decay products of tau+
  Float_t tauPlus_sumPhotonEn_gen_;  // energy sum of photons radiated from decay products of tau+
  Float_t tauPlus_mT_gen_;           // transverse mass of particles produced in decay of tau+
  Int_t tauMinus_nChargedKaons_gen_; // number of charged Kaons produced in decay of tau-
  Int_t tauMinus_nNeutralKaons_gen_; // number of neutral Kaons produced in decay of tau-
  Int_t tauMinus_nPhotons_gen_;      // number of photons radiated from decay products of tau-
  Float_t tauMinus_sumPhotonEn_gen_; // energy sum of photons radiated from decay products of tau-
  Float_t tauMinus_mT_gen_;          // transverse mass of particles produced in decay of tau-

  branchType_KinematicEvent branches_KinematicEvent_gen_smeared_;

  branchType_KinematicEvent branches_KinematicEvent_startPos_;
  Int_t   startPos_isCorrectSign_;   // +1 if correct solution was taken, -1 if wrong solution was taken, 0 if solution is unique

  branchType_KinematicEvent branches_KinematicEvent_kinFit_;
  Int_t   kinFit_status_;            // status of kinematic fit: +1 = converged, -1 = failed
  Float_t kinFit_cov_[kinFit::numParameters][kinFit::numParameters]; // covariance matrix V_alpha of kinematic fit
  Float_t kinFit_chi2_;              // chi^2 of kinematic fit (computed according to Eq. (4) in https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf)

  branchType_KinematicEvent branches_KinematicEvent_svFit_;
  Int_t   svFit_status_;             // status returned by ClassicSVfit algorithm: +1 = valid solution, -1 = no valid solution

  branchType_KinematicEvent branches_KinematicEvent_zmf_;

  Bool_t  passesAcceptanceCuts_;     // flag indicating whether or not all charged hadrons and all photons 
                                     // produced in the decays of tau+ and tau- pass pT/energy and eta/theta cuts

  Float_t evtWeight_;                // event weight of Monte Carlo generator
};

#endif // TauAnalysis_Entanglement_EntanglementNtuple_h
