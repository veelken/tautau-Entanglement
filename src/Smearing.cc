#include "TauAnalysis/Entanglement/interface/Smearing.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // Candidate::LorentzVector, Candidate::Point, Candidate::Vector

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/comp_visP4.h"                // comp_visP4()
#include "TauAnalysis/Entanglement/interface/constants.h"                 // kLHC, kSuperKEKB, mChargedPion, mTau
#include "TauAnalysis/Entanglement/interface/fixMass.h"                   // fixTauMass()
#include "TauAnalysis/Entanglement/interface/get_leadTrack.h"             // get_leadTrack()
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/KinematicParticle.h"         // KinematicParticle
#include "TauAnalysis/Entanglement/interface/rotateVector.h"              // rotateVector()
#include "TauAnalysis/Entanglement/interface/square.h"                    // square()

#include <cmath>                                                          // std::sqrt()
#include <iostream>                                                       // std::cout

Smearing::Smearing(const edm::ParameterSet& cfg)
  : resolutions_(nullptr)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<Smearing::Smearing>:\n";
  }

  edm::ParameterSet cfg_smearing = cfg.getParameterSet("smearing");
  applySmearing_recoil_px_     = cfg_smearing.getParameter<bool>("applySmearing_recoil_px");
  applySmearing_recoil_py_     = cfg_smearing.getParameter<bool>("applySmearing_recoil_py");
  applySmearing_recoil_pz_     = cfg_smearing.getParameter<bool>("applySmearing_recoil_pz");
  applySmearing_recoil_mass_   = cfg_smearing.getParameter<bool>("applySmearing_recoil_mass");

  applySmearing_pv_xy_         = cfg_smearing.getParameter<bool>("applySmearing_pv_xy");
  applySmearing_pv_z_          = cfg_smearing.getParameter<bool>("applySmearing_pv_z");

  applySmearing_track_pt_      = cfg_smearing.getParameter<bool>("applySmearing_track_pt");
  applySmearing_track_theta_   = cfg_smearing.getParameter<bool>("applySmearing_track_theta");
  applySmearing_track_phi_     = cfg_smearing.getParameter<bool>("applySmearing_track_phi");

  applySmearing_ecal_energy_   = cfg_smearing.getParameter<bool>("applySmearing_ecal_energy");
  applySmearing_ecal_theta_    = cfg_smearing.getParameter<bool>("applySmearing_ecal_theta");
  applySmearing_ecal_phi_      = cfg_smearing.getParameter<bool>("applySmearing_ecal_phi");

  applySmearing_sv_perp_       = cfg_smearing.getParameter<bool>("applySmearing_sv_perp");
  applySmearing_sv_parl_       = cfg_smearing.getParameter<bool>("applySmearing_sv_parl");

  applySmearing_tip_perp_      = cfg_smearing.getParameter<bool>("applySmearing_tip_perp");

  rndSeed_                     = cfg_smearing.getParameter<unsigned long long>("rndSeed");
  rnd_.SetSeed(rndSeed_);

  if ( verbosity_ >= 1 )
  {
    std::cout << "applySmearing_recoil:\n";
    std::cout << "Recoil:" 
              << " Px = " << applySmearing_recoil_px_ << ","
              << " Py = " << applySmearing_recoil_py_ << ","
              << " Pz = " << applySmearing_recoil_pz_ << ","
              << " mass = " << applySmearing_recoil_mass_ << "\n";
    std::cout << "PV:" 
              << " xy = " << applySmearing_pv_xy_ << ","
              << " z = " << applySmearing_pv_z_ << "\n";
    std::cout << "Track:" 
              << " pT = " << applySmearing_track_pt_ << ","
              << " theta = " << applySmearing_track_theta_ << ","
              << " phi = " << applySmearing_track_phi_ << "\n";
    std::cout << "ECAL:" 
              << " E = " << applySmearing_ecal_energy_ << ","
              << " theta = " << applySmearing_ecal_theta_ << ","
              << " phi = " << applySmearing_ecal_phi_ << "\n";
    std::cout << "SV:" 
              << " perp = " << applySmearing_sv_perp_ << ","
              << " parl = " << applySmearing_sv_parl_ << "\n";
    std::cout << "TIP:" 
              << " perp = " << applySmearing_tip_perp_ << "\n";
    std::cout << "rndSeed = " << rndSeed_ << "\n";
  }

  edm::ParameterSet cfg_resolutions = cfg.getParameterSet("resolutions");
  resolutions_ = new Resolutions(cfg_resolutions);

  std::string collider = cfg.getParameter<std::string>("collider");
  if      ( collider == "LHC"       ) collider_ = kLHC;
  else if ( collider == "SuperKEKB" ) collider_ = kSuperKEKB;
  else throw cmsException("Smearing", __LINE__)
    << "Invalid Configuration parameter 'collider' = " << collider << " !!\n";
}

Smearing::~Smearing()
{
  delete resolutions_;
}

KinematicEvent
Smearing::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<Smearing::operator()>:\n";
    std::cout << " rndSeed = " << rnd_.GetSeed() << "\n";
  }

  KinematicEvent kineEvt_smeared(kineEvt);

  kineEvt_smeared.pv_ = smear_pv(kineEvt.pv());

  kineEvt_smeared.recoilP4_ = smear_recoil_p4(kineEvt.recoilP4());
  double drecoilPx = kineEvt_smeared.recoilP4_.px()     - kineEvt.recoilP4().px();
  double drecoilPy = kineEvt_smeared.recoilP4_.py()     - kineEvt.recoilP4().py();
  double drecoilPz = kineEvt_smeared.recoilP4_.pz()     - kineEvt.recoilP4().pz();
  double drecoilE  = kineEvt_smeared.recoilP4_.energy() - kineEvt.recoilP4().energy();

  double u = 0.5;
  if ( kineEvt.tauPlusP4_isValid() && kineEvt.tauMinusP4_isValid() )
  {
    u = kineEvt.tauPlusP4().energy()/(kineEvt.tauPlusP4().energy() + kineEvt.tauMinusP4().energy());
  }
  else
  {
    u = rnd_.Uniform(0., 1.);
  }

  if ( kineEvt.tauPlusP4_isValid() )
  {
    const reco::Candidate::LorentzVector& tauPlusP4 = kineEvt.tauPlusP4();
    double tauPlusPx_smeared = tauPlusP4.px()      + u*drecoilPx;
    double tauPlusPy_smeared = tauPlusP4.py()      + u*drecoilPy;
    double tauPlusPz_smeared = tauPlusP4.pz()      + u*drecoilPz;
    double tauPlusE_smeared  = tauPlusP4.energy()  + u*drecoilE;
    reco::Candidate::LorentzVector tauPlusP4_smeared(tauPlusPx_smeared, tauPlusPy_smeared, tauPlusPz_smeared, tauPlusE_smeared);
    kineEvt_smeared.tauPlusP4_ = fixTauMass(tauPlusP4_smeared);
  }
  if ( kineEvt.svTauPlus_isValid() )
  {
    kineEvt_smeared.svTauPlus_ = smear_sv(kineEvt.visTauPlusP4(), kineEvt.svTauPlus());
    if ( verbosity_ >= 3 )
    {
      auto svTauPlus_dxyz = kineEvt_smeared.svTauPlus_ - kineEvt.svTauPlus();
      auto svTauPlus_drnk = rotateVector(svTauPlus_dxyz, kineEvt.tauPlus_rotMatrix_xyz2rnk());
      std::cout << "svTauPlus: dx = " << svTauPlus_dxyz.x() << ", dy = " << svTauPlus_dxyz.y() << ", dz = " << svTauPlus_dxyz.z() << "\n";
      std::cout << " (dr = " << svTauPlus_drnk.x() << ", dn = " << svTauPlus_drnk.y() << ", dk = " << svTauPlus_drnk.z() << ")\n";
    }
  }
  const std::vector<KinematicParticle>& daughtersTauPlus = kineEvt.daughtersTauPlus();
  kineEvt_smeared.daughtersTauPlus_.clear();
  for ( const KinematicParticle& daughter : daughtersTauPlus )
  {
    KinematicParticle daughter_smeared = daughter;
    daughter_smeared.p4_ = smear_daughter_p4(daughter);
    daughter_smeared.vertex_ = kineEvt_smeared.svTauPlus_;
    kineEvt_smeared.daughtersTauPlus_.push_back(daughter_smeared);
  }
  kineEvt_smeared.visTauPlusP4_ = comp_visP4(kineEvt_smeared.daughtersTauPlus_);
  kineEvt_smeared.tipPCATauPlus_ = smear_tipPCA(daughtersTauPlus, kineEvt.tipPCATauPlus());

  if ( kineEvt.tauMinusP4_isValid() )
  {
    const reco::Candidate::LorentzVector& tauMinusP4 = kineEvt.tauMinusP4();
    double tauMinusPx_smeared = tauMinusP4.px()      + (1. - u)*drecoilPx;
    double tauMinusPy_smeared = tauMinusP4.py()      + (1. - u)*drecoilPy;
    double tauMinusPz_smeared = tauMinusP4.pz()      + (1. - u)*drecoilPz;
    double tauMinusE_smeared  = tauMinusP4.energy()  + (1. - u)*drecoilE;
    reco::Candidate::LorentzVector tauMinusP4_smeared(tauMinusPx_smeared, tauMinusPy_smeared, tauMinusPz_smeared, tauMinusE_smeared);
    kineEvt_smeared.tauMinusP4_ = fixTauMass(tauMinusP4_smeared);
  }
  if ( kineEvt.svTauMinus_isValid() )
  {
    kineEvt_smeared.svTauMinus_ = smear_sv(kineEvt.visTauMinusP4(), kineEvt.svTauMinus());
    if ( verbosity_ >= 3 )
    {
      auto svTauMinus_dxyz = kineEvt_smeared.svTauMinus_ - kineEvt.svTauMinus();
      auto svTauMinus_drnk = rotateVector(svTauMinus_dxyz, kineEvt.tauMinus_rotMatrix_xyz2rnk());
      std::cout << "svTauMinus: dx = " << svTauMinus_dxyz.x() << ", dy = " << svTauMinus_dxyz.y() << ", dz = " << svTauMinus_dxyz.z() << "\n";
      std::cout << " (dr = " << svTauMinus_drnk.x() << ", dn = " << svTauMinus_drnk.y() << ", dk = " << svTauMinus_drnk.z() << ")\n";
    }
  }
  const std::vector<KinematicParticle>& daughtersTauMinus = kineEvt.daughtersTauMinus();
  kineEvt_smeared.daughtersTauMinus_.clear();
  for ( const KinematicParticle& daughter : daughtersTauMinus )
  {
    KinematicParticle daughter_smeared = daughter;
    daughter_smeared.p4_ = smear_daughter_p4(daughter);
    daughter_smeared.vertex_ = kineEvt_smeared.svTauPlus_;
    kineEvt_smeared.daughtersTauMinus_.push_back(daughter_smeared);
  }
  kineEvt_smeared.visTauMinusP4_ = comp_visP4(kineEvt_smeared.daughtersTauMinus_);
  kineEvt_smeared.tipPCATauMinus_ = smear_tipPCA(daughtersTauMinus, kineEvt.tipPCATauMinus());

  kineEvt_smeared.hPlus_ = reco::Candidate::Vector(0.,0.,0.);
  kineEvt_smeared.hPlus_isValid_ = false;
  kineEvt_smeared.hMinus_ = reco::Candidate::Vector(0.,0.,0.);
  kineEvt_smeared.hMinus_isValid_ = false;

  return kineEvt_smeared;
}

reco::Candidate::Point
Smearing::smear_pv(const reco::Candidate::Point& pv)
{
  double smeared_pvX = pv.x();
  double smeared_pvY = pv.y();
  double smeared_pvZ = pv.z();
  if ( applySmearing_pv_xy_ )
  {
    smeared_pvX += rnd_.Gaus(0., resolutions_->pvResolution_xy());
    smeared_pvY += rnd_.Gaus(0., resolutions_->pvResolution_xy());
  }
  if ( applySmearing_pv_z_ )
  {
    smeared_pvZ += rnd_.Gaus(0., resolutions_->pvResolution_z());
  }
  return reco::Candidate::Point(smeared_pvX, smeared_pvY, smeared_pvZ);
}

reco::Candidate::LorentzVector
Smearing::smear_recoil_p4(const reco::Candidate::LorentzVector& recoilP4)
{
  double smeared_recoilPx = recoilP4.px();
  double smeared_recoilPy = recoilP4.py();
  double smeared_recoilPz = recoilP4.pz();
  double smeared_recoilM  = recoilP4.mass();
  if ( applySmearing_recoil_px_ )
  {
    smeared_recoilPx += rnd_.Gaus(0., resolutions_->recoilResolution_px());
  }
  if ( applySmearing_recoil_py_ )
  {
    smeared_recoilPy += rnd_.Gaus(0., resolutions_->recoilResolution_py());
  }
  if ( applySmearing_recoil_pz_ )
  {
    smeared_recoilPz += rnd_.Gaus(0., resolutions_->recoilResolution_pz());
  }
  if ( applySmearing_recoil_mass_ )
  {
    smeared_recoilM  += rnd_.Gaus(0., resolutions_->recoilResolution_mass());
    //double min_recoilM = 0.;
    double min_recoilM = 2.*mTau;
    if ( smeared_recoilM < min_recoilM ) smeared_recoilM = min_recoilM;
  }
  double smeared_recoilE  = std::sqrt(square(smeared_recoilPx) + square(smeared_recoilPy) + square(smeared_recoilPz) + square(smeared_recoilM));
  return reco::Candidate::LorentzVector(smeared_recoilPx, smeared_recoilPy, smeared_recoilPz, smeared_recoilE);
}

reco::Candidate::Point
Smearing::smear_sv(const reco::Candidate::LorentzVector& p4, const reco::Candidate::Point& sv)
{
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(p4, nullptr, nullptr, kBeam, collider_, r, n, k);
  double dr = rnd_.Gaus(0., resolutions_->svResolution_perp());
  double dn = rnd_.Gaus(0., resolutions_->svResolution_perp());
  double dk = rnd_.Gaus(0., resolutions_->svResolution_parl());
  double smeared_svX = sv.x();
  double smeared_svY = sv.y();
  double smeared_svZ = sv.z();
  if ( applySmearing_sv_perp_ )
  {
    smeared_svX += dr*r.x() + dn*n.x();
    smeared_svY += dr*r.y() + dn*n.y();
    smeared_svZ += dr*r.z() + dn*n.z();
  }
  if ( applySmearing_sv_parl_ )
  {
    smeared_svX += dk*k.x();
    smeared_svY += dk*k.y();
    smeared_svZ += dk*k.z();
  }
  reco::Candidate::Point smeared_sv(smeared_svX, smeared_svY, smeared_svZ);
  return smeared_sv;
}

reco::Candidate::LorentzVector
Smearing::smear_daughter_p4(const KinematicParticle& daughter)
{
  double sigma_pt    = 0.;
  double sigma_theta = 0.;
  double sigma_phi   = 0.;      
  if ( std::fabs(daughter.charge()) > 0.5 )
  {
    sigma_pt         = get_trackResolution_pt(daughter.p4(), *resolutions_);
    sigma_theta      = resolutions_->trackResolution_theta();
    sigma_phi        = resolutions_->trackResolution_phi();
  }
  else if ( daughter.pdgId() == 111 )
  {
    sigma_pt         = get_ecalResolution_pt(daughter.p4(), *resolutions_);
    sigma_theta      = get_ecalResolution_theta(daughter.p4(), *resolutions_);
    sigma_phi        = get_ecalResolution_phi(daughter.p4(), *resolutions_);
  }
  const reco::Candidate::LorentzVector& daughterP4 = daughter.p4();
  double smeared_daughterPt    = rnd_.Gaus(daughterP4.pt(), sigma_pt);
  double smeared_daughterTheta = rnd_.Gaus(daughterP4.theta(), sigma_theta);
  double smeared_daughterPhi   = rnd_.Gaus(daughterP4.phi(), sigma_phi);
  double smeared_daughterPx    = smeared_daughterPt*cos(smeared_daughterPhi);
  double smeared_daughterPy    = smeared_daughterPt*sin(smeared_daughterPhi);
  double smeared_daughterPz    = smeared_daughterPt/tan(smeared_daughterTheta);
  double smeared_daughterE     = std::sqrt(square(smeared_daughterPx) + square(smeared_daughterPy) + square(smeared_daughterPz) + square(daughterP4.mass()));
  return reco::Candidate::LorentzVector(smeared_daughterPx, smeared_daughterPy, smeared_daughterPz, smeared_daughterE);
}

reco::Candidate::Point
Smearing::smear_tipPCA(const std::vector<KinematicParticle>& daughters, const reco::Candidate::Point& tipPCA)
{
  const KinematicParticle* leadTrack = get_leadTrack(daughters);
  if ( !leadTrack )
  {
    std::cerr << "WARNING: Failed to find leading track of tau --> returning null vector !!\n";
    return reco::Candidate::Point(0.,0.,0.);
  }
  reco::Candidate::Vector r, n, k;
  get_localCoordinateSystem(leadTrack->p4(), nullptr, nullptr, kBeam, collider_, r, n, k);
  double dr = rnd_.Gaus(0., resolutions_->tipResolution_perp());
  double dn = rnd_.Gaus(0., resolutions_->tipResolution_perp());
  double smeared_pcaX = tipPCA.x();
  double smeared_pcaY = tipPCA.y();
  double smeared_pcaZ = tipPCA.z();
  if ( applySmearing_tip_perp_ )
  {
    smeared_pcaX += dr*r.x() + dn*n.x();
    smeared_pcaY += dr*r.y() + dn*n.y();
    smeared_pcaZ += dr*r.z() + dn*n.z();
  }
  reco::Candidate::Point smeared_tipPCA(smeared_pcaX, smeared_pcaY, smeared_pcaZ);
  return smeared_tipPCA;
}

