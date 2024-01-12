import FWCore.ParameterSet.Config as cms

resolutions_LHC = cms.PSet(
    # settings for CLIC/FCC-ee
    #recoilResolution_px     = cms.double(1.),     # [GeV]
    #recoilResolution_py     = cms.double(1.),     # [GeV]
    #recoilResolution_pz     = cms.double(5.),     # [GeV]
    #recoilResolution_mass   = cms.double(5.),     # [GeV]
    # settings for LHC
    # N.B.: MET resolution taken from SVfit paper 
    #         Nucl.Instrum.Meth.A 862 (2017) 54 
    recoilResolution_px     = cms.double(10.),     # [GeV]
    recoilResolution_py     = cms.double(10.),     # [GeV]
    recoilResolution_pz     = cms.double(1.e+3),   # [GeV]
    recoilResolution_mass   = cms.double(1.e+3),   # [GeV]

    pvResolution_xy         = cms.double(0.0015),  # [cm]
    pvResolution_z          = cms.double(0.0030),  # [cm]
    
    trackResolution_pt      = cms.string("2.e-5*pow(pT, 2)"),                                                     # resolution on pT in units of GeV
    trackResolution_theta   = cms.double(1.e-3),   # [rad]
    trackResolution_phi     = cms.double(1.e-3),   # [rad]

    ecalResolution_energy   = cms.string("sqrt(pow(0.166*sqrt(E), 2) + pow(0.011, 2))"),                          # resolution on energy in units of GeV
    ecalResolution_theta    = cms.string("5.e-3"), # [rad]
    ecalResolution_phi      = cms.string("5.e-3"), # [rad]

    svResolution_parl       = cms.double(0.1000),  # [cm]
    svResolution_perp       = cms.double(0.0020),  # [cm]

    tipResolution_perp      = cms.double(0.0020)   # [cm]
)

resolutions_SuperKEKB = cms.PSet(
    recoilResolution_px     = cms.double(0.010),   # [GeV]
    recoilResolution_py     = cms.double(0.010),   # [GeV]
    recoilResolution_pz     = cms.double(0.100),   # [GeV]
    recoilResolution_mass   = cms.double(0.100),   # [GeV]

    pvResolution_xy         = cms.double(0.0010),  # [cm]
    pvResolution_z          = cms.double(0.0020),  # [cm]
    
    trackResolution_pt      = cms.string("sqrt(pow(1.e-3*pow(pT, 2), 2) + pow(3.e-3*pT/beta, 2))"),               # resolution on pT in units of GeV
    trackResolution_theta   = cms.double(1.e-3),   # [rad]
    trackResolution_phi     = cms.double(1.e-3),   # [rad]

    ecalResolution_energy   = cms.string("sqrt(pow(2.e-3, 2) + pow(1.6e-2*pow(E, 0.75), 2) + pow(1.2e-2*E, 2))"), # resolution on energy in units of GeV
    ecalResolution_theta    = cms.string("4.e-3/sqrt(E)"), # [rad]
    ecalResolution_phi      = cms.string("4.e-3/sqrt(E)"), # [rad]

    svResolution_parl       = cms.double(0.0500),  # [cm]
    svResolution_perp       = cms.double(0.0010),  # [cm]

    tipResolution_perp      = cms.double(0.0010)   # [cm]
)
