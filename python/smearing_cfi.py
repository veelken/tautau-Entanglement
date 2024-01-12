import FWCore.ParameterSet.Config as cms

smearing = cms.PSet(
    # settings for SuperKEKB, CLIC/FCC-ee
    #applySmearing_recoil_px   = cms.bool(True),
    #applySmearing_recoil_py   = cms.bool(True),
    #applySmearing_recoil_pz   = cms.bool(True),
    #applySmearing_recoil_mass = cms.bool(True),
    # settings for LHC 
    applySmearing_recoil_px   = cms.bool(True),
    applySmearing_recoil_py   = cms.bool(True),
    applySmearing_recoil_pz   = cms.bool(False),
    applySmearing_recoil_mass = cms.bool(False),

    applySmearing_pv_xy       = cms.bool(True),
    applySmearing_pv_z        = cms.bool(True),

    applySmearing_track_pt    = cms.bool(False),
    applySmearing_track_theta = cms.bool(False),
    applySmearing_track_phi   = cms.bool(False),

    applySmearing_ecal_energy = cms.bool(False),
    applySmearing_ecal_theta  = cms.bool(False),
    applySmearing_ecal_phi    = cms.bool(False),

    applySmearing_sv_perp     = cms.bool(True),
    applySmearing_sv_parl     = cms.bool(True),

    applySmearing_tip_perp    = cms.bool(True),

    rndSeed                   = cms.uint64(0)
)
