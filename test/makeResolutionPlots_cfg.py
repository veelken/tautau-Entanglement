# NB! This is not a regulard python file but a jinja2 template.
# However, the file can be turned into a functional python config file if all variables
# between double curly braces (e.g., {{ variable_name }}) are substituted with something that's
# compatible with the rest of the code.

import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring({{ inputFileNames }}),
    maxEvents_beforeCuts = cms.int32(-1),
    maxEvents_afterCuts = cms.int32(-1),
    outputEvery = cms.uint32(10000)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('{{ outputFileName }}')
)

process.makeResolutionPlots = cms.PSet(
    treeName = cms.string('ntupleProducer/{{ decayMode }}'),
    mode = cms.string('{{ mode }}'),
    collider = cms.string('{{ collider }}'),

    minVisTauPt = cms.double(-1.),
    maxAbsVisTauEta = cms.double(1.e+3),
    minTauTIP = cms.double(-1.),
    maxNumChargedKaons = cms.int32(0),
    maxNumNeutralKaons = cms.int32(0),
    maxNumPhotons = cms.int32(-1),
    maxSumPhotonEn = cms.double({{ maxSumPhotonEn }}),

    maxChi2 = cms.double(1.e+2),
    statusSelection = cms.vint32(0,1,2),
    apply_statusSelection = cms.bool(True),

    branchName_evtWeight = cms.string('evtWeight'),
    apply_evtWeight = cms.bool({{ apply_evtWeight }}),

    isDEBUG = cms.bool({{ is_debug }})
)
