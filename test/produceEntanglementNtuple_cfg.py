# NB! This is not a regulard python file but a jinja2 template.
# However, the file can be turned into a functional python config file if all variables
# between double curly braces (e.g., {{ variable_name }}) are substituted with something that's
# compatible with the rest of the code.

import FWCore.ParameterSet.Config as cms

process = cms.Process("produceEntanglementNtuple")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2018_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring({{ inputFileNames }}),
    eventsToProcess = cms.untracked.VEventRange(),
)

{% if collider == 'LHC' %}
srcGenParticles = 'prunedGenParticles'
tauPairMassCut = "mass > 120. & mass < 130."
from TauAnalysis.Entanglement.resolutions_cfi import resolutions_LHC as resolutions
from TauAnalysis.Entanglement.acceptanceCuts_cfi import acceptanceCuts_LHC as acceptanceCuts
startPosFinder_applyHiggsMassConstraint = True
svFit_applyHiggsMassConstraint = True
{% elif collider == 'SuperKEKB' %}
srcGenParticles = 'genParticles'
tauPairMassCut = "mass > 0."
from TauAnalysis.Entanglement.resolutions_cfi import resolutions_SuperKEKB as resolutions
from TauAnalysis.Entanglement.acceptanceCuts_cfi import acceptanceCuts_SuperKEKB as acceptanceCuts
startPosFinder_applyHiggsMassConstraint = False
svFit_applyHiggsMassConstraint = False
{% else %}
raise ValueError("Invalid Configuration parameter 'collider' = '{{ collider }}' !!")
{% endif %}

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun2_v2', '')

process.analysisSequence = cms.Sequence()

# process.dumpGenParticles = cms.EDAnalyzer("ParticleListDrawer",
#     src = cms.InputTag(srcGenParticles),
#     maxEventsToPrint = cms.untracked.int32(10)
# )
# process.analysisSequence += process.dumpGenParticles

#--------------------------------------------------------------------------------
# apply event selection

{% if apply_inv_mass %}
process.load("TauAnalysis.Entanglement.filterByTauPairMass_cff")
process.genTaus.src = cms.InputTag(srcGenParticles)
process.genTauPair.cut = cms.string(tauPairMassCut)
process.analysisSequence += process.filterByTauPairMass
{% endif %}

{% if apply_decay_mode %}
process.load("TauAnalysis.Entanglement.filterByTauDecayMode_cff")
process.analysisSequence += process.filterByTauDecayMode
{% endif %}

{% if apply_pt_eta %}
{% if not apply_decay_mode %}
process.load("TauAnalysis.Entanglement.filterByTauDecayMode_cff")
process.analysisSequence += process.tauGenJets
process.analysisSequence += process.tauGenJetsSelectorAllHadrons
{% endif %}
process.load("TauAnalysis.Entanglement.filterByVisPtAndEta_cff")
process.analysisSequence += process.filterByVisPtAndEta
{% endif %}
#--------------------------------------------------------------------------------

process.genWeight = cms.EDProducer("GenWeightProducer",
    src = cms.InputTag({%- if genWeight_includeSource -%}'source', {%- endif -%}'generator')
)
process.analysisSequence += process.genWeight
 
process.load("TauAnalysis.Entanglement.EntanglementNtupleProducer_cfi")
process.ntupleProducer.src = cms.InputTag(srcGenParticles)
process.ntupleProducer.collider = cms.string('{{ collider }}')
process.ntupleProducer.hAxis = cms.string('{{ hAxis }}')
process.ntupleProducer.resolutions = resolutions
process.ntupleProducer.smearing.rndSeed = cms.uint64({{ rndSeed }})
process.ntupleProducer.startPosFinder.applyHiggsMassConstraint = cms.bool(startPosFinder_applyHiggsMassConstraint)
process.ntupleProducer.startPosFinder.skip = cms.bool(False)
process.ntupleProducer.kinematicFit.skip = cms.bool(False)
process.ntupleProducer.svFit.applyHiggsMassConstraint = cms.bool(svFit_applyHiggsMassConstraint)
process.ntupleProducer.svFit.skip = cms.bool(False)
process.ntupleProducer.acceptanceCuts = acceptanceCuts
process.ntupleProducer.verbosity = cms.untracked.int32({{ verbosity }})
process.analysisSequence += process.ntupleProducer

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string('{{ outputFileName }}'),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.analysisSequence)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
