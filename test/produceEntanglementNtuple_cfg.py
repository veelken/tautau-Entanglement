
import FWCore.ParameterSet.Config as cms

process = cms.Process("produceEntanglementNtuple")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2018_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    #input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/100000/4FC3731A-E9C4-DD47-B222-83083ECF5684.root' 
    ),
    #eventsToProcess = cms.untracked.VEventRange(
    # tau+ tau- -> pi- nu pi+ nu event for synchronization with Luca
    #    '1:97:96091' 
    # tau+ -> pi+ pi0 nu event in which Px of neutrino is 3 GeV off when running in 'rec' mode
    #    '1:97:96038'
    # tau+ -> pi+ pi0 nu event in which kinematic fit modifies svTauPlus by large amount and fails to converge
    #    '1:97:96065'
    #)
)

inputFilePath = '/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/100000/'
inputFileNames = None
processName = "qqH_htt_pythia8"
hAxis = "beam"
rndSeed = 0
outputFileName = "entanglementNtuple_%s_DEBUG.root" % processName

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##hAxis = "$hAxis"
##rndSeed = $rndSeed
##outputFileName = "$outputFileName"

inputFile_regex = r"[a-zA-Z0-9-_]+.root"

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
    print("Found %i input files." % len(inputFileNames))
    #process.source.fileNames = cms.untracked.vstring(inputFileNames)
else:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v24', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# CV: generator-level tau decay mode selection
#     
# Note that the KinematicFit will likely fail in case the event contains >= leptonic tau decay 
# and the spin analyzer vectors will be zero,
# so it is recommended to ALWAYS apply a tau decay mode selection 
# and require that BOTH taus decay hadronically and to the selected decay modes
#
process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.tauGenJets.GenParticles = cms.InputTag('prunedGenParticles')
process.analysisSequence += process.tauGenJets

process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.tauGenJetsSelectorAllHadrons.select = cms.vstring(
    'oneProng0Pi0', 
    'oneProng1Pi0', 
##    #'oneProng2Pi0', 
##    #'oneProngOther',
    'threeProng0Pi0', 
##    #'threeProng1Pi0', 
##    #'threeProngOther', 
##    #'rare'
)
process.analysisSequence += process.tauGenJetsSelectorAllHadrons

process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    cut = cms.string('pt > 20. & abs(eta) < 2.3'),
    filter = cms.bool(False)
)
process.analysisSequence += process.selectedGenHadTaus

process.selectedGenHadTauFilter = cms.EDFilter("CandViewCountFilter",
    ##src = cms.InputTag('selectedGenHadTaus'),
    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    minNumber = cms.uint32(2)
)
process.analysisSequence += process.selectedGenHadTauFilter
#--------------------------------------------------------------------------------

process.dumpGenParticles = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag('prunedGenParticles'),
    maxEventsToPrint = cms.untracked.int32(10) 
)
#process.analysisSequence += process.dumpGenParticles

process.genWeight = cms.EDProducer("GenWeightProducer",
    src = cms.InputTag('generator')
)
process.analysisSequence += process.genWeight

from TauAnalysis.Entanglement.resolutions_cfi import resolutions
from TauAnalysis.Entanglement.smearing_cfi import smearing
process.ntupleProducer = cms.EDAnalyzer("EntanglementNtupleProducer",
    src = cms.InputTag('prunedGenParticles'),
    hAxis = cms.string(hAxis),
    resolutions = resolutions,
    smearing = smearing.clone(
        rndSeed = cms.uint64(rndSeed)
    ),
    applySmearing = cms.bool(False),
    srcEvtWeights = cms.VInputTag('genWeight'),
    # CV: 0 = "regular" tau mass constraint, 1 = constraint on Gottfried-Jackson angle
    applyTauMassConstraint = cms.int32(0),
    applyLifetimeConstraint = cms.bool(False),
    verbosity = cms.untracked.int32(-1),
    #verbosity = cms.untracked.int32(1),
    cartesian = cms.untracked.bool(True)
)
process.analysisSequence += process.ntupleProducer

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(outputFileName),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.analysisSequence)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)
