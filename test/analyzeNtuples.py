#!/usr/bin/env python3

import datetime
import argparse
import getpass
import sys
import os
import re

# Example usage:
# ./test/analyzeNtuples.py -v 2023Nov16_a1PolVectorSignFix_wSmearing -s dy_lo_pythia8_ext -M inclusive scan -j cluster
# ./test/analyzeNtuples.py -v 2023Nov16_a1PolVectorSignFix_wSmearing -s dy_lo_pythia8_ext -S by_differentialXsec1d by_differentialXsec2d by_asymmetry -M inclusive scan -j cluster
# ./test/analyzeNtuples.py -v 2023Nov16_a1PolVectorSignFix_wSmearing_wAcceptanceCuts -V 2023Nov16_a1PolVectorSignFix_wSmearing -s dy_lo_pythia8_ext -S by_differentialXsec1d by_differentialXsec2d by_asymmetry -M inclusive scan --apply-acceptance-cuts -j cluster
# ./test/analyzeNtuples.py -v 2023Nov16_a1PolVectorSignFix_wSmearing_wCorrectStartPosSignCut -V 2023Nov16_a1PolVectorSignFix_wSmearing -s dy_lo_pythia8_ext -S by_summation -m kinFit -M inclusive -j cluster
# For debugging:
# ./test/analyzeNtuples.py -v 2024Jan26_wHiggsMassConstraint -c LHC -s ggH_htt_pythia8 -d pi_pi pi_rho pi_a1 rho_rho rho_a1 a1_a1 had_had -M inclusive -A analyzeEntanglementNtuple makeResolutionPlots makeControlPlots -j cluster

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames, build_Makefile, query_yes_no, \
  build_cfg, mkdir, read_contents, save_cmd, positive_int_type, build_sbatchSubmission

mode_choices = [ "gen", "gen_smeared", "startPos", "kinFit", "svFit", "zmf" ]
hAxes_choices = [ "beam", "higgs" ]
collider_choices = [ "LHC", "SuperKEKB" ]
decayMode_choices = [ "pi_pi", "pi_rho", "pi_a1", "rho_rho", "rho_a1", "a1_a1" ]
decayMode_choices_all = decayMode_choices + [ "had_had" ]
spinAnalyzer_choices = [ "by_summation", "by_mlfit", "by_differentialXsec1d", "by_differentialXsec2d", "by_asymmetry" ]
analysis_modes = [ "inclusive", "scan" ]
analysis_choices = [ "analyzeEntanglementNtuple", "makeResolutionPlots", "makeControlPlots" ]
measurement_choices = [
  "zPlus", "zMinus", "cosThetaStar", "zPlus_vs_cosThetaStar", "zMinus_vs_cosThetaStar", "zPlus_vs_zMinus", "visPlusPt_vs_visMinusPt",
]

parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--version', type = str, required = True, help = f'Version, e.g. {datetime.datetime.now().strftime("%Y%b%d")}')
parser.add_argument('-V', '--ntuple-version', type = str, default = '', help = 'Version of input Ntuples (not needed if the same as --version)')
parser.add_argument('-c', '--collider', type = str, choices = collider_choices, default = 'SuperKEKB', help = 'Collider')
parser.add_argument('-C', '--apply-acceptance-cuts', dest = "apply_acceptance_cuts", action = "store_true", default = False, help = 'Apply acceptance cuts')
parser.add_argument('-a', '--axes', nargs = '*', type = str, choices = collider_choices, default = [ 'beam' ], help = 'Axes')
parser.add_argument('-S', '--spin-analyzers', nargs = '*', type = str, choices = spinAnalyzer_choices, default = spinAnalyzer_choices, help = 'Spin analyzers')
parser.add_argument('-d', '--decay-modes', nargs = '*', type = str, choices = decayMode_choices_all, default = decayMode_choices, help = 'Tau decay modes')
parser.add_argument('-m', '--modes', nargs = '*', type = str, choices = mode_choices, default = mode_choices, help = 'Input data')
parser.add_argument('-M', '--analysis-modes', nargs = '*', type = str, choices = analysis_modes, default = [ 'inclusive' ], help = 'Binned analysis')
parser.add_argument('-N', '--nbins', type = positive_int_type, default = 20, help = 'Number of bins in |cos(theta)|')
parser.add_argument('-s', '--samples', nargs = '*', default = [], help = 'Whitelisted samples')
parser.add_argument('-n', '--max-events', type = int, default = -1, help = 'Max number of events considered')
parser.add_argument('-b', '--bootstrap-size', type = int, default = -1, help = 'Size of bootstrap dataset (use -1 to consider all events from the input sample)')
parser.add_argument('-B', '--bootstrap-count', type = positive_int_type, default = 100, help = 'Number of bootstrap datasets')
parser.add_argument('-j', '--job-type', type = str, choices = ['local', 'cluster'], required = True, help = 'Job type')
parser.add_argument('-w', '--verbosity', type = int, default = 1, help = 'Verbosity level')
parser.add_argument('-A', '--analysis', nargs = '*', type = str, choices = analysis_choices, default = [ "analyzeEntanglementNtuple" ], help = 'Analysis type')
parser.add_argument('-x', '--binned-measurements', nargs = '*', type = str, choices = measurement_choices, default = [], help = "Binned measurements")
args = parser.parse_args()

version = args.version
ntuple_version = args.ntuple_version
if not ntuple_version:
  ntuple_version = version

hAxes = args.axes
collider = args.collider
modes = args.modes
analysis_modes = args.analysis_modes
nbins = args.nbins
decayModes = args.decay_modes
spinAnalyzers = args.spin_analyzers
whitelist = args.samples
max_events = args.max_events
bootstrap_size = args.bootstrap_size
bootstrap_count = args.bootstrap_count
run_makefile = args.job_type == 'local'
verbosity = args.verbosity
analysis = args.analysis
binned_measurements = args.binned_measurements
apply_acceptanceCuts = args.apply_acceptance_cuts

if 0 < max_events < bootstrap_size:
  parser.error("Max events cannot be smaller than bootstrap size if both are specified!")

absCosTheta_acceptanceCut = -1
if collider == "LHC":
    from TauAnalysis.Entanglement.samples import samples_LHC as samples, PAR_GEN_LHC as par_gen, \
      MAX_SUM_PHOTON_EN_LHC as maxSumPhotonEn
elif collider == "SuperKEKB":
    from TauAnalysis.Entanglement.samples import samples_SuperKEKB as samples, PAR_GEN_SUPERKEKB as par_gen, \
      MAX_SUM_PHOTON_EN_SUPERKEKB as maxSumPhotonEn
    absCosTheta_acceptanceCut = 0.92 # obtained with acceptanceCalculator
else:
    assert(False)

absCosTheta_bins = None
if "scan" in analysis_modes:
  absCosTheta_bins =  [ (nbins - binIdx) / nbins for binIdx in range(nbins) ]
  if absCosTheta_acceptanceCut > 0:
    # Impose cuts starting from the acceptance cut
    assert(absCosTheta_acceptanceCut < 1)
    absCosTheta_bins_new = [ binEdge for binEdge in absCosTheta_bins if binEdge < absCosTheta_acceptanceCut ]
    absCosTheta_bins_new.insert(0, absCosTheta_acceptanceCut)
    if len(absCosTheta_bins_new) < len(absCosTheta_bins):
      print(f"Switched from {len(absCosTheta_bins)} bins to {len(absCosTheta_bins_new)} bins")
    absCosTheta_bins = absCosTheta_bins_new

if not whitelist:
  run_all_samples = query_yes_no(
    "Do you really want to process all those samples: {}?".format(', '.join(samples.keys())),
    default = "no",
  )
  if not run_all_samples:
    sys.exit(0)

if not whitelist:
  run_all_samples = query_yes_no(
    "Do you really want to process all those samples: {}?".format(', '.join(samples.keys())),
    default = "no",
  )
  if not run_all_samples:
    sys.exit(0)

inputFilePath = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/ntuples", collider, ntuple_version)

configDir  = os.path.join("/home",               getpass.getuser(), "Entanglement/analysis", collider, version)
outputDir  = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/analysis", collider, version)
testDir    = os.path.dirname(os.path.abspath(__file__))
cmsswDir   = os.getenv('CMSSW_BASE')

outputDir_plots = os.path.join(configDir, "plots")

analyzeNtuple_template = read_contents(os.path.join(testDir, "analyzeEntanglementNtuple_cfg.py"))
makeControlPlot_template = read_contents(os.path.join(testDir, "makeControlPlots_cfg.py"))
makeResolutionPlot_template = read_contents(os.path.join(testDir, "makeResolutionPlots_cfg.py"))

analyzeEntanglementNtuple_cmd = 'analyzeEntanglementNtuple'
makeControlPlots_cmd = 'makeControlPlots'
makeResolutionPlots_cmd = 'makeResolutionPlots'

mkdir(configDir)
mkdir(os.path.join(configDir, "plots"))
mkdir(outputDir)

if not save_cmd(os.path.join(configDir, "cmd.txt")):
  sys.exit(0)

def init_dict(dictionary, keys):
    """Auxiliary function to initialize dictionary for access with multiple levels of keys
    """
    dictionary_at_keylevel = dictionary
    numKeys = len(keys)
    for idxKey in range(numKeys - 1):
        key = keys[idxKey]
        if key not in dictionary_at_keylevel.keys():
            dictionary_at_keylevel[key] = {}
        dictionary_at_keylevel = dictionary_at_keylevel[key]

def get_all_values(dictionary):
    """Auxiliary function to get all values stored in dictionary with multiple levels of keys
    """
    all_values = []
    for key in dictionary.keys():
        if isinstance(dictionary[key], dict):
            values = get_all_values(dictionary[key])
            all_values.extend(values)
        else:
            all_values.append(dictionary[key])
    return all_values

jobOptions_analysis              = {} # key = job_key_analysis
jobOptions_ctrlPlots             = {} # key = job_key_ctrlPlots
jobOptions_resPlots              = {} # key = job_key_resPlots

# use list of all previous writeSelBiasCorrection jobs as dependencies for next writeSelBiasCorrection job, to ensure that
# only one writeSelBiasCorrection job runs at a time; also use all writeSelBiasCorrection jobs as dependencies for all analysis jobs
outputFiles_writeSelBiasCorrection = []
for sampleName, sample in samples.items():
  if whitelist and sampleName not in whitelist:
    continue
  for hAxis in hAxes:
    print(f"processing sample = '{sampleName}', hAxis = '{hAxis}'")
    print(f" inputFilePath = '{inputFilePath}'")
    inputFile_regex = rf"entanglementNtuple_{sampleName}_{hAxis}Axis_[0-9]+.root"
    inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
    numInputFiles = len(inputFileNames)
    print("Found %i input files." % numInputFiles)
    for mode in modes:
      for decayMode in decayModes:
        if "analyzeEntanglementNtuple" in analysis and mode != "gen_smeared":
          for spinAnalyzer in spinAnalyzers:
            for analysis_mode in analysis_modes:
              absCosTheta_cuts = None
              if analysis_mode == 'inclusive':
                absCosTheta_cuts = [ -1 ]
              elif analysis_mode == 'scan':
                absCosTheta_cuts = absCosTheta_bins
              else:
                assert(False)
              for absCosTheta_cut in absCosTheta_cuts:
                if analysis_mode != "inclusive":
                  assert(absCosTheta_cut > 0)
                suffix = ""
                if analysis_mode == "scan":
                  suffix = f"absCosTheta{absCosTheta_cut:.2f}".replace(".", "p")
                cfg_baseName = f"analyzeEntanglementNtuple_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode_{spinAnalyzer}"
                if suffix:
                  cfg_baseName += f"_{suffix}"
                cfgFileName_analysis_modified = os.path.join(configDir, f"{cfg_baseName}_cfg.py")
                outputFileName_analysis = f"{cfg_baseName}.root"
                jsonOutputFileName_analysis = os.path.join(
                  configDir, re.sub(r'.root$', '.json', outputFileName_analysis)
                )
                args_analysis = {
                  'inputFileNames'       : inputFileNames,
                  'par_gen'              : par_gen,
                  'mode'                 : mode,
                  'collider'             : collider,
                  'decayMode'            : decayMode,
                  'apply_evtWeight'      : sample['apply_evtWeight'],
                  'spinAnalyzer'         : spinAnalyzer,
                  'maxSumPhotonEn'       : maxSumPhotonEn,
                  'bootstrapSize'        : bootstrap_size,
                  'numBootstrapSamples'  : bootstrap_count,
                  'outputFileName'       : outputFileName_analysis,
                  'verbosity'            : verbosity,
                  'jsonOutputFileName'   : jsonOutputFileName_analysis,
                  'max_events'           : max_events,
                  'absCosTheta_cut'      : absCosTheta_cut,
                  'binned_measurements'  : binned_measurements,
                  'apply_acceptanceCuts' : apply_acceptanceCuts
                }
                build_cfg(analyzeNtuple_template, cfgFileName_analysis_modified, args_analysis)

                logFileName_analysis = cfgFileName_analysis_modified.replace("_cfg.py", ".log")
                dependencies_analysis = [ cfgFileName_analysis_modified ]
                dependencies_analysis.extend(inputFileNames)
                sbatch_options = ''
                if decayMode == 'had_had':
                  if 0 < absCosTheta_cut < 0.3:
                    pass
                  elif 0.3 <= absCosTheta_cut < 0.8:
                    sbatch_options = '--mem 4000M'
                  else:
                    sbatch_options = '--mem 6000M'
                jobOptions_analysis[cfg_baseName] = {
                  'inputFileNames' : dependencies_analysis,
                  'cfgFileName'    : cfgFileName_analysis_modified,
                  'outputFilePath' : outputDir,
                  'outputFileName' : outputFileName_analysis,
                  'logFileName'    : logFileName_analysis,
                  'cmd'            : analyzeEntanglementNtuple_cmd,
                  'options'        : sbatch_options,
                }
        if "makeControlPlots" in analysis:
          cfgFileName_ctrlPlots_modified = os.path.join(
            configDir, f"makeControlPlots_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode_cfg.py"
          )
          outputFileName_ctrlPlots = os.path.join(
            outputDir_plots, f"makeControlPlots_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode.root"
          )
          args_ctrlPlots = {
            'inputFileNames'  : inputFileNames,
            'mode'            : mode,
            'collider'        : collider,
            'decayMode'       : decayMode,
            'apply_evtWeight' : sample['apply_evtWeight'],
            'maxSumPhotonEn'  : maxSumPhotonEn,
            'outputFileName'  : outputFileName_ctrlPlots,
            'is_debug'        : False, #verbosity >= 0,
          }
          build_cfg(makeControlPlot_template, cfgFileName_ctrlPlots_modified, args_ctrlPlots)

          logFileName_ctrlPlots = cfgFileName_ctrlPlots_modified.replace("_cfg.py", ".log")
          job_key_ctrlPlots = f'{sampleName}_{mode}_{hAxis}_{decayMode}_ctrlPlots'
          dependencies_ctrlPlots = [ cfgFileName_ctrlPlots_modified ]
          dependencies_ctrlPlots.extend(inputFileNames)
          jobOptions_ctrlPlots[job_key_ctrlPlots] = {
            'inputFileNames' : dependencies_ctrlPlots,
            'cfgFileName'    : cfgFileName_ctrlPlots_modified,
            'outputFilePath' : outputDir_plots,
            'outputFileName' : outputFileName_ctrlPlots,
            'logFileName'    : logFileName_ctrlPlots,
            'cmd'            : makeControlPlots_cmd,
          }
        if "makeResolutionPlots" in analysis and mode != "gen":
          cfgFileName_resPlots_modified = os.path.join(
            configDir, f"makeResolutionPlots_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode_cfg.py"
          )
          outputFileName_resPlots = os.path.join(
            outputDir_plots, f"makeResolutionPlots_{sampleName}_{mode}Mode_{hAxis}Axis_{decayMode}DecayMode.root"
          )
          args_resPlots = {
            'inputFileNames'  : inputFileNames,
            'mode'            : mode,
            'collider'        : collider,
            'decayMode'       : decayMode,
            'apply_evtWeight' : sample['apply_evtWeight'],
            'maxSumPhotonEn'  : maxSumPhotonEn,
            'outputFileName'  : outputFileName_resPlots,
            'is_debug'        : False, #verbosity >= 0,
          }
          build_cfg(makeResolutionPlot_template, cfgFileName_resPlots_modified, args_resPlots)

          logFileName_resPlots = cfgFileName_resPlots_modified.replace("_cfg.py", ".log")
          job_key_resPlots = f'{sampleName}_{mode}_{hAxis}_{decayMode}_resPlots'
          dependencies_resPlots = [ cfgFileName_resPlots_modified ]
          dependencies_resPlots.extend(inputFileNames)
          jobOptions_resPlots[job_key_resPlots] = {
            'inputFileNames' : dependencies_resPlots,
            'cfgFileName'    : cfgFileName_resPlots_modified,
            'outputFilePath' : outputDir,
            'outputFileName' : outputFileName_resPlots,
            'logFileName'    : logFileName_resPlots,
            'cmd'            : makeResolutionPlots_cmd,
          }

jobOptions_Makefile = []
for job_key, job in jobOptions_analysis.items():
  commands = []
  commands.append('rm -f {}'.format(job['outputFileName']))
  commands.append('/usr/bin/time --verbose {} {} &> {}'.format(analyzeEntanglementNtuple_cmd, job['cfgFileName'], job['logFileName']))
  commands.append('cp -v {} {}'.format(job['outputFileName'], os.path.join(outputDir, job['outputFileName'])))
  commands.append('rm -f {}'.format(job['outputFileName']))
  jobOptions_Makefile.append({
    'target'          : os.path.join(outputDir, job['outputFileName']),
    'dependencies'    : [ inputFileName.replace("file:", "") for inputFileName in job['inputFileNames'] ],
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, job['outputFileName']) ],
  })
for job_key, job in jobOptions_ctrlPlots.items():
  commands = []
  commands.append('rm -f {}'.format(job['outputFileName']))
  commands.append('/usr/bin/time --verbose {} {} &> {}'.format(makeControlPlots_cmd, job['cfgFileName'], job['logFileName']))
  jobOptions_Makefile.append({
    'target'          : os.path.join(outputDir, job['outputFileName']),
    'dependencies'    : [ inputFileName.replace("file:", "") for inputFileName in job['inputFileNames'] ],
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, job['outputFileName']) ],
  })
for job_key, job in jobOptions_resPlots.items():
  commands = []
  commands.append('rm -f {}'.format(job['outputFileName']))
  commands.append('/usr/bin/time --verbose {} {} &> {}'.format(makeResolutionPlots_cmd, job['cfgFileName'], job['logFileName']))
  jobOptions_Makefile.append({
    'target'          : os.path.join(outputDir, job['outputFileName']),
    'dependencies'    : [ inputFileName.replace("file:", "") for inputFileName in job['inputFileNames'] ],
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, job['outputFileName']) ],
  })

jobOptions = { **jobOptions_analysis, **jobOptions_ctrlPlots, **jobOptions_resPlots }
message = f"Finished building config files for {len(jobOptions)} job(s)."
if run_makefile:
  makeFileName = os.path.join(configDir, "Makefile")
  build_Makefile(makeFileName, jobOptions_Makefile)
  message += f" Now execute 'make -j 12 -f {makeFileName}' to start the jobs."
else:
  sbatchSubmissionFileName = os.path.join(configDir, "sbatch_submission.sh")
  build_sbatchSubmission(sbatchSubmissionFileName, jobOptions, 'analyzeEntanglementNtuple', version)
  os.chmod(sbatchSubmissionFileName, 0o755)
  message += f" Now execute '{sbatchSubmissionFileName}' to submit the jobs to SLURM."

print(message)
