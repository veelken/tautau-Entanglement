#!/usr/bin/env python3

# Example usage:
# ./test/produceNtuples.py -v 2023Nov16_a1PolVectorSignFix_wSmearing -s dy_lo_pythia8_ext -m -1 -j cluster
# For debugging:
# ./test/produceNtuples.py -v 2023Nov16_a1PolVectorSignFix_wSmearing -s dy_lo_pythia8_ext dy_lo_taudecay dy_lo_kkmc_orig -m 100 -j cluster
# ./test/produceNtuples.py -v 2024Jan26_wHiggsMassConstraint -c LHC -s ggH_htt_pythia8 -m 100 -j cluster

import datetime
import argparse
import getpass
import os
import sys

from TauAnalysis.Entanglement.tools.jobTools import getInputFileNames, build_Makefile, build_sbatchSubmission, query_yes_no, \
  build_cfg, mkdir, read_contents, save_cmd

parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--version', type = str, required = True, help = f'Version, e.g. {datetime.datetime.now().strftime("%Y%b%d")}')
parser.add_argument('-c', '--collider', type = str, choices = ['LHC', 'SuperKEKB'], default = 'SuperKEKB', help = 'Collider')
parser.add_argument('-a', '--axes', nargs = '*', type = str, choices = ['beam', 'higgs'], default = ['beam'], help = 'Axes')
parser.add_argument('-s', '--samples', nargs = '*', default = [], help = 'Whitelisted samples')
parser.add_argument('-m', '--max-input-files', type = int, default = -1, help = 'Maximum number of input files to be processed (use -1 to process all files)')
parser.add_argument('-j', '--job-type', type = str, choices = ['local', 'cluster'], required = True, help = 'Job type')
parser.add_argument('-f', '--filter', nargs = '*', choices = [ 'inv_mass', 'decay_mode', 'pt_eta' ], default = [], help = 'Event selection filters')
parser.add_argument('-w', '--verbosity', type = int, default = -1, help = 'Verbosity level')
args = parser.parse_args()

hAxes = args.axes
collider = args.collider
version = args.version
whitelist = args.samples
apply_inv_mass = 'inv_mass' in args.filter
apply_decay_mode = 'decay_mode' in args.filter
apply_pt_eta = 'pt_eta' in args.filter
run_makefile = args.job_type == 'local'
verbosity = args.verbosity

samples = None
if collider == "LHC":
    from TauAnalysis.Entanglement.samples import samples_LHC as samples
elif collider == "SuperKEKB":
    from TauAnalysis.Entanglement.samples import samples_SuperKEKB as samples
else:
    assert(False)

if not whitelist:
  run_all_samples = query_yes_no(
    "Do you really want to process all those samples: {}?".format(', '.join(samples.keys())),
    default = "no",
  )
  if not run_all_samples:
    sys.exit(0)

configDir  = os.path.join("/home",               getpass.getuser(), "Entanglement/ntuples", collider, version)
outputDir  = os.path.join("/scratch/persistent", getpass.getuser(), "Entanglement/ntuples", collider, version)
testDir    = os.path.dirname(os.path.abspath(__file__))
cmsswDir   = os.getenv('CMSSW_BASE')

produceNtuple_template = read_contents(os.path.join(testDir, "produceEntanglementNtuple_cfg.py"))

mkdir(configDir)
mkdir(outputDir)

if not save_cmd(os.path.join(configDir, "cmd.txt")):
  sys.exit(0)

jobOptions = {} # key = process, hAxis, jobId
for sampleName, sample in samples.items():
  if whitelist and sampleName not in whitelist:
    continue
  print(f"processing sample = '{sampleName}'")
  inputFilePath = sample['inputFilePath']
  print(f" inputFilePath = '{inputFilePath}'")
  inputFileNames = getInputFileNames(inputFilePath)
  numInputFiles = len(inputFileNames)
  print(f"Found {numInputFiles} input files.")
  numJobs = sample['numJobs']
  if args.max_input_files > 0 and args.max_input_files < numInputFiles:
    numJobs = int((sample['numJobs'] * args.max_input_files)/numInputFiles)
    print(f"Processing the first {args.max_input_files} input files, using {numJobs} jobs.")
    inputFileNames = inputFileNames[:args.max_input_files]
    numInputFiles = len(inputFileNames)
  is_kkmc = 'kkmc' in sampleName.lower()
  for hAxis in hAxes:
    for jobId in range(numJobs):
      idxFirstFile = int(jobId * numInputFiles / numJobs)
      idxLastFile = int((jobId + 1) * numInputFiles / numJobs - 1)
      inputFileNames_job = inputFileNames[idxFirstFile:idxLastFile + 1]
      cfgFileName_modified = os.path.join(
        configDir, f"produceEntanglementNtuple_{sampleName}_{hAxis}Axis_{jobId}_cfg.py"
      )
      rndSeed = jobId + 1
      outputFileName = f"entanglementNtuple_{sampleName}_{hAxis}Axis_{jobId}.root"
      cfg_args = {
        'inputFileNames'          : inputFileNames_job,
        'outputFileName'          : outputFileName,
        'collider'                : collider,
        'hAxis'                   : hAxis,
        'rndSeed'                 : rndSeed,
        'genWeight_includeSource' : is_kkmc,
        'apply_inv_mass'          : apply_inv_mass,
        'apply_decay_mode'        : apply_decay_mode,
        'apply_pt_eta'            : apply_pt_eta,
        'verbosity'               : verbosity,
      }
      build_cfg(produceNtuple_template, cfgFileName_modified, cfg_args)
      logFileName = cfgFileName_modified.replace("_cfg.py", ".log")
      job_key = '_'.join([sampleName, hAxis, str(jobId)])
      dependencies = [ cfgFileName_modified ]
      dependencies.extend(inputFileNames_job)
      jobOptions[job_key] = {
        'inputFileNames' : dependencies,
        'cfgFileName'    : cfgFileName_modified,
        'outputFilePath' : outputDir,
        'outputFileName' : outputFileName,
        'logFileName'    : logFileName,
        'cmd'            : 'cmsRun',
      }
if whitelist and len(jobOptions.keys()) == 0:
  raise ValueError("Invalid Configuration parameter 'whitelist' = '%s', no samples selected !!" % whitelist)

message = f"Finished building config files for {len(jobOptions)} job(s)."
if run_makefile:
  jobOptions_Makefile = []
  for job_key, job in jobOptions.items():
    commands = []
    commands.append('rm -f {}'.format(job['outputFileName']))
    commands.append('/usr/bin/time --verbose cmsRun {} &> {}'.format(job['cfgFileName'], job['logFileName']))
    commands.append('cp -v {} {}'.format(job['outputFileName'], os.path.join(outputDir, job['outputFileName'])))
    commands.append('rm -f {}'.format(job['outputFileName']))
    jobOptions_Makefile.append({
      'target'          : os.path.join(outputDir, job['outputFileName']),
      'dependencies'    : [ inputFileName.replace("file:", "") for inputFileName in job['inputFileNames'] ],
      'commands'        : commands,
      'outputFileNames' : [ os.path.join(outputDir, job['outputFileName']) ],
    })
  makeFileName = os.path.join(configDir, "Makefile")
  build_Makefile(makeFileName, jobOptions_Makefile)
  message += f" Now execute 'make -j 12 -f {makeFileName}' to start the jobs."
else:
  sbatchSubmissionFileName = os.path.join(configDir, "sbatch_submission.sh")
  build_sbatchSubmission(sbatchSubmissionFileName, jobOptions, 'produceEntanglementNtuple', version)
  os.chmod(sbatchSubmissionFileName, 0o755)
  message += f" Now execute '{sbatchSubmissionFileName}' to submit the jobs to SLURM."
print(message)
