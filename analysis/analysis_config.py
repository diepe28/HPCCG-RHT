"""Name of the database file on disk to visualize, or if it does not exist
   create it and populate with injection information
"""
database = "campaign.db"

"""Boolean that specifies if we should rebuild the database on each execution
"""
rebuild_database = True

"""Path to where the LLVM log files, (*.LLVM.bin), generated by FlipIt exist.

    Notes
    ----
    analysis scripts recursilvely search this directory for the logs
"""

MAIN_PATH = "../../reliabilityExperiments/Friday27March_reliability-results/7_NOPdcDrt_arithmetic_control_runWangJV/"
LLVM_log_path =  MAIN_PATH + "runWangJV/"

"""Determines what type of log files to look for and parse.

    Notes
    ----
    Options are "Binary" or "ASCII"
"""
LLVM_log_type = "ASCII"

"""Path to where the output files of the run(s) are stored. Each output file
    should a seperate fault injection run of the application.

    Notes
    -----
        1) assumes files from a fault injection campaign end in _# or _#.txt
            where # is the trial number e.g. foo_1 or foo_1.txt
        2) stdout and stderr are in the same file

    See Also
    --------
    migrate.py to suffix _# or _#.txt to existing run output files
"""
#trial_path = "../../App-Executer/Execution-Wang/Output"
#trial_path = "../../reliabilityExperiments/results-6-August/0_MPI_NOPdcDrt_allErrorsrunNotReplicated/runNotReplicated/bck/readyFiles/output"
trial_path = LLVM_log_path + "Output"

"""Begining of the trial file name e.g. foo for trials foo_1 or foo_1.txt
"""
trial_prefix = "Output"

"""Path to source code.

    Notes
    -----
    assumes stdout and stderr are in the same file
"""
srcPath = "../"

"""Number of fault injection trials to read
"""
#numTrials = 20.
numTrials = 300.

"""Names of functions that a more detailed analysis should be conducted for

    Notes
    -----
    A more detailed analysis includes color coding source code based on
    injection frequency
"""
#more_detail_funcs = ["ddot", "ddot_producer"]
more_detail_funcs = ["ddot_producer"]

"""Snipits of system generated messages. These SHOULD be changed based on
    your system
"""
busError = "BAD TERMINATION"
assertMessage = "EXIT CODE: 11"
segError = "signal 11"

"""Snipit of a detection message. This SHOULD be changed based on your
    detection scheme
"""
detectMessage = "SOFT ERROR DETECTED"
hungMessage = "Application got stuck!"
corruptedMessage = "Corrupted execution!"

"""FlipIt Injection markers. You SHOULD NOT need to change these.
"""
bitMessage = "Bit position"
siteMessage = "/*********************************Start**************************************/"
siteEndMessage = "/*********************************End**************************************/"
