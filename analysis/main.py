from matplotlib import pyplot as plot
from database import init, finalize
from visualize import initVis, visClassifications, visFunctions, visCrashes,\
    visSignals, visDetections, visDetectedInjections, visDetectionLatency
from custom import customInit, customParser
#import analysis_config
from analysis_config import *
import os

def visualize (c, more_detail_funcs = None):
    """Driver function to visualize a fault injection campaign.
    Parameters
    ----------
    c : object
        sqlite3 database handle that is open to a valid filled database

    more_detail_funcs : list of str
        function names to generate extra analysis of injections inside them
    Notes
    ----------
    Allows for user custom visualzations after first 'plot.show()'
    """

    try:
        os.mkdir(MAIN_PATH + CAMPAIGN_NAME + "/Plots/")
    except OSError:
        print ("Failed creating plot directory in %s " % MAIN_PATH)

    initVis(c)
    visClassifications(c, more_detail_funcs)
    visFunctions(c, more_detail_funcs)
    visCrashes(c)
    visSignals(c)
    visDetections(c, more_detail_funcs)
    visDetectedInjections(c)
    visDetectionLatency(c)
    #dperez, plot.show()

    # add custom plotting function below
    #dperez, plot.show()


if __name__ == "__main__":
    #c = init(database, LLVM_log_path, trial_path +"/"+ trial_prefix,\
    #    customFuncs=(customInit, customParser))
    ################################## hello, this line is not working, that is why the above (original) is commented
    c = init(database, LLVM_log_path, trial_path ,\
        customFuncs=(customInit, customParser))
    visualize(c, more_detail_funcs)
    finalize()
