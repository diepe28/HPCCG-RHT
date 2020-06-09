############### Injector Parameters ##################
#
#    config - config file used by the compiler pass
#    funcList - list of functions that are faulty
#    prob - probability that instuction is faulty
#    byte - which byte is faulty (0-7) -1 random
#    singleInj - one injection per active rank (0 or 1)
#    ptr - add code to inject into pointers (0 or 1)
#    arith - add code to inject into mathematics (0 or 1)
#    ctrl - add code to inject into control (0 or 1)
#    stateFile - unique counter for fault site index;
#                should differ based on application
#
#####################################################

config = "HPCCG-RHT.config"
#funcList = "\"\""
#funcList = "\"ddot ddot_producer RHT_Produce Wang_Produce VG_Produce\""
funcList = "\"ddot ddot_producer waxpby waxpby_producer HPC_sparsemv HPC_sparsemv_producer HPCCG HPCCG_producer\""
prob = 1e-7
byte = -1
bit = -1
ptr = 1
arith = 1
ctrl = 1
singleInj = 1
stateFile = "HPCCG-RHT"

############# Library Parameters #####################
#
#    FLIPIT_PATH - Path to FlipIt repo
#    SHOW - libraries and path wraped by mpicc
#
#####################################################
import os
FLIPIT_PATH = os.environ['FLIPIT_PATH']
LLVM_BUILD_PATH = os.environ['LLVM_BUILD_PATH']
# set with include and library paths from mpicc -show; or mpixx -showme
SHOW = "-Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpich"
CPP_LIB = "" # not needed for this example (C program)


########### Files to NOT inject inside ###############
notInject = [
			 "main.cpp",
			 "generate_matrix.cpp",
			 "HPC_Sparse_Matrix.cpp",
			 "read_HPC_row.cpp",
			 "compute_residual.cpp",
			 "mytimer.cpp",
			 "dump_matlab_matrix.cpp",
			 #"HPC_sparsemv.cpp",
			 #"HPCCG.cpp",
			 #"waxpby.cpp",
			 #"ddot.cpp",
			 "make_local_matrix.cpp",
			 "exchange_externals.cpp",
			 "YAML_Element.cpp",
			 "YAML_Doc.cpp",
			 "RHT.cpp",
			 "RHT.h",
			 "QueueStressTest.cpp"
			 ]

############ Default Compiler #################
#cc = "gcc"
#cc = "g++"
cc = "mpicxx"

############ Verbose compiler output ##############
verbose = True

############ Generate a histogram of fault site traversals #########
histogram = True
