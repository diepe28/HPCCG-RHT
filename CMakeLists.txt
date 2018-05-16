cmake_minimum_required(VERSION 3.5)
project(HPCCG_RHT)

set(CMAKE_BUILD_TYPE Release) # optimizations, but debuging is not allowed
#set(CMAKE_BUILD_TYPE Debug) # NO optimizations, debuging allowed

#set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_MPI=1 -lm") # for clang
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ftree-vectorize -ftree-vectorizer-verbose=7 -lm -DUSING_MPI=1") #testing for vectorization

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2") #--profile-generate, run tests, recompile with --profile-use

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_MPI")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_OMP -fopenmp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -lm") #just for gcc

#print output to cmd
set(FLAGS_CURRENT "${CMAKE_CXX_FLAGS} -DPRINT_OUTPUT=0 -DPERCENTAGE_OF_REPLICATION=0")
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DAPPROACH_NEW_ENQ_INDEX=0") # use with APPROACH_NEW_LIMIT
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DCOUNT_QUEUE_DESYNC=0")
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DBRANCH_HINT=1") # should always be true

# These are other approaches of the literature
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DAPPROACH_USING_POINTERS=0")
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DAPPROACH_ALREADY_CONSUMED=0 ")
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DAPPROACH_SRMT=1")

# These are our approaches
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DAPPROACH_CONSUMER_NO_SYNC=0")
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DAPPROACH_NEW_LIMIT=0") # new limit includes the consumer no sync
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DAPPROACH_WRITE_INVERTED_NEW_LIMIT=0")
set(FLAGS_CURRENT "${FLAGS_CURRENT} -DVAR_GROUPING=0") # it can be used with any of the below approaches

set(SOURCE_FILES main.cpp
        generate_matrix.cpp generate_matrix.h
        HPC_Sparse_Matrix.cpp HPC_Sparse_Matrix.h
        compute_residual.cpp compute_residual.h
        ddot.cpp ddot.h
        dump_matlab_matrix.cpp dump_matlab_matrix.h
        exchange_externals.cpp exchange_externals.h
        HPCCG.cpp HPCCG.h
        HPC_sparsemv.cpp HPC_sparsemv.h
        make_local_matrix.cpp make_local_matrix.h
        mytimer.cpp mytimer.h
        read_HPC_row.cpp read_HPC_row.h
        waxpby.cpp waxpby.h
        YAML_Doc.cpp YAML_Doc.h
        YAML_Element.cpp YAML_Element.h
        RHT.h RHT.cpp atomicops.h readerwriterqueue.h QueueStressTest.cpp QueueStressTest.h)

function(createExecutable exeName compFlags)
    add_executable(${exeName} ${SOURCE_FILES})
    set_target_properties(${exeName} PROPERTIES COMPILE_FLAGS "${compFlags}")
    target_link_libraries(${exeName} ${MPI_LIBRARIES})
    target_link_libraries(${exeName} pthread)
    target_link_libraries(${exeName} m) #math
    if(MPI_COMPILE_FLAGS)
        set_target_properties(${exeName} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    endif()

    if(MPI_LINK_FLAGS)
        set_target_properties(${exeName} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()
endfunction(createExecutable)


#current, for testing
#createExecutable("HPCCG-current" ${FLAGS_CURRENT})

# using pointer flags
set(FLAGS_USING_POINTERS "-DAPPROACH_USING_POINTERS=1")
set(FLAGS_USING_POINTERS_VG " ${FLAGS_USING_POINTERS} -DVAR_GROUPING=1")

# using already consumed value flags
set(FLAGS_ALREADY_CONSUMED "-DAPPROACH_ALREADY_CONSUMED=1")
set(FLAGS_ALREADY_CONSUMED_VG " ${FLAGS_ALREADY_CONSUMED} -DVAR_GROUPING=1")

# using SRMT value flags
set(FLAGS_SRMT "-DAPPROACH_SRMT=1")
set(FLAGS_SRMT_VG " ${FLAGS_SRMT} -DVAR_GROUPING=1")

# using no sync consumer approachh
set(FLAGS_NO_SYNC_CONSUMER "-DAPPROACH_CONSUMER_NO_SYNC=1")
set(FLAGS_NO_SYNC_CONSUMER_VG " ${FLAGS_NO_SYNC_CONSUMER} -DVAR_GROUPING=1")

# using new limit approach (includes the previous)
set(FLAGS_NEW_LIMIT "-DAPPROACH_NEW_LIMIT=1")
set(FLAGS_NEW_LIMIT_VG " ${FLAGS_NEW_LIMIT} -DVAR_GROUPING=1")

# using new limit write inverted approach
set(FLAGS_NEW_LIMIT_INVERTED "-DAPPROACH_WRITE_INVERTED_NEW_LIMIT=1")
set(FLAGS_NEW_LIMIT_INVERTED_VG " ${FLAGS_NEW_LIMIT_INVERTED} -DVAR_GROUPING=1")

# un comment the next lines to create multiple executables, with different configurations
#using pointers
#createExecutable("HPCCG-UP" ${FLAGS_USING_POINTERS})
#createExecutable("HPCCG-UP-VG" ${FLAGS_USING_POINTERS_VG})

# already consumed
createExecutable("HPCCG-AC" ${FLAGS_ALREADY_CONSUMED})
#createExecutable("HPCCG-AC-VG" ${FLAGS_ALREADY_CONSUMED_VG})

# SRMT
createExecutable("HPCCG-SRMT" ${FLAGS_SRMT})
#createExecutable("HPCCG-SRMT-VG" ${FLAGS_SRMT_VG})

# no sync consumer
#createExecutable("HPCCG-NSC" ${FLAGS_NO_SYNC_CONSUMER})
#createExecutable("HPCCG-NSC-VG" ${FLAGS_NO_SYNC_CONSUMER_VG})

# new limit
#createExecutable("HPCCG-NL" ${FLAGS_NEW_LIMIT})
#createExecutable("HPCCG-NL-VG" ${FLAGS_NEW_LIMIT_VG})

# new limit inverted
#createExecutable("HPCCG-NLI" ${FLAGS_NEW_LIMIT_INVERTED})
#createExecutable("HPCCG-NLI-VG" ${FLAGS_NEW_LIMIT_INVERTED_VG})

# Best approach
createExecutable("HPCCG-RHT" ${FLAGS_NEW_LIMIT_VG})