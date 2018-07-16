cmake_minimum_required(VERSION 3.7)
project(HPCCG_RHT)

set(CMAKE_BUILD_TYPE Release) # optimizations, but debuging is not allowed
#set(CMAKE_BUILD_TYPE Debug) # NO optimizations, debuging allowed

#set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_MPI=1 -lm") # for clang
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ftree-vectorize -ftree-vectorizer-verbose=7") #testing for vectorization
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --profile-generate") #--profile-generate, run tests, recompile with --profile-use
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_OMP -fopenmp")

# General compilation flags
add_definitions(-DUSING_MPI)
add_definitions(-ftree-vectorize -flto -funroll-all-loops -mavx2 -lm)
add_definitions(-DPRINT_OUTPUT=0 -DPERCENTAGE_OF_REPLICATION=0)
add_definitions(-DBRANCH_HINT=1 -DCOUNT_QUEUE_DESYNC=0 -DAPPROACH_NEW_ENQ_INDEX=0)

# These are other approaches of the literature
#add_definitions(-DAPPROACH_USING_POINTERS=1)
#add_definitions(-DAPPROACH_ALREADY_CONSUMED=1)
#add_definitions(-DAPPROACH_WANG=1)

# These are our approaches
#add_definitions(-DJUST_VOLATILES=1)
#add_definitions(-DAPPROACH_MIX_WANG=1)
#add_definitions(-DAPPROACH_MIX_IMPROVED=1)
#add_definitions(-DVAR_GROUPING=1) # it can be used with any of the above approaches

# no longer used
#add_definitions(-DSKIP_VOLATILE=1)
#add_definitions(-DAPPROACH_NEW_LIMIT=1)
#add_definitions(-DAPPROACH_WRITE_INVERTED_NEW_LIMIT=1)

#Init flags current
set(FLAGS_CURRENT "anything")

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

function(createExecutable exeName compFlags varFlag)
    add_executable(${exeName} ${SOURCE_FILES})
    #this line causes problems when used with cmake 3.7 (version of g5k machines)
    #set_target_properties(${exeName} PROPERTIES COMPILE_FLAGS "${compFlags}")
    set_property(TARGET ${exeName} APPEND PROPERTY COMPILE_DEFINITIONS ${compFlags})
    if(varFlag EQUAL 1)
        set_property(TARGET ${exeName} APPEND PROPERTY COMPILE_DEFINITIONS "VAR_GROUPING=1")
    endif()
    if(varFlag EQUAL 2)
        set_property(TARGET ${exeName} APPEND PROPERTY COMPILE_DEFINITIONS "JUST_VOLATILES=1")
    endif()
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
#createExecutable("HPCCG-current" ${FLAGS_CURRENT} 0)

# using pointer flags

set(FLAGS_USING_POINTERS "APPROACH_USING_POINTERS=1")

# using already consumed value flags
set(FLAGS_ALREADY_CONSUMED "APPROACH_ALREADY_CONSUMED=1")

# using WANG value flags
set(FLAGS_WANG "APPROACH_WANG=1")

# using no sync consumer approachh
set(FLAGS_NO_SYNC_CONSUMER "APPROACH_CONSUMER_NO_SYNC=1")

# using new limit approach (includes the previous)
set(FLAGS_NEW_LIMIT "APPROACH_NEW_LIMIT=1")

# using new limit write inverted approach
set(FLAGS_NEW_LIMIT_INVERTED "APPROACH_WRITE_INVERTED_NEW_LIMIT=1")

# using mix-wang approach
set(FLAGS_MIX_WANG "APPROACH_MIX_WANG=1")

# uncomment the next lines to create multiple executables, with different configurations
#using pointers
#createExecutable("HPCCG-UP" ${FLAGS_USING_POINTERS})

# already consumed
#createExecutable("HPCCG-AC" ${FLAGS_ALREADY_CONSUMED} 0)
#createExecutable("HPCCG-AC-VG" ${FLAGS_ALREADY_CONSUMED} 1)

# WANG approach
createExecutable("HPCCG-WANG" ${FLAGS_WANG} 0)
createExecutable("HPCCG-WANG-VG" ${FLAGS_WANG} 1)
createExecutable("HPCCG-WANG-JV" ${FLAGS_WANG} 2)

# no sync consumer
#createExecutable("HPCCG-NSC" ${FLAGS_NO_SYNC_CONSUMER})
#createExecutable("HPCCG-NSC-VG" ${FLAGS_NO_SYNC_CONSUMER_VG})

# # new limit
#createExecutable("HPCCG-NL" ${FLAGS_NEW_LIMIT} 0)
#createExecutable("HPCCG-NL-VG" ${FLAGS_NEW_LIMIT} 1)

# new limit inverted
#createExecutable("HPCCG-NLI" ${FLAGS_NEW_LIMIT_INVERTED} 0)
#createExecutable("HPCCG-NLI-VG" ${FLAGS_NEW_LIMIT_INVERTED} 1)

# # mix wang approach
#createExecutable("HPCCG-MIX" ${FLAGS_MIX_WANG} 0)
#createExecutable("HPCCG-MIX-VG" ${FLAGS_MIX_WANG} 1)

# # Best approach
#createExecutable("HPCCG-NL" ${FLAGS_NEW_LIMIT} 0)
#createExecutable("HPCCG-RHT" ${FLAGS_NEW_LIMIT} 1)
