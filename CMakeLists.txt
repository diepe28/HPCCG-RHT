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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -lm -DUSING_MPI=1") #just for gcc

#print output to cmd
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPRINT_OUTPUT=0")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DAPPROACH_USING_POINTERS=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DAPPROACH_ALREADY_CONSUMED=1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DAPPROACH_NEW_LIMIT=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DAPPROACH_WRITE_INVERTED_NEW_LIMIT=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DAPPROACH_NEW_ENQ_INDEX=0") # use with APPROACH_NEW_LIMIT

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCOUNT_QUEUE_DESYNC=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBRANCH_HINT=1") # should always be true
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVAR_GROUPING=1")

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
add_executable(HPCCG_RHT ${SOURCE_FILES})

target_link_libraries(HPCCG_RHT ${MPI_LIBRARIES})

target_link_libraries(HPCCG_RHT pthread)

if(MPI_COMPILE_FLAGS)
    set_target_properties(HPCCG_RHT PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(HPCCG_RHT PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()