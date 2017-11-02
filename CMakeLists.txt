cmake_minimum_required(VERSION 3.8)
project(HPCCG_1_0)

#set(CMAKE_BUILD_TYPE Release) # optimizations, but debuging is not allowed
set(CMAKE_BUILD_TYPE Debug) # optimizations, but debuging is not allowed

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=2 -DUSING_MPI=1 -lm")

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
        SyncQueue.cpp SyncQueue.h)
add_executable(HPCCG_1_0 ${SOURCE_FILES})

target_link_libraries(HPCCG_1_0 ${MPI_LIBRARIES})

target_link_libraries(HPCCG_1_0 pthread)

if(MPI_COMPILE_FLAGS)
    set_target_properties(HPCCG_1_0 PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(HPCCG_1_0 PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()