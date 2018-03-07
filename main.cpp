
//@HEADER
// ************************************************************************
// 
//               HPCCG: Simple Conjugate Gradient Benchmark Code
//                 Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

// Changelog
//
// Version 0.3
// - Added timing of setup time for sparse MV
// - Corrected percentages reported for sparse MV with overhead
//
/////////////////////////////////////////////////////////////////////////

// Main routine of a program that reads a sparse matrix, right side
// vector, solution vector and initial guess from a file  in HPC
// format.  This program then calls the HPCCG conjugate gradient
// solver to solve the problem, and then prints results.

// Calling sequence:

// test_HPCCG linear_system_file

// Routines called:

// read_HPC_row - Reads in linear system

// mytimer - Timing routine (compile with -DWALL to get wall clock
//           times

// HPCCG - CG Solver

// compute_residual - Compares HPCCG solution to known solution.

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include <string>
#include <cmath>
#ifdef USING_MPI
#include <mpi.h> // If this routine is compiled with -DUSING_MPI
                 // then include mpi.h
#include "make_local_matrix.h" // Also include this function
#include "readerwriterqueue.h"
#endif
#ifdef USING_OMP
#include <omp.h>
#endif
#include "generate_matrix.h"
#include "read_HPC_row.h"
#include "mytimer.h"
#include "HPC_sparsemv.h"
#include "compute_residual.h"
#include "HPCCG.h"
#include "HPC_Sparse_Matrix.h"
#include "dump_matlab_matrix.h"

#include "YAML_Element.h"
#include "YAML_Doc.h"
#include "QueueStressTest.h"

#undef DEBUG

using namespace moodycamel;

#define NUM_RUNS 300

//-D CMAKE_C_COMPILER=/usr/bin/clang-5.0 -D CMAKE_CXX_COMPILER=/usr/bin/clang++-5.0
//-D CMAKE_C_COMPILER=/usr/bin/gcc-7 -D CMAKE_CXX_COMPILER=/usr/bin/g++-7

// a trick to remove ALREADY_CONSUMED VALUE, have another field indicating the times it has been read, something like that
// try instead of ASM(pause), a inner loop
//batching, try to update the deqPtr locally, and when it reaches a threshold we update the shared variable, backoff

typedef struct {
    HPC_Sparse_Matrix *A;
    double *x, *b;
    int max_iter, niters, executionCore;
    double tolerance, normr;
    double *times;
} ConsumerParams;

void
PrintSummary(const HPC_Sparse_Matrix *A, const double *times, int nx, int ny, int nz, int size, int rank, int niters,
             double normr, double t4min, double t4max, double t4avg);

void freeMemory(HPC_Sparse_Matrix *sparseMatrix, double * x, double * b, double * xexact);

void consumer_thread_func(void * args);

#define mySize 999999
#define myTimes 100

void testVectorization() throw(){
    v4sf v1, v2 ,v3;
    float n1[mySize], n2[mySize];
    double t1, t2;

    int i, k, m;
    float vector1[mySize], vector2[mySize];


    for (i = 0; i < mySize; i++) {
        vector1[i] = vector2[i] = i;
    }

    t1 = mytimer();

    while(k++ < myTimes) {
        for (i = 0; i < mySize; i++) {
            n1[i] = vector1[i] * 2;
            n2[i] = vector2[i] * 2;

            if (n1[i] != n2[i]) {
                printf("SoftError...\n");
                exit(0);
            }
        }
    }
    k = 0;
    t1 = mytimer() - t1;

    printf("Normal comparison: %f seconds. \n", t1);

    t1 = mytimer();
    while(k++ < myTimes) {
        for (i = 0; i < mySize; i+=4) {
            printf("v1 %f vs v2 %f\n", vector1[i], vector2[i]);

            n1[i] = vector1[i] * 2;
            n1[i+1] = vector1[i+1] * 2;
            n1[i+2] = vector1[i+2] * 2;
            n1[i+3] = vector1[i+3] * 2;

            n2[i] = vector2[i] * 2;
            n2[i+1] = vector2[i+1] * 2;
            n2[i+2] = vector2[i+2] * 2;
            n2[i+3] = vector2[i+3] * 2;

            if (n1[i] != n2[i] ||
                n1[i+1] != n2[i+1] ||
                n1[i+2] != n2[i+2] ||
                n1[i+3] != n2[i+3]){
                printf("Soft Error\n");
                exit(0);
            }


        }
    }

    t1 = mytimer() - t1;
    printf("Vectorized comparison: %f seconds. \n", t1);

    exit(0);
}

int main(int argc, char *argv[]) {
    //testVectorization();

    HPC_Sparse_Matrix *sparseMatrix;

    double *x, *b, *xexact, *x2, *b2, *xexact2;
    double norm, d;
    int ierr = 0;
    int i, j, iterator;
    int ione = 1;
    double times[7];
    double t6 = 0.0;
    int nx, ny, nz;
    int replicated, numCores, producerCore, consumerCore, *coreNumbers;
    double t4 = times[4];
    double t4min = 0.0;
    double t4max = 0.0;
    double t4avg = 0.0;
    double timesBaseline[NUM_RUNS], meanBaseline, sdBaseline, producerMean, consumerMean;
    double timesRHT[NUM_RUNS], meanRHT, sdRHT;

#ifdef USING_MPI

    MPI_Init(&argc, &argv);
    int size, rank; // Number of MPI processes, My process ID
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //  if (size < 100) cout << "Process "<<rank<<" of "<<size<<" is alive." <<endl;

#else

    int size = 1; // Serial case (not using MPI)
    int rank = 0;

#endif

#ifdef DEBUG
    if (rank==0)
     {
      int junk = 0;
      cout << "Press enter to continue"<< endl;
      cin >> junk;
     }

    MPI_Barrier(MPI_COMM_WORLD);
#endif

//    TestQueues();
//    return 0;

    if (argc < 3) { // dperez, original argc != 4
        if (rank == 0)
            cerr << "Usage:" << endl
                 << "Mode 1: " << argv[0] << " nx ny nz" << endl
                 << "     where nx, ny and nz are the local sub-block dimensions, or" << endl
                 << "Mode 2: " << argv[0] << " HPC_data_file " << endl
                 << "     where HPC_data_file is a globally accessible file containing matrix data." << endl;
        exit(1);
    }

    if (argc == 4) { // original number
        nx = atoi(argv[1]);
        ny = atoi(argv[2]);
        nz = atoi(argv[3]);
        replicated = 0;
    } else {
        if (argc > 4) { // dperez, for our purposes
            replicated = 1;
            nx = atoi(argv[1]);
            ny = atoi(argv[2]);
            nz = atoi(argv[3]);
            numCores = atoi(argv[4]); // should be the same as -np of MPI
            coreNumbers = (int *) malloc(sizeof(int) * numCores * 2);
            for (int i = 0; i < numCores * 2; i++) {
                // Each pair is the producer and consumer core for each MPI process.
                // Example 2 0 2 1 3: means 2 MPI processes the first one runs on core 0,2 and
                // the second one is cores 1,3. Depending on the machine config it may be HT or not.
                coreNumbers[i] = atoi(argv[5 + i]);
            }
        } else {
            read_HPC_row(argv[1], &sparseMatrix, &x, &b, &xexact);
        }
    }

    generate_matrix(nx, ny, nz, &sparseMatrix, &x, &b, &xexact);

    bool dump_matrix = false;
    if (dump_matrix && size <= 4) dump_matlab_matrix(sparseMatrix, rank);

#ifdef USING_MPI
    // Transform matrix indices from global to local values.
    // Define number of columns for the local matrix.

    t6 = mytimer();
    make_local_matrix(sparseMatrix);
    t6 = mytimer() - t6;
    times[6] = t6;

#endif

    double t1 = mytimer();   // Initialize it (if needed)
    int niters = 0;
    double normr = 0.0;
    int max_iter = 150;
    double tolerance = 0.0; // Set tolerance to zero to make all runs do max_iter iterations

    if (replicated) {
        // --- Protected runs ---
        ConsumerParams *consumerParams = (ConsumerParams *) (malloc(sizeof(ConsumerParams)));
        int local_nrow = nx * ny * nz; // This is the size of our subblock

        x2 = new double[local_nrow];

        consumerParams->max_iter = max_iter;
        consumerParams->tolerance = tolerance;
        consumerParams->niters = niters;
        consumerParams->normr = normr;
        consumerParams->times = times;

        // Assigning core numbers to each thread of each MPI process
        producerCore = coreNumbers[rank * 2];
        consumerCore = coreNumbers[rank * 2 + 1];

        consumerParams->executionCore = consumerCore;

        printf("\n--- REPLICATED VERSION ON RANK %d WITH CORES %d, %d\n", rank, producerCore, consumerCore);

        for (iterator = meanRHT = 0; iterator < NUM_RUNS; iterator++) {
            RHT_Replication_Init(1);

            // Parameters that need to be reset every run
            for (int m = 0; m < local_nrow; m++) {
                x2[m] = x[m];
            }
            consumerParams->A = sparseMatrix;
            consumerParams->b = b;
            consumerParams->x = x2;

            int err = pthread_create(consumerThreads[0], NULL, (void *(*)(void *)) consumer_thread_func,
                                     (void *) consumerParams);
            if (err) {
                fprintf(stderr, "Fail creating thread %d\n", 1);
                exit(1);
            }

            SetThreadAffinity(producerCore);
#if APPROACH_NEW_LIMIT == 1 || APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
            ierr = HPCCG_producer_newLimit(sparseMatrix, b, x, max_iter, tolerance, niters, normr, times);
#else
            ierr = HPCCG_producer(sparseMatrix, b, x, max_iter, tolerance, niters, normr, times);
#endif

            /*-- RHT -- */ pthread_join(*consumerThreads[0], NULL);

            timesRHT[iterator] = times[0];
            meanRHT += times[0];
            consumerMean += consumerCount;
            producerMean += producerCount;

            RHT_Replication_Finish();

            if (rank == 0) {
                printf("RHT approach: ");

#if APPROACH_USING_POINTERS == 1
                printf("USING POINTERS");
#elif APPROACH_ALREADY_CONSUMED == 1
                printf("ALREADY CONSUMED");
#elif APPROACH_NEW_LIMIT == 1
                printf("NEW LIMIT & ALREADY CONSUMED");
#elif APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
                printf("NEW LIMIT_WRITE_INVERTED");
#endif
                printf(" [%d]: %f seconds, on cores: %d, %d --- ProducerWaiting: %ld, ConsumerWaiting: %ld\n",
                       iterator, timesRHT[iterator], producerCore, consumerCore, producerCount, consumerCount);
            }
            freeMemory(sparseMatrix, x, b, xexact);
            generate_matrix(nx, ny, nz, &sparseMatrix, &x, &b, &xexact);
#ifdef USING_MPI
            // Transform matrix indices from global to local values.
            // Define number of columns for the local matrix.
            t6 = mytimer();
            make_local_matrix(sparseMatrix);
            t6 = mytimer() - t6;
            times[6] = t6;
#endif
        }

        meanRHT /= NUM_RUNS;
        consumerMean /= NUM_RUNS;
        producerMean /= NUM_RUNS;

        for (iterator = sdRHT = 0; iterator < NUM_RUNS; iterator++) {
            sdRHT += fabs(meanRHT - timesRHT[iterator]);
        }

        sdRHT /= NUM_RUNS;
        delete x2;
    } else {
        // Not replicated (normal) --- Unprotected runs ---
        for (iterator = meanBaseline = 0; iterator < NUM_RUNS; iterator++) {
            ierr = HPCCG(sparseMatrix, b, x, max_iter, tolerance, niters, normr, times);
            timesBaseline[iterator] = times[0];
            meanBaseline += times[0];

            freeMemory(sparseMatrix, x, b, xexact);
            generate_matrix(nx, ny, nz, &sparseMatrix, &x, &b, &xexact);

#ifdef USING_MPI
            // Transform matrix indices from global to local values.
            // Define number of columns for the local matrix.
            t6 = mytimer();
            make_local_matrix(sparseMatrix);
            t6 = mytimer() - t6;
            times[6] = t6;
#endif
            if (rank == 0) {
                //PrintSummary(sparseMatrix, times, nx, ny, nz, size, rank, niters, normr, t4min, t4max, t4avg);
                printf("Baseline[%d]: %f seconds --- \n", iterator, timesBaseline[iterator]);
            }
        }

        meanBaseline /= NUM_RUNS;

        for (iterator = sdBaseline = 0; iterator < NUM_RUNS; iterator++) {
            sdBaseline += fabs(meanBaseline - timesBaseline[iterator]);
        }

        sdBaseline /= NUM_RUNS;
    }

    if (ierr) cerr << "Error in call to CG: " << ierr << ".\n" << endl;

#ifdef USING_MPI
//    t4 = times[4];
//    t4min = 0.0;
    t4max = 0.0;
    t4avg = 0.0;

    MPI_Allreduce(&t4, &t4min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&t4, &t4max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&t4, &t4avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t4avg = t4avg / ((double) size);
#endif

// Finish up
#ifdef USING_MPI
    MPI_Finalize();
#endif

    if (rank == 0) {
        printf("\n-------------------------- Summary --------------------------\n");
        if (replicated == 0) {
            printf("Mean baseline %f , SD baseline %f \n\n", meanBaseline, sdBaseline);
        } else {
            printf("Mean RHT approach: ");
#if APPROACH_USING_POINTERS == 1
            printf("USING POINTERS");
#elif APPROACH_ALREADY_CONSUMED == 1
            printf("ALREADY CONSUMED");
#elif APPROACH_NEW_LIMIT == 1
            printf("NEW LIMIT & ALREADY CONSUMED");
#elif APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
            printf("NEW LIMIT_WRITE_INVERTED");
#endif

            printf(": %f , SD RHT %f --- PWaiting: %lf, CWaiting: %lf \n\n", meanRHT, sdRHT, producerMean, consumerMean);
        }
    }

    return 0;
}

void PrintSummary(const HPC_Sparse_Matrix *A, const double *times, int nx, int ny, int nz, int size, int rank, int niters,
             double normr, double t4min, double t4max, double t4avg) {// initialize YAML doc

    if (rank == 0)  // Only PE 0 needs to compute and report timing results
    {
        double fniters = niters;
        double fnrow = A->total_nrow;
        double fnnz = A->total_nnz;
        double fnops_ddot = fniters * 4 * fnrow;
        double fnops_waxpby = fniters * 6 * fnrow;
        double fnops_sparsemv = fniters * 2 * fnnz;
        double fnops = fnops_ddot + fnops_waxpby + fnops_sparsemv;

        YAML_Doc doc("hpccg", "1.0");

        doc.add("Parallelism", "");

#ifdef USING_MPI
        doc.get("Parallelism")->add("Number of MPI ranks",size);
#else
        doc.get("Parallelism")->add("MPI not enabled", "");
#endif

#ifdef USING_OMP
        int nthreads = 1;
#pragma omp parallel
        nthreads = omp_get_num_threads();
        doc.get("Parallelism")->add("Number of OpenMP threads",nthreads);
#else
        doc.get("Parallelism")->add("OpenMP not enabled", "");
#endif

        doc.add("Dimensions", "");
        doc.get("Dimensions")->add("nx", nx);
        doc.get("Dimensions")->add("ny", ny);
        doc.get("Dimensions")->add("nz", nz);


        doc.add("Number of iterations: ", niters);
        doc.add("Final residual: ", normr);
        doc.add("********** Performance Summary (times in sec) ***********", "");

        doc.add("Time Summary", "");
        doc.get("Time Summary")->add("Total   ", times[0]);
        doc.get("Time Summary")->add("DDOT    ", times[1]);
        doc.get("Time Summary")->add("WAXPBY  ", times[2]);
        doc.get("Time Summary")->add("SPARSEMV", times[3]);

        doc.add("FLOPS Summary", "");
        doc.get("FLOPS Summary")->add("Total   ", fnops);
        doc.get("FLOPS Summary")->add("DDOT    ", fnops_ddot);
        doc.get("FLOPS Summary")->add("WAXPBY  ", fnops_waxpby);
        doc.get("FLOPS Summary")->add("SPARSEMV", fnops_sparsemv);

        doc.add("MFLOPS Summary", "");
        doc.get("MFLOPS Summary")->add("Total   ", fnops / times[0] / 1.0E6);
        doc.get("MFLOPS Summary")->add("DDOT    ", fnops_ddot / times[1] / 1.0E6);
        doc.get("MFLOPS Summary")->add("WAXPBY  ", fnops_waxpby / times[2] / 1.0E6);
        doc.get("MFLOPS Summary")->add("SPARSEMV", fnops_sparsemv / (times[3]) / 1.0E6);

#ifdef USING_MPI
        doc.add("DDOT Timing Variations","");
        doc.get("DDOT Timing Variations")->add("Min DDOT MPI_Allreduce time",t4min);
        doc.get("DDOT Timing Variations")->add("Max DDOT MPI_Allreduce time",t4max);
        doc.get("DDOT Timing Variations")->add("Avg DDOT MPI_Allreduce time",t4avg);

        double totalSparseMVTime = times[3] + times[5]+ times[6];
        doc.add("SPARSEMV OVERHEADS","");
        doc.get("SPARSEMV OVERHEADS")->add("SPARSEMV MFLOPS W OVERHEAD",fnops_sparsemv/(totalSparseMVTime)/1.0E6);
        doc.get("SPARSEMV OVERHEADS")->add("SPARSEMV PARALLEL OVERHEAD Time", (times[5]+times[6]));
        doc.get("SPARSEMV OVERHEADS")->add("SPARSEMV PARALLEL OVERHEAD Pct", (times[5]+times[6])/totalSparseMVTime*100.0);
        doc.get("SPARSEMV OVERHEADS")->add("SPARSEMV PARALLEL OVERHEAD Setup Time", (times[6]));
        doc.get("SPARSEMV OVERHEADS")->add("SPARSEMV PARALLEL OVERHEAD Setup Pct", (times[6])/totalSparseMVTime*100.0);
        doc.get("SPARSEMV OVERHEADS")->add("SPARSEMV PARALLEL OVERHEAD Bdry Exch Time", (times[5]));
        doc.get("SPARSEMV OVERHEADS")->add("SPARSEMV PARALLEL OVERHEAD Bdry Exch Pct", (times[5])/totalSparseMVTime*100.0);
#endif

        if (rank == 0) { // only PE 0 needs to compute and report timing results
            std::__cxx11::string yaml = doc.generateYAML();
            cout << yaml;
        }
    }

    // Compute difference between known exact solution and computed solution
    // All processors are needed here.

    double residual = 0;
    //  if ((ierr = compute_residual(A->local_nrow, x, xexact, &residual)))
    //  cerr << "Error in call to compute_residual: " << ierr << ".\n" << endl;

    // if (rank==0)
    //   cout << "Difference between computed and exact  = "
    //        << residual << ".\n" << endl;
}

void freeMemory(HPC_Sparse_Matrix *sparseMatrix, double * x, double * b, double * xexact){
    destroyMatrix(sparseMatrix);
    delete x;
    delete b;
    delete xexact;
}

void consumer_thread_func(void * args) {
    ConsumerParams *params = (ConsumerParams *) args;

    SetThreadAffinity(params->executionCore);

    HPCCG_consumer(params->A, params->b, params->x, params->max_iter,
                   params->tolerance, params->niters, params->normr, params->times);
}