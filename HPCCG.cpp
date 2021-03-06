
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

/////////////////////////////////////////////////////////////////////////

// Routine to compute an approximate solution to Ax = b where:

// A - known matrix stored as an HPC_Sparse_Matrix struct

// b - known right hand side vector

// x - On entry is initial guess, on exit new approximate solution

// max_iter - Maximum number of iterations to perform, even if
//            tolerance is not met.

// tolerance - Stop and assert convergence if norm of residual is <=
//             to tolerance.

// niters - On output, the number of iterations actually performed.

/////////////////////////////////////////////////////////////////////////

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <cmath>
#include "mytimer.h"
#include "HPCCG.h"
#include "RHT.h"

#define TICK()  clock_gettime(CLOCK_MONOTONIC, &newStart); // Use TICK and TOCK to time a code section
#define TOCK(t)                                                 \
    clock_gettime(CLOCK_MONOTONIC, &newEnd);                    \
    t += (newEnd.tv_sec - newStart.tv_sec) + (newEnd.tv_nsec - newStart.tv_nsec) / 1000000000.0;

int HPCCG(HPC_Sparse_Matrix * hpc_sparse_matrix,
	  const double * const b, double * const x,
	  const int max_iter, const double tolerance, int &niters, double & normr,
	  double * times) {

    struct timespec startAll, newStart, newEnd;
    clock_gettime(CLOCK_MONOTONIC, &startAll);
    //double t_begin = mytimer();  // Start timing right away

    double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;

#ifdef USING_MPI
    double t5 = 0.0;
#endif
    int nrow = hpc_sparse_matrix->local_nrow;
    int ncol = hpc_sparse_matrix->local_ncol;

    double *r = new double[nrow];
    double *p = new double[ncol]; // In parallel case, hpc_sparse_matrix is rectangular
    double *Ap = new double[nrow];

    normr = 0.0;
    double rtrans = 0.0;
    double oldrtrans = 0.0;

#ifdef USING_MPI
    int rank; // Number of MPI processes, My process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    int rank = 0; // Serial case (not using MPI)
#endif

    int print_freq = max_iter / 10;
    if (print_freq > 50) print_freq = 50;
    if (print_freq < 1) print_freq = 1;

    // p is of length ncols, copy x to p for sparse MV operation
    TICK();
    waxpby(nrow, 1.0, x, 0.0, x, p);
    TOCK(t2);
#ifdef USING_MPI
    TICK();
    exchange_externals(hpc_sparse_matrix, p);
    TOCK(t5);
#endif
    TICK();
    HPC_sparsemv(hpc_sparse_matrix, p, Ap);
    TOCK(t3);
    TICK();
    waxpby(nrow, 1.0, b, -1.0, Ap, r);
    TOCK(t2);
    TICK();
    ddot(nrow, r, r, &rtrans, t4);
    TOCK(t1);
    normr = sqrt(rtrans);

#if PRINT_OUTPUT == 1
    if (rank == 0) cout << "Initial Residual = " << normr << endl;
#endif
    for (int k = 1; k < max_iter && normr > tolerance; k++) {
        if (k == 1) {
            TICK();
            waxpby(nrow, 1.0, r, 0.0, r, p);
            TOCK(t2);
        } else {
            oldrtrans = rtrans;
            TICK();
            ddot(nrow, r, r, &rtrans, t4);
            TOCK(t1);// 2*nrow ops
            double beta = rtrans / oldrtrans;
            TICK();
            waxpby(nrow, 1.0, r, beta, p, p);
            TOCK(t2);// 2*nrow ops
        }
        normr = sqrt(rtrans);
#if PRINT_OUTPUT == 1
        if (rank == 0 && (k % print_freq == 0 || k + 1 == max_iter))
            cout << "Iteration = " << k << "   Residual = " << normr << endl;
#endif

#ifdef USING_MPI
        TICK();
        exchange_externals(hpc_sparse_matrix, p);
        TOCK(t5);
#endif
        TICK();
        HPC_sparsemv(hpc_sparse_matrix, p, Ap);
        TOCK(t3); // 2*nnz ops
        double alpha = 0.0;
        TICK();
        ddot(nrow, p, Ap, &alpha, t4);
        TOCK(t1); // 2*nrow ops
        alpha = rtrans / alpha;
        TICK();
        waxpby(nrow, 1.0, x, alpha, p, x);// 2*nrow ops
        waxpby(nrow, 1.0, r, -alpha, Ap, r);
        TOCK(t2);// 2*nrow ops
        niters = k;
    }

    // Store times
    times[1] = t1; // ddot time
    times[2] = t2; // waxpby time
    times[3] = t3; // sparsemv time
    times[4] = t4; // AllReduce time
#ifdef USING_MPI
    times[5] = t5; // exchange boundary time
#endif
    delete[] p;
    delete[] Ap;
    delete[] r;
    // Total time. All done...
    GetTimeSince(startAll, times[0] =) //times[0] = mytimer() - t_begin;
    return (0);
}

int HPCCG_producer(HPC_Sparse_Matrix *hpc_sparse_matrix,
                            const double *const b, double *const x,
                            const int max_iter, const double tolerance, int &niters, double &normr,
                            double *times) {
    struct timespec startAll, newStart, newEnd;
    clock_gettime(CLOCK_MONOTONIC, &startAll);
    //double t_begin = mytimer();  // Start timing right away

    double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
    FLIPIT_SetInjector(FLIPIT_OFF);
    /*-- RHT -- */ RHT_Produce(t0);
    /*-- RHT -- */ RHT_Produce(t1);
    /*-- RHT -- */ RHT_Produce(t2);
    /*-- RHT -- */ RHT_Produce(t3);
    /*-- RHT -- */ RHT_Produce(t4);
    FLIPIT_SetInjector(FLIPIT_ON);

#ifdef USING_MPI
    double t5 = 0.0;
    FLIPIT_SetInjector(FLIPIT_OFF);
    /*-- RHT -- */ RHT_Produce(t5);
    FLIPIT_SetInjector(FLIPIT_ON);
#endif

    int nrow = hpc_sparse_matrix->local_nrow;
    int ncol = hpc_sparse_matrix->local_ncol;

    FLIPIT_SetInjector(FLIPIT_OFF);
#if JUST_VOLATILES == 1
    /*-- RHT -- */ RHT_Produce_Volatile(nrow);
#else
    /*-- RHT -- */ RHT_Produce(nrow);
#endif
    /*-- RHT -- */ RHT_Produce_Volatile(ncol);
    FLIPIT_SetInjector(FLIPIT_ON);

    double *r = new double[nrow];
    double *p = new double[ncol]; // In parallel case, hpc_sparse_matrix is rectangular
    double *Ap = new double[nrow];

    normr = 0.0;
    double rtrans = 0.0;
    double oldrtrans = 0.0;
    FLIPIT_SetInjector(FLIPIT_OFF);
    /*-- RHT -- */ RHT_Produce(normr);
    /*-- RHT -- */ RHT_Produce(rtrans);
    /*-- RHT -- */ RHT_Produce(oldrtrans);
    FLIPIT_SetInjector(FLIPIT_ON);

#ifdef USING_MPI
    int rank; // Number of MPI processes, My process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    FLIPIT_SetInjector(FLIPIT_OFF);
		/*-- RHT -- */ RHT_Produce_NoCheck(rank);
    FLIPIT_SetInjector(FLIPIT_ON);

#else
    int rank = 0; // Serial case (not using MPI)
#endif

    int print_freq = max_iter / 10;
    FLIPIT_SetInjector(FLIPIT_OFF);
    /*-- RHT -- */ RHT_Produce(print_freq);
    FLIPIT_SetInjector(FLIPIT_ON);

    if (print_freq > 50) print_freq = 50;
    if (print_freq < 1) print_freq = 1;

    FLIPIT_SetInjector(FLIPIT_OFF);
    /*-- RHT -- */ RHT_Produce(print_freq);
    FLIPIT_SetInjector(FLIPIT_ON);
    // p is of length ncols, copy x to p for sparse MV operation
    TICK();
    waxpby_producer(nrow, 1.0, x, 0.0, x, p);
    TOCK(t2);

//    printf("Producer so far so good\n");

#ifdef USING_MPI
    TICK();
    exchange_externals_producer(hpc_sparse_matrix, p);
    TOCK(t5);
#endif

    TICK();
    HPC_sparsemv_producer(hpc_sparse_matrix, p, Ap);
    TOCK(t3);

    TICK();
    waxpby_producer(nrow, 1.0, b, -1.0, Ap, r);
    TOCK(t2);

    TICK();
    ddot_producer(nrow, r, r, &rtrans, t4);
    TOCK(t1);

    normr = sqrt(rtrans);

    FLIPIT_SetInjector(FLIPIT_OFF);
    /*-- RHT -- */ RHT_Produce(rank);
    /*-- RHT Volatile -- */ RHT_Produce_Volatile(normr);
    FLIPIT_SetInjector(FLIPIT_ON);
#if PRINT_OUTPUT == 1
    if (rank == 0) cout << "Initial Residual = " << normr << endl;
#endif

    TICK();
    waxpby_producer(nrow, 1.0, r, 0.0, r, p);
    TOCK(t2);

    double beta = 0;
    int k = 1;

    goto inFor;
    for (; k < max_iter && normr > tolerance; k++) {
        oldrtrans = rtrans;
        FLIPIT_SetInjector(FLIPIT_OFF);
        /*-- RHT -- */ RHT_Produce(print_freq);
        FLIPIT_SetInjector(FLIPIT_ON);

        TICK();
        ddot_producer(nrow, r, r, &rtrans, t4);
        TOCK(t1);// 2*nrow ops

        beta = rtrans / oldrtrans;
        FLIPIT_SetInjector(FLIPIT_OFF);
        /*-- RHT -- */ RHT_Produce(beta);
        FLIPIT_SetInjector(FLIPIT_ON);

        TICK();
        waxpby_producer(nrow, 1.0, r, beta, p, p);
        TOCK(t2);// 2*nrow ops

    inFor:

        normr = sqrt(rtrans);

        FLIPIT_SetInjector(FLIPIT_OFF);
#if JUST_VOLATILES == 1
        /*-- RHT -- */ RHT_Produce_Volatile(k);
        /*-- RHT -- */ RHT_Produce_Volatile(rank);
        /*-- RHT -- */ RHT_Produce_Volatile(print_freq);
        /*-- RHT -- */ RHT_Produce_Volatile(max_iter);
#else
        /*-- RHT -- */ RHT_Produce(k);
        /*-- RHT -- */ RHT_Produce(rank);
        /*-- RHT -- */ RHT_Produce(print_freq);
        /*-- RHT -- */ RHT_Produce(max_iter);
#endif

        /*-- RHT Volatile -- */ RHT_Produce_Volatile(normr);
        FLIPIT_SetInjector(FLIPIT_ON);

#if PRINT_OUTPUT == 1
        if (rank == 0 && (k % print_freq == 0 || k + 1 == max_iter))
            cout << "Iteration = " << k << "   Residual = " << normr << endl;
#endif


#ifdef USING_MPI
        TICK();
        exchange_externals_producer(hpc_sparse_matrix, p);
        TOCK(t5);
#endif
        TICK();
        HPC_sparsemv_producer(hpc_sparse_matrix, p, Ap);
        TOCK(t3); // 2*nnz ops
        double alpha = 0.0;
        FLIPIT_SetInjector(FLIPIT_OFF);
        /*-- RHT -- */ RHT_Produce(alpha);
        FLIPIT_SetInjector(FLIPIT_ON);

        TICK();
        ddot_producer(nrow, p, Ap, &alpha, t4);
        TOCK(t1); // 2*nrow ops

        alpha = rtrans / alpha;
        FLIPIT_SetInjector(FLIPIT_OFF);
        /*-- RHT -- */ RHT_Produce(alpha);
        FLIPIT_SetInjector(FLIPIT_ON);

        TICK();
        waxpby_producer(nrow, 1.0, x, alpha, p, x);// 2*nrow ops
        waxpby_producer(nrow, 1.0, r, -alpha, Ap, r);
        TOCK(t2);// 2*nrow ops

        niters = k;
        FLIPIT_SetInjector(FLIPIT_OFF);
        /*-- RHT -- */ RHT_Produce(niters);
        FLIPIT_SetInjector(FLIPIT_ON);
    }

    #if APPROACH_WANG == 1
        // done replication but UNIT might not have been reached
        wangQueue.enqPtr = wangQueue.enqPtrLocal;
    #endif

    /// TODO, what to do with times? should we exchange them, I mean it is not necessary and since we are doing this
    /// manually we can decide what is worth replicating or not...
    // Store times
    times[1] = t1; // ddot time
    times[2] = t2; // waxpby time
    times[3] = t3; // sparsemv time
    times[4] = t4; // AllReduce time:
#ifdef USING_MPI
    times[5] = t5; // exchange boundary time
#endif
    delete[] p;
    delete[] Ap;
    delete[] r;

    // Total time. All done...
    GetTimeSince(startAll, times[0] = ) //times[0] = mytimer() - t_begin;

    return (0);
}

int HPCCG_consumer(HPC_Sparse_Matrix * hpc_sparse_matrix,
                   const double * const b, double * const x,
                   const int max_iter, const double tolerance, int &niters, double & normr,
                   double * times) {
    //-- RHT -- Not replicated double t_begin = mytimer();  // Start timing right away

    double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
    /*-- RHT -- */ RHT_Consume_Check(t0);
    /*-- RHT -- */ RHT_Consume_Check(t1);
    /*-- RHT -- */ RHT_Consume_Check(t2);
    /*-- RHT -- */ RHT_Consume_Check(t3);
    /*-- RHT -- */ RHT_Consume_Check(t4);

#ifdef USING_MPI
    double t5 = 0.0;
    /*-- RHT -- */ RHT_Consume_Check(t5);
#endif

    int nrow = hpc_sparse_matrix->local_nrow;
    int ncol = hpc_sparse_matrix->local_ncol;

#if JUST_VOLATILES == 1
    /*-- RHT -- */ RHT_Consume_Volatile(nrow);
#else
    /*-- RHT -- */ RHT_Consume_Check(nrow);
#endif
    /*-- RHT -- */ RHT_Consume_Volatile(ncol);

    double *r = new double[nrow];
    double *p = new double[ncol]; // In parallel case, hpc_sparse_matrix is rectangular
    double *Ap = new double[nrow];

    normr = 0.0;
    double rtrans = 0.0;
    double oldrtrans = 0.0;
    /*-- RHT -- */ RHT_Consume_Check(normr);
    /*-- RHT -- */ RHT_Consume_Check(rtrans);
    /*-- RHT -- */ RHT_Consume_Check(oldrtrans);

#ifdef USING_MPI
    int rank; // Number of MPI processes, My process ID
    /*-- RHT -- */ rank = (int) RHT_Consume();
#else
    int rank = 0; // Serial case (not using MPI)
#endif

    int print_freq = max_iter / 10;
    /*-- RHT -- */ RHT_Consume_Check(print_freq);

    if (print_freq > 50) print_freq = 50;
    if (print_freq < 1) print_freq = 1;

    /*-- RHT -- */ RHT_Consume_Check(print_freq);

    // p is of length ncols, copy x to p for sparse MV operation
    //-- RHT Not replicated -- TICK();
    waxpby_consumer(nrow, 1.0, x, 0.0, x, p);
    //-- RHT Not replicated --TOCK(t2);

//    printf("Consumer so far so good\n");

#ifdef USING_MPI
    //-- RHT Not replicated --TICK();
    exchange_externals_consumer(hpc_sparse_matrix,p);
    //-- RHT Not replicated --TOCK(t5);
#endif

    //-- RHT Not replicated --TICK();
    HPC_sparsemv_consumer(hpc_sparse_matrix, p, Ap);
    //-- RHT Not replicated --TOCK(t3);

    //-- RHT Not replicated --TICK();
    waxpby_consumer(nrow, 1.0, b, -1.0, Ap, r);
    //-- RHT Not replicated --TOCK(t2);

    //-- RHT Not replicated --TICK();
    ddot_consumer(nrow, r, r, &rtrans, t4);
    //-- RHT Not replicated --TOCK(t1);

    normr = sqrt(rtrans);

    /*-- RHT -- */ RHT_Consume_Check(rank);
    /*-- RHT Volatile -- */ RHT_Consume_Volatile(normr);
    /*-- RHT Not Replicated -- */// if (rank == 0) cout << "Initial Residual = " << normr << endl;

    //-- RHT Not replicated --TICK();
    waxpby_consumer(nrow, 1.0, r, 0.0, r, p);
    //-- RHT Not replicated --TOCK(t2);

    double beta = 0;
    int k = 1;

    goto inFor;

    for (; k < max_iter && normr > tolerance; k++) {
        oldrtrans = rtrans;
        /*-- RHT -- */ RHT_Consume_Check(print_freq);

        //-- RHT Not replicated --TICK();
        ddot_consumer(nrow, r, r, &rtrans, t4);
        //-- RHT Not replicated --TOCK(t1);// 2*nrow ops

        beta = rtrans / oldrtrans;
        /*-- RHT -- */ RHT_Consume_Check(beta);

        //-- RHT Not replicated --TICK();
        waxpby_consumer(nrow, 1.0, r, beta, p, p);
        //-- RHT Not replicated --TOCK(t2);// 2*nrow ops

      inFor:

        normr = sqrt(rtrans);
#if JUST_VOLATILES == 1
        /*-- RHT -- */ RHT_Consume_Volatile(k);
        /*-- RHT -- */ RHT_Consume_Volatile(rank);
        /*-- RHT -- */ RHT_Consume_Volatile(print_freq);
        /*-- RHT -- */ RHT_Consume_Volatile(max_iter);
#else
        /*-- RHT -- */ RHT_Consume_Check(k);
        /*-- RHT -- */ RHT_Consume_Check(rank);
        /*-- RHT -- */ RHT_Consume_Check(print_freq);
        /*-- RHT -- */ RHT_Consume_Check(max_iter);
#endif
        /*-- RHT Volatile -- */ RHT_Consume_Volatile(normr);
        /*-- RHT Not replicated -- */// if (rank == 0 && (k % print_freq == 0 || k + 1 == max_iter))
        // cout << "Iteration = " << k << "   Residual = " << normr << endl;

#ifdef USING_MPI
        //-- RHT Not replicated --TICK();
        exchange_externals_consumer(hpc_sparse_matrix, p);
        //-- RHT Not replicated --TOCK(t5);
#endif
        //-- RHT Not replicated --TICK();
        HPC_sparsemv_consumer(hpc_sparse_matrix, p, Ap);
        //-- RHT Not replicated --TOCK(t3); // 2*nnz ops

        double alpha = 0.0;
        /*-- RHT -- */ RHT_Consume_Check(alpha);

        //-- RHT Not replicated --TICK();
        ddot_consumer(nrow, p, Ap, &alpha, t4);
        //-- RHT Not replicated --TOCK(t1); // 2*nrow ops

        alpha = rtrans / alpha;
        /*-- RHT -- */ RHT_Consume_Check(alpha);

        //-- RHT Not replicated --TICK();
        waxpby_consumer(nrow, 1.0, x, alpha, p, x);// 2*nrow ops
        waxpby_consumer(nrow, 1.0, r, -alpha, Ap, r);
        //-- RHT Not replicated --TOCK(t2);// 2*nrow ops

        niters = k;
        RHT_Consume_Check(niters);
    }

    /// dperez, Times not needed to be replicated
    // Store times
    times[1] = t1; // ddot time
    times[2] = t2; // waxpby time
    times[3] = t3; // sparsemv time
    times[4] = t4; // AllReduce time
#ifdef USING_MPI
    times[5] = t5; // exchange boundary time
#endif
    delete[] p;
    delete[] Ap;
    delete[] r;
    //-- RHT -- Not replicated times[0] = mytimer() - t_begin;  // Total time. All done...
    return (0);
}
