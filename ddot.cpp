
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

// Routine to compute the dot product of two vectors where:

// n - number of vector elements (on this processor)

// x, y - input vectors

// residual - pointer to scalar value, on exit will contain result.

/////////////////////////////////////////////////////////////////////////

#include "ddot.h"

int ddot (const int n, const double * const x, const double * const y,
	  double * const result, double & time_allreduce) {
    double local_result = 0.0;
    if (y == x)
#ifdef USING_OMP
#pragma omp parallel for reduction (+:local_result)
#endif
        for (int i = 0; i < n; i++) local_result += x[i] * x[i];
    else
#ifdef USING_OMP
#pragma omp parallel for reduction (+:local_result)
#endif
        for (int i = 0; i < n; i++) local_result += x[i] * y[i];

#ifdef USING_MPI
    // Use MPI's reduce function to collect all partial sums
    double t0 = mytimer();
    double global_result = 0.0;
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    *result = global_result;
    time_allreduce += mytimer() - t0;
#else
    *result = local_result;
#endif

    return (0);
}

int ddot_producer (const int n, const double * const x, const double * const y,
          double * const result, double & time_allreduce) {
    double local_result = 0.0;
    if (y == x)
#ifdef USING_OMP
#pragma omp parallel for reduction (+:local_result)
#endif
        for (int i = 0; i < n; i++) {
            local_result += x[i] * x[i];
            /*-- RHT -- */ SyncQueue_Produce_Simple(local_result);
        }
    else
#ifdef USING_OMP
#pragma omp parallel for reduction (+:local_result)
#endif
        for (int i = 0; i < n; i++) {
            local_result += x[i] * y[i];
            /*-- RHT -- */ SyncQueue_Produce_Simple(local_result);
        }

#ifdef USING_MPI
    // Use MPI's reduce function to collect all partial sums
    double t0 = mytimer();
    double global_result = 0.0;

    /*-- RHT Volatile -- */ SyncQueue_Produce_Volatile(global_result);
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    /*-- RHT -- */ SyncQueue_Produce_Simple(global_result);

    *result = global_result;

    time_allreduce += mytimer() - t0;
#else
    *result = local_result;
#endif

    /*-- RHT -- */ SyncQueue_Produce_Simple(*result);
    return (0);
}

int ddot_consumer (const int n, const double * const x, const double * const y,
                   double * const result, double & time_allreduce) {
    double local_result = 0.0;
    if (y == x)
#ifdef USING_OMP
#pragma omp parallel for reduction (+:local_result)
#endif
        for (int i = 0; i < n; i++) {
            local_result += x[i] * x[i];
            /*-- RHT -- */ SyncQueue_Consume_Check(local_result);
        }
    else
#ifdef USING_OMP
#pragma omp parallel for reduction (+:local_result)
#endif
        for (int i = 0; i < n; i++) {
            local_result += x[i] * y[i];
            /*-- RHT -- */ SyncQueue_Consume_Check(local_result);
        }

#ifdef USING_MPI
    // Use MPI's reduce function to collect all partial sums
    double t0 = mytimer();
    double global_result = 0.0;

    /*-- RHT Volatile -- */ SyncQueue_Consume_Volatile(global_result);
    /*-- RHT Not replicated -- */// MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    global_result = SyncQueue_Consume();

    *result = global_result;
    time_allreduce += mytimer() - t0;
#else
    *result = local_result;
#endif

    /*-- RHT -- */ SyncQueue_Consume_Check(*result);
    return (0);
}
