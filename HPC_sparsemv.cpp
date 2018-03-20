
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

// Routine to compute matrix vector product y = Ax where:
// First call exchange_externals to get off-processor values of x

// A - known matrix 
// x - known vector
// y - On exit contains Ax.

/////////////////////////////////////////////////////////////////////////

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
#include "HPC_sparsemv.h"
#include "RHT.h"

int HPC_sparsemv( HPC_Sparse_Matrix *hpc_sparse_matrix,
		 const double * const x, double * const y) {

    const int nrow = (const int) hpc_sparse_matrix->local_nrow;

#ifdef USING_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nrow; i++) {
        double sum = 0.0;
        const double *const cur_vals =
                (const double *const) hpc_sparse_matrix->ptr_to_vals_in_row[i];

        const int *const cur_inds =
                (const int *const) hpc_sparse_matrix->ptr_to_inds_in_row[i];

        const int cur_nnz = (const int) hpc_sparse_matrix->nnz_in_row[i];

        for (int j = 0; j < cur_nnz; j++)
            sum += cur_vals[j] * x[cur_inds[j]];
        y[i] = sum;
    }
    return (0);
}

int HPC_sparsemv_producer_no_sync( HPC_Sparse_Matrix *hpc_sparse_matrix,
                  const double * const x, double * const y) {

    const int nrow = (const int) hpc_sparse_matrix->local_nrow;
    /*-- RHT -- */ RHT_Produce_Secure(nrow);
#ifdef USING_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nrow; i++) {
        double sum = 0.0;
        /*-- RHT -- */ RHT_Produce_Secure(sum);

        const double *const cur_vals = (const double *const) hpc_sparse_matrix->ptr_to_vals_in_row[i];
        /*-- RHT -- */ RHT_Produce_Secure(*cur_vals);

        const int *const cur_inds = (const int *const) hpc_sparse_matrix->ptr_to_inds_in_row[i];
        /*-- RHT -- */ RHT_Produce_Secure(*cur_inds);

        const int cur_nnz = (const int) hpc_sparse_matrix->nnz_in_row[i];
        /*-- RHT -- */ RHT_Produce_Secure(cur_nnz);

        int j;
        replicate_loop_producer(cur_nnz, j, sum, sum += cur_vals[j] * x[cur_inds[j]])

        y[i] = sum;
        /*-- RHT -- */ RHT_Produce_Secure(y[i]);
    }
    return (0);
}

int HPC_sparsemv_producer( HPC_Sparse_Matrix *hpc_sparse_matrix,
                           const double * const x, double * const y) {

    const int nrow = (const int) hpc_sparse_matrix->local_nrow;
    /*-- RHT -- */ RHT_Produce(nrow);
#ifdef USING_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nrow; i++) {
        double sum = 0.0;
        /*-- RHT -- */ RHT_Produce(sum);

        const double *const cur_vals = (const double *const) hpc_sparse_matrix->ptr_to_vals_in_row[i];
        /*-- RHT -- */ RHT_Produce(*cur_vals);

        const int *const cur_inds = (const int *const) hpc_sparse_matrix->ptr_to_inds_in_row[i];
        /*-- RHT -- */ RHT_Produce(*cur_inds);

        const int cur_nnz = (const int) hpc_sparse_matrix->nnz_in_row[i];
        /*-- RHT -- */ RHT_Produce(cur_nnz);

        for (int j = 0; j < cur_nnz; j++) {
            sum += cur_vals[j] * x[cur_inds[j]];
            /*-- RHT -- */ RHT_Produce(sum);
        }

        y[i] = sum;
        /*-- RHT -- */ RHT_Produce(y[i]);
    }
    return (0);
}

int HPC_sparsemv_consumer( HPC_Sparse_Matrix *hpc_sparse_matrix,
                           const double * const x, double * const y) {

    //printf("nrow in consumer %d: \n", hpc_sparse_matrix->local_nrow);
    const int nrow = (const int) hpc_sparse_matrix->local_nrow;
    /*-- RHT -- */ RHT_Consume_Check(nrow);
#ifdef USING_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nrow; i++) {
        double sum = 0.0;
        /*-- RHT -- */ RHT_Consume_Check(sum);

        const double *const cur_vals = (const double *const) hpc_sparse_matrix->ptr_to_vals_in_row[i];
        /*-- RHT -- */ RHT_Consume_Check(*cur_vals);

        const int *const cur_inds = (const int *const) hpc_sparse_matrix->ptr_to_inds_in_row[i];
        /*-- RHT -- */ RHT_Consume_Check(*cur_inds);

        const int cur_nnz = (const int) hpc_sparse_matrix->nnz_in_row[i];
        /*-- RHT -- */ RHT_Consume_Check(cur_nnz);

        int j;
#if VAR_GROUPING == 1
        replicate_loop_consumer(cur_nnz, j, sum, sum += cur_vals[j] * x[cur_inds[j]])
#else
        for (j = 0; j < cur_nnz; j++) {
            sum += cur_vals[j] * x[cur_inds[j]];
            /*-- RHT -- */ RHT_Consume_Check(sum);
        }
#endif
        y[i] = sum;
        /*-- RHT -- */ RHT_Consume_Check(y[i]);
    }
    return (0);
}


