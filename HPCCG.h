
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

#ifndef HPCCG_H
#define HPCCG_H
#include "HPC_sparsemv.h"
#include "ddot.h"
#include "waxpby.h"
#include "HPC_Sparse_Matrix.h"

#ifdef USING_MPI
#include "exchange_externals.h"
#include <mpi.h> // If this routine is compiled with -DUSING_MPI
                 // then include mpi.h
#endif
#include "FlipIt/corrupt/corrupt.h"
int HPCCG(HPC_Sparse_Matrix * hpc_sparse_matrix,
	  const double * const b, double * const x,
	  const int max_iter, const double tolerance, int & niters, double & normr, double * times);

// this function will compute the Conjugate Gradient...
// A <=> Matrix
// b <=> constant
// xnot <=> initial guess
// max_iter <=> how many times we iterate
// tolerance <=> specifies how "good"of a value we would like
// x <=> used for return value

// A is known
// x is unknown vector
// b is known vector
// xnot = 0
// niters is the number of iterations

int HPCCG_producer(HPC_Sparse_Matrix *hpc_sparse_matrix,
                            const double *const b, double *const x,
                            const int max_iter, const double tolerance, int &niters, double &normr, double *times);

int HPCCG_consumer(HPC_Sparse_Matrix * hpc_sparse_matrix,
		  const double * const b, double * const x,
		  const int max_iter, const double tolerance, int & niters, double & normr, double * times);
#endif
