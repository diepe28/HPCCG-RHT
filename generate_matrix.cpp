
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

// Routine to read a sparse matrix, right hand side, initial guess, 
// and exact solution (as computed by a direct solver).

/////////////////////////////////////////////////////////////////////////

// nrow - number of rows of matrix (on this processor)

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include "generate_matrix.h"

void generate_matrix(int nx, int ny, int nz, HPC_Sparse_Matrix **A, double **x, double **b, double **xexact) {
#ifdef DEBUG
    int debug = 1;
#else
    int debug = 0;
#endif

#ifdef USING_MPI
    int size, rank; // Number of MPI processes, My process ID
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    int size = 1; // Serial case (not using MPI)
int rank = 0;
#endif

    *A = new HPC_Sparse_Matrix; // Allocate matrix struct and fill it
    (*A)->title = 0;


    // Set this bool to true if you want a 7-pt stencil instead of a 27 pt stencil
    bool use_7pt_stencil = false;

    int local_nrow = nx * ny * nz; // This is the size of our subblock
    assert(local_nrow > 0); // Must have something to work with
    int local_nnz = 27 * local_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)

    int total_nrow = local_nrow * size; // Total number of grid points in mesh
    long long total_nnz =
            27 * (long long) total_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)

    int start_row = local_nrow * rank; // Each processor gets a section of a chimney stack domain
    int stop_row = start_row + local_nrow - 1;


    // Allocate arrays that are of length local_nrow
    (*A)->nnz_in_row = new int[local_nrow];
    (*A)->ptr_to_vals_in_row = new double *[local_nrow];
    (*A)->ptr_to_inds_in_row = new int *[local_nrow];
    (*A)->ptr_to_diags = new double *[local_nrow];

    *x = new double[local_nrow];
    *b = new double[local_nrow];
    *xexact = new double[local_nrow];


    // Allocate arrays that are of length local_nnz
    (*A)->list_of_vals = new double[local_nnz];
    (*A)->list_of_inds = new int[local_nnz];

    double *curvalptr = (*A)->list_of_vals;
    int *curindptr = (*A)->list_of_inds;

    long long nnzglobal = 0;
    for (int iz = 0; iz < nz; iz++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                int curlocalrow = iz * nx * ny + iy * nx + ix;
                int currow = start_row + iz * nx * ny + iy * nx + ix;
                int nnzrow = 0;
                (*A)->ptr_to_vals_in_row[curlocalrow] = curvalptr;
                (*A)->ptr_to_inds_in_row[curlocalrow] = curindptr;
                for (int sz = -1; sz <= 1; sz++) {
                    for (int sy = -1; sy <= 1; sy++) {
                        for (int sx = -1; sx <= 1; sx++) {
                            int curcol = currow + sz * nx * ny + sy * nx + sx;
//            Since we have a stack of nx by ny by nz domains , stacking in the z direction, we check to see
//            if sx and sy are reaching outside of the domain, while the check for the curcol being valid
//            is sufficient to check the z values
                            if ((ix + sx >= 0) && (ix + sx < nx) && (iy + sy >= 0) && (iy + sy < ny) &&
                                (curcol >= 0 && curcol < total_nrow)) {
                                if (!use_7pt_stencil || (sz * sz + sy * sy + sx * sx <=
                                                         1)) { // This logic will skip over point that are not part of a 7-pt stencil
                                    if (curcol == currow) {
                                        (*A)->ptr_to_diags[curlocalrow] = curvalptr;
                                        *curvalptr++ = 27.0;
                                    } else {
                                        *curvalptr++ = -1.0;
                                    }
                                    *curindptr++ = curcol;
                                    nnzrow++;
                                }
                            }
                        } // end sx loop
                    } // end sy loop
                } // end sz loop
                (*A)->nnz_in_row[curlocalrow] = nnzrow;
                nnzglobal += nnzrow;
                (*x)[curlocalrow] = 0.0;
                (*b)[curlocalrow] = 27.0 - ((double) (nnzrow - 1));
                (*xexact)[curlocalrow] = 1.0;
            } // end ix loop
        } // end iy loop
    } // end iz loop
    if (debug) cout << "Process " << rank << " of " << size << " has " << local_nrow;

    if (debug)
        cout << " rows. Global rows " << start_row
             << " through " << stop_row << endl;

    if (debug)
        cout << "Process " << rank << " of " << size
             << " has " << local_nnz << " nonzeros." << endl;

    (*A)->start_row = start_row;
    (*A)->stop_row = stop_row;
    (*A)->total_nrow = total_nrow;
    (*A)->total_nnz = total_nnz;
    (*A)->local_nrow = local_nrow;
    (*A)->local_ncol = local_nrow;
    (*A)->local_nnz = local_nnz;

    return;
}

void generate_matrix_producer(int nx, int ny, int nz, HPC_Sparse_Matrix **A, double **x, double **b, double **xexact) {
#ifdef DEBUG
    int debug = 1;
#else
    int debug = 0;
#endif

#ifdef USING_MPI
    int size, rank; // Number of MPI processes, My process ID
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    SyncQueue_Simple_Produce(size);
    SyncQueue_Simple_Produce(rank);
#else
    int size = 1; // Serial case (not using MPI)
    int rank = 0;
#endif

    *A = new HPC_Sparse_Matrix; // Allocate matrix struct and fill it
    (*A)->title = 0;

    ///*-- RHT -- */ SyncQueue_Simple_Produce((int) (*A)->title);

    // Set this bool to true if you want a 7-pt stencil instead of a 27 pt stencil
    bool use_7pt_stencil = false;
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) use_7pt_stencil);

    int local_nrow = nx * ny * nz; // This is the size of our subblock
    assert(local_nrow > 0); // Must have something to work with
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) local_nrow);

    int local_nnz = 27 * local_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) local_nnz);

    int total_nrow = local_nrow * size; // Total number of grid points in mesh
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) total_nrow);

    long long total_nnz = 27 * (long long) total_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) total_nnz);

    int start_row = local_nrow * rank; // Each processor gets a section of a chimney stack domain
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) start_row);

    int stop_row = start_row + local_nrow - 1;
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) stop_row);

    // TODO, what about these?
    // Allocate arrays that are of length local_nrow
    (*A)->nnz_in_row = new int[local_nrow];
    (*A)->ptr_to_vals_in_row = new double *[local_nrow];
    (*A)->ptr_to_inds_in_row = new int *[local_nrow];
    (*A)->ptr_to_diags = new double *[local_nrow];

    *x = new double[local_nrow];
    *b = new double[local_nrow];
    *xexact = new double[local_nrow];

    // Allocate arrays that are of length local_nnz
    (*A)->list_of_vals = new double[local_nnz];
    (*A)->list_of_inds = new int[local_nnz];

    double *curvalptr = (*A)->list_of_vals;
    int *curindptr = (*A)->list_of_inds;

    long long nnzglobal = 0;
    for (int iz = 0; iz < nz; iz++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                int curlocalrow = iz * nx * ny + iy * nx + ix;
                /*-- RHT -- */ SyncQueue_Simple_Produce((double) curlocalrow);

                int currow = start_row + iz * nx * ny + iy * nx + ix;
                /*-- RHT -- */ SyncQueue_Simple_Produce((double) currow);

                int nnzrow = 0;
                /*-- RHT -- */ SyncQueue_Simple_Produce((double) nnzrow);

                // TODO, should it be checked?, they are arrays
                (*A)->ptr_to_vals_in_row[curlocalrow] = curvalptr;
                (*A)->ptr_to_inds_in_row[curlocalrow] = curindptr;

                for (int sz = -1; sz <= 1; sz++) {
                    for (int sy = -1; sy <= 1; sy++) {
                        for (int sx = -1; sx <= 1; sx++) {
                            int curcol = currow + sz * nx * ny + sy * nx + sx;
                            /*-- RHT -- */ SyncQueue_Simple_Produce((double) curcol);

//            Since we have a stack of nx by ny by nz domains , stacking in the z direction, we check to see
//            if sx and sy are reaching outside of the domain, while the check for the curcol being valid
//            is sufficient to check the z values
                            if ((ix + sx >= 0) && (ix + sx < nx) && (iy + sy >= 0) && (iy + sy < ny) &&
                                (curcol >= 0 && curcol < total_nrow)) {
                                if (!use_7pt_stencil || (sz * sz + sy * sy + sx * sx <=
                                                         1)) { // This logic will skip over point that are not part of a 7-pt stencil
                                    if (curcol == currow) {
                                        (*A)->ptr_to_diags[curlocalrow] = curvalptr;
                                        *curvalptr++ = 27.0;
                                    } else {
                                        *curvalptr++ = -1.0;
                                    }
                                    // TODO, what about the ++, also should we check for the result of the if?
                                    /*-- RHT -- */ SyncQueue_Simple_Produce((double) *curvalptr-1);

                                    *curindptr++ = curcol;
                                    /*-- RHT -- */ SyncQueue_Simple_Produce((double) *curindptr-1);

                                    nnzrow++;
                                    /*-- RHT -- */ SyncQueue_Simple_Produce((double) *curindptr-1);
                                }
                            }
                        } // end sx loop
                    } // end sy loop
                } // end sz loop
                (*A)->nnz_in_row[curlocalrow] = nnzrow;
                /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*A)->nnz_in_row[curlocalrow]);

                nnzglobal += nnzrow;
                /*-- RHT -- */ SyncQueue_Simple_Produce((double) nnzglobal);

                (*x)[curlocalrow] = 0.0;
                /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*x)[curlocalrow]);

                (*b)[curlocalrow] = 27.0 - ((double) (nnzrow - 1));
                /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*b)[curlocalrow]);

                (*xexact)[curlocalrow] = 1.0;
                /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*xexact)[curlocalrow]);
            } // end ix loop
        } // end iy loop
    } // end iz loop
    if (debug) cout << "Process " << rank << " of " << size << " has " << local_nrow;

    if (debug)
        cout << " rows. Global rows " << start_row
             << " through " << stop_row << endl;

    if (debug)
        cout << "Process " << rank << " of " << size
             << " has " << local_nnz << " nonzeros." << endl;

    (*A)->start_row = start_row;
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*A)->start_row);

    (*A)->stop_row = stop_row;
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*A)->stop_row);

    (*A)->total_nrow = total_nrow;
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*A)->total_nrow);

    (*A)->total_nnz = total_nnz;
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*A)->total_nnz);

    (*A)->local_nrow = local_nrow;
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*A)->local_nrow);

    (*A)->local_ncol = local_nrow;
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*A)->local_ncol);

    (*A)->local_nnz = local_nnz;
    /*-- RHT -- */ SyncQueue_Simple_Produce((double) (*A)->local_nnz);

    return;
}

void generate_matrix_consumer(int nx, int ny, int nz, HPC_Sparse_Matrix **A, double **x, double **b, double **xexact, int size, int rank){

#ifdef DEBUG
    int debug = 1;
#else
    int debug = 0;
#endif

    *A = new HPC_Sparse_Matrix; // Allocate matrix struct and fill it
    (*A)->title = 0;

    ///*-- RHT -- */ SyncQueue_Simple_Produce((double) (*A)->title);

    // Set this bool to true if you want a 7-pt stencil instead of a 27 pt stencil
    bool use_7pt_stencil = false;
    /*-- RHT -- */ SyncQueue_Consume_Check((double) use_7pt_stencil);

    int local_nrow = nx * ny * nz; // This is the size of our subblock
    assert(local_nrow > 0); // Must have something to work with
    /*-- RHT -- */ SyncQueue_Consume_Check((double) local_nrow);

    int local_nnz = 27 * local_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)
    /*-- RHT -- */ SyncQueue_Consume_Check((double) local_nnz);

    int total_nrow = local_nrow * size; // Total number of grid points in mesh
    /*-- RHT -- */ SyncQueue_Consume_Check((double) total_nrow);

    long long total_nnz = 27 * (long long) total_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)
    /*-- RHT -- */ SyncQueue_Consume_Check((double) total_nnz);

    int start_row = local_nrow * rank; // Each processor gets a section of a chimney stack domain
    /*-- RHT -- */ SyncQueue_Consume_Check((double) start_row);

    int stop_row = start_row + local_nrow - 1;
    /*-- RHT -- */ SyncQueue_Consume_Check((double) stop_row);

    // TODO, what about these?
    // Allocate arrays that are of length local_nrow
    (*A)->nnz_in_row = new int[local_nrow];
    (*A)->ptr_to_vals_in_row = new double *[local_nrow];
    (*A)->ptr_to_inds_in_row = new int *[local_nrow];
    (*A)->ptr_to_diags = new double *[local_nrow];

    *x = new double[local_nrow];
    *b = new double[local_nrow];
    *xexact = new double[local_nrow];

    // Allocate arrays that are of length local_nnz
    (*A)->list_of_vals = new double[local_nnz];
    (*A)->list_of_inds = new int[local_nnz];

    double *curvalptr = (*A)->list_of_vals;
    int *curindptr = (*A)->list_of_inds;

    long long nnzglobal = 0;
    for (int iz = 0; iz < nz; iz++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                int curlocalrow = iz * nx * ny + iy * nx + ix;
                /*-- RHT -- */ SyncQueue_Consume_Check((double) curlocalrow);

                int currow = start_row + iz * nx * ny + iy * nx + ix;
                /*-- RHT -- */ SyncQueue_Consume_Check((double) currow);

                int nnzrow = 0;
                /*-- RHT -- */ SyncQueue_Consume_Check((double) nnzrow);

                // TODO, should it be checked?, they are arrays
                (*A)->ptr_to_vals_in_row[curlocalrow] = curvalptr;
                (*A)->ptr_to_inds_in_row[curlocalrow] = curindptr;

                for (int sz = -1; sz <= 1; sz++) {
                    for (int sy = -1; sy <= 1; sy++) {
                        for (int sx = -1; sx <= 1; sx++) {
                            int curcol = currow + sz * nx * ny + sy * nx + sx;
                            /*-- RHT -- */ SyncQueue_Consume_Check((double) curcol);

//            Since we have a stack of nx by ny by nz domains , stacking in the z direction, we check to see
//            if sx and sy are reaching outside of the domain, while the check for the curcol being valid
//            is sufficient to check the z values
                            if ((ix + sx >= 0) && (ix + sx < nx) && (iy + sy >= 0) && (iy + sy < ny) &&
                                (curcol >= 0 && curcol < total_nrow)) {
                                if (!use_7pt_stencil || (sz * sz + sy * sy + sx * sx <=
                                                         1)) { // This logic will skip over point that are not part of a 7-pt stencil
                                    if (curcol == currow) {
                                        (*A)->ptr_to_diags[curlocalrow] = curvalptr;
                                        *curvalptr++ = 27.0;
                                    } else {
                                        *curvalptr++ = -1.0;
                                    }
                                    // TODO, what about the ++, also should we check for the result of the if?
                                    /*-- RHT -- */ SyncQueue_Consume_Check((double) *curvalptr-1);

                                    *curindptr++ = curcol;
                                    /*-- RHT -- */ SyncQueue_Consume_Check((double) *curindptr-1);

                                    nnzrow++;
                                    /*-- RHT -- */ SyncQueue_Consume_Check((double) *curindptr-1);
                                }
                            }
                        } // end sx loop
                    } // end sy loop
                } // end sz loop
                (*A)->nnz_in_row[curlocalrow] = nnzrow;
                /*-- RHT -- */ SyncQueue_Consume_Check((double) (*A)->nnz_in_row[curlocalrow]);

                nnzglobal += nnzrow;
                /*-- RHT -- */ SyncQueue_Consume_Check((double) nnzglobal);

                (*x)[curlocalrow] = 0.0;
                /*-- RHT -- */ SyncQueue_Consume_Check((double) (*x)[curlocalrow]);

                (*b)[curlocalrow] = 27.0 - ((double) (nnzrow - 1));
                /*-- RHT -- */ SyncQueue_Consume_Check((double) (*b)[curlocalrow]);

                (*xexact)[curlocalrow] = 1.0;
                /*-- RHT -- */ SyncQueue_Consume_Check((double) (*xexact)[curlocalrow]);
            } // end ix loop
        } // end iy loop
    } // end iz loop
    if (debug) cout << "Process " << rank << " of " << size << " has " << local_nrow;

    if (debug)
        cout << " rows. Global rows " << start_row
             << " through " << stop_row << endl;

    if (debug)
        cout << "Process " << rank << " of " << size
             << " has " << local_nnz << " nonzeros." << endl;

    (*A)->start_row = start_row;
    /*-- RHT -- */ SyncQueue_Consume_Check((double) (*A)->start_row);

    (*A)->stop_row = stop_row;
    /*-- RHT -- */ SyncQueue_Consume_Check((double) (*A)->stop_row);

    (*A)->total_nrow = total_nrow;
    /*-- RHT -- */ SyncQueue_Consume_Check((double) (*A)->total_nrow);

    (*A)->total_nnz = total_nnz;
    /*-- RHT -- */ SyncQueue_Consume_Check((double) (*A)->total_nnz);

    (*A)->local_nrow = local_nrow;
    /*-- RHT -- */ SyncQueue_Consume_Check((double) (*A)->local_nrow);

    (*A)->local_ncol = local_nrow;
    /*-- RHT -- */ SyncQueue_Consume_Check((double) (*A)->local_ncol);

    (*A)->local_nnz = local_nnz;
    /*-- RHT -- */ SyncQueue_Consume_Check((double) (*A)->local_nnz);

    return;
}

