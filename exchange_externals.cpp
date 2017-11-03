
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

#ifdef USING_MPI  // Compile this routine only if running in parallel
#include <iostream>
using std::cerr;
using std::endl;
#include <cstdlib>
#include <cstdio>
#include "exchange_externals.h"

#undef DEBUG

void exchange_externals(HPC_Sparse_Matrix * A, const double *x) {
    int i, j, k;
    int num_external = 0;

    // Extract Matrix pieces

    int local_nrow = A->local_nrow;
    int num_neighbors = A->num_send_neighbors;
    int *recv_length = A->recv_length;
    int *send_length = A->send_length;
    int *neighbors = A->neighbors;
    double *send_buffer = A->send_buffer;
    int total_to_be_sent = A->total_to_be_sent;
    int *elements_to_send = A->elements_to_send;

    int size, rank; // Number of MPI processes, My process ID
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //
    //  first post receives, these are immediate receives
    //  Do not wait for result to come, will do that at the
    //  wait call below.
    //

    int MPI_MY_TAG = 99;

    MPI_Request *request = new MPI_Request[num_neighbors];

    //
    // Externals are at end of locals
    //
    double *x_external = (double *) x + local_nrow;

    // Post receives first
    for (i = 0; i < num_neighbors; i++) {
        int n_recv = recv_length[i];
        MPI_Irecv(x_external, n_recv, MPI_DOUBLE, neighbors[i], MPI_MY_TAG,
                  MPI_COMM_WORLD, request + i);
        x_external += n_recv;
    }


    //
    // Fill up send buffer
    //

    for (i = 0; i < total_to_be_sent; i++) send_buffer[i] = x[elements_to_send[i]];

    //
    // Send to each neighbor
    //

    for (i = 0; i < num_neighbors; i++) {
        int n_send = send_length[i];
        MPI_Send(send_buffer, n_send, MPI_DOUBLE, neighbors[i], MPI_MY_TAG,
                 MPI_COMM_WORLD);
        send_buffer += n_send;
    }

    //
    // Complete the reads issued above
    //

    MPI_Status status;
    for (i = 0; i < num_neighbors; i++) {
        if (MPI_Wait(request + i, &status)) {
            cerr << "MPI_Wait error\n" << endl;
            exit(-1);
        }
    }

    delete[] request;

    return;
}


void exchange_externals_producer(HPC_Sparse_Matrix * A, const double *x) {
    int i, j, k;
    int num_external = 0;

    // Extract Matrix pieces

    int local_nrow = A->local_nrow;
    /*-- RHT -- */ SyncQueue_Produce_Simple(local_nrow);

    int num_neighbors = A->num_send_neighbors;
    /*-- RHT -- */ SyncQueue_Produce_Simple(num_neighbors);

    /// TODO what to do with arrays
    int *recv_length = A->recv_length;
    int *send_length = A->send_length;
    int *neighbors = A->neighbors;
    double *send_buffer = A->send_buffer;

    int total_to_be_sent = A->total_to_be_sent;
    /*-- RHT -- */ SyncQueue_Produce_Simple(total_to_be_sent);

    int *elements_to_send = A->elements_to_send;

    int size, rank; // Number of MPI processes, My process ID
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /*-- RHT -- */ SyncQueue_Produce_Simple(size);
    /*-- RHT -- */ SyncQueue_Produce_Simple(rank);

    //
    //  first post receives, these are immediate receives
    //  Do not wait for result to come, will do that at the
    //  wait call below.
    //

    int MPI_MY_TAG = 99;
    /*-- RHT -- */ SyncQueue_Produce_Simple(MPI_MY_TAG);

    MPI_Request *request = new MPI_Request[num_neighbors];

    //
    // Externals are at end of locals
    //
    double *x_external = (double *) x + local_nrow;
    /*-- RHT -- */ SyncQueue_Produce_Simple(*x_external);

    // Post receives first
    for (i = 0; i < num_neighbors; i++) {
        int n_recv = recv_length[i];
        /*-- RHT -- */ SyncQueue_Produce_Simple(n_recv);
        /*-- RHT Volatile -- */ SyncQueue_Produce_Volatile(neighbors[i]);
        /// TODO what to with last parameter? is a user def type
        MPI_Irecv(x_external, n_recv, MPI_DOUBLE, neighbors[i], MPI_MY_TAG,
                  MPI_COMM_WORLD, request + i);


        x_external += n_recv;
        /*-- RHT -- */ SyncQueue_Produce_Simple(*x_external);
    }


    //
    // Fill up send buffer
    //

    for (i = 0; i < total_to_be_sent; i++) send_buffer[i] = x[elements_to_send[i]];

    //
    // Send to each neighbor
    //

    for (i = 0; i < num_neighbors; i++) {
        int n_send = send_length[i];
        /*-- RHT -- */ SyncQueue_Produce_Simple(n_send);
        /*-- RHT Volatile -- */ SyncQueue_Produce_Volatile(neighbors[i]);
        MPI_Send(send_buffer, n_send, MPI_DOUBLE, neighbors[i], MPI_MY_TAG,
                 MPI_COMM_WORLD);
        send_buffer += n_send;
        /*-- RHT -- */ SyncQueue_Produce_Simple(n_send);
    }

    //
    // Complete the reads issued above
    //

    MPI_Status status;
    for (i = 0; i < num_neighbors; i++) {
        /// TODO, what to with this parameters, should we check every single value
        if (MPI_Wait(request + i, &status)) {
            cerr << "MPI_Wait error\n" << endl;
            exit(-1);
        }
    }

    delete[] request;

    return;
}

void exchange_externals_consumer(HPC_Sparse_Matrix * A, const double *x) {
    int i, j, k;
    int num_external = 0;

    // Extract Matrix pieces

    int local_nrow = A->local_nrow;
    /*-- RHT -- */ SyncQueue_Consume_Check(local_nrow);

    int num_neighbors = A->num_send_neighbors;
    /*-- RHT -- */ SyncQueue_Consume_Check(num_neighbors);

    /// TODO what to do with arrays
    int *recv_length = A->recv_length;
    int *send_length = A->send_length;
    int *neighbors = A->neighbors;
    double *send_buffer = A->send_buffer;

    int total_to_be_sent = A->total_to_be_sent;
    /*-- RHT -- */ SyncQueue_Consume_Check(total_to_be_sent);

    int *elements_to_send = A->elements_to_send;

    int size, rank; // Number of MPI processes, My process ID
    /*-- RHT -- */ size = SyncQueue_Consume();
    /*-- RHT -- */ rank = SyncQueue_Consume();

    //
    //  first post receives, these are immediate receives
    //  Do not wait for result to come, will do that at the
    //  wait call below.
    //

    int MPI_MY_TAG = 99;
    /*-- RHT -- */ SyncQueue_Consume_Check(MPI_MY_TAG);

    MPI_Request *request = new MPI_Request[num_neighbors];

    //
    // Externals are at end of locals
    //
    double *x_external = (double *) x + local_nrow;
    /*-- RHT -- */ SyncQueue_Consume_Check(*x_external);

    // Post receives first
    for (i = 0; i < num_neighbors; i++) {
        int n_recv = recv_length[i];
        /*-- RHT -- */ SyncQueue_Consume_Check(n_recv);
        /*-- RHT Volatile -- */ SyncQueue_Consume_Volatile(neighbors[i]);
        /// TODO what to with last parameter? is a user def type
        /*-- RHT Volatile Not replicated -- */// MPI_Irecv(x_external, n_recv, MPI_DOUBLE, neighbors[i], MPI_MY_TAG,
        // MPI_COMM_WORLD, request + i);


        x_external += n_recv;
        /*-- RHT -- */ SyncQueue_Consume_Check(*x_external);
    }


    //
    // Fill up send buffer
    //

    for (i = 0; i < total_to_be_sent; i++) send_buffer[i] = x[elements_to_send[i]];

    //
    // Send to each neighbor
    //

    for (i = 0; i < num_neighbors; i++) {
        int n_send = send_length[i];
        /*-- RHT -- */ SyncQueue_Consume_Check(n_send);
        /*-- RHT Volatile -- */ SyncQueue_Consume_Volatile(neighbors[i]);
        /*-- RHT Volatile Not replicated -- */// MPI_Send(send_buffer, n_send, MPI_DOUBLE, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD);
        send_buffer += n_send;
        /*-- RHT -- */ SyncQueue_Consume_Check(n_send);
    }

    //
    // Complete the reads issued above
    //

    /*-- RHT Volatile Not replicated -- */
//    MPI_Status status;
//    for (i = 0; i < num_neighbors; i++) {
//        /// TODO, what to with this parameters, should we check every single value
//
//        if (MPI_Wait(request + i, &status)) {
//            cerr << "MPI_Wait error\n" << endl;
//            exit(-1);
//        }
//    }

    delete[] request;

    return;
}

#endif // USING_MPI
