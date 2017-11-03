//
// Created by diego on 30/10/17.
//

#ifndef HPCCG_1_0_SYNCQUEUE_H
#define HPCCG_1_0_SYNCQUEUE_H

#include <cstdlib>
#include <cstdio>
#include <pthread.h>

#define SYNC_QUEUE_SIZE 1024
#define ALREADY_CONSUMED -2

//void SetThreadAffinity(int threadId);

typedef struct{
    volatile int deqPtr;
    double padding0[15];
    volatile int enqPtr;
    double padding1[15];
    volatile double* content;
    double padding2[15];
    volatile int checkState;
    double padding3[15];
    volatile double volatileValue;
}SyncQueue;

SyncQueue SyncQueue_Init();
void SyncQueue_Produce_Simple(double value);
void SyncQueue_Produce_Volatile(double value);

double SyncQueue_Consume();
void SyncQueue_Consume_Check(double currentValue);
void SyncQueue_Consume_Volatile(double currentValue);


void Replication_Init(int numThreads);
void Replication_Finish();

//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// GLOBAL VARS ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

extern SyncQueue globaQueue;

extern volatile long producerCount;
extern volatile long consumerCount;

extern pthread_t ** consumerThreads;
extern int NUM_CONSUMER_THREADS;

#endif //HPCCG_1_0_SYNCQUEUE_H
