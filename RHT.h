//
// Created by diego on 30/10/17.
//

#ifndef HPCCG_1_0_SYNCQUEUE_H
#define HPCCG_1_0_SYNCQUEUE_H

#include <cstdlib>
#include <cstdio>
#include <pthread.h>

#define RHT_QUEUE_SIZE 1024
#define ALREADY_CONSUMED -2

void SetThreadAffinity(int threadId);

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
}RHT_Queue;

void RHT_Produce(double value);
void RHT_Produce_Volatile(double value);

double RHT_Consume();
void RHT_Consume_Check(double currentValue);
void RHT_Consume_Volatile(double currentValue);

void RHT_Replication_Init(int numThreads);
void RHT_Replication_Finish();

//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// GLOBAL VARS ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

extern RHT_Queue globaQueue;

extern volatile long producerCount;
extern volatile long consumerCount;

extern pthread_t ** consumerThreads;
extern int NUM_CONSUMER_THREADS;

#endif //HPCCG_1_0_SYNCQUEUE_H
