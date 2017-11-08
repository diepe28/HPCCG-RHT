    //
// Created by diego on 30/10/17.
//

#ifndef HPCCG_1_0_SYNCQUEUE_H
#define HPCCG_1_0_SYNCQUEUE_H

#include <cstdlib>
#include <cstdio>
#include <pthread.h>

// -------- Macros ----------

extern int nextEnq, localDeq, newLimit, diff;
extern double otherValue, thisValue;


#define RHT_QUEUE_SIZE 1024
#define MIN_PTR_DIST 200
#define ALREADY_CONSUMED -2

//((RHT_QUEUE_SIZE-globalQueue.enqPtr + localDeq) % RHT_QUEUE_SIZE)-1; \
/// with the -1 we make sure that the producer never really catches up to the consumer, but it might still happen
/// the other way around
#define replicate_forLoop_no_sync(numIters, iterator, value, operation)         \
    do{                                                                         \
        asm("pause");                                                           \
        localDeq = globalQueue.deqPtr;                                          \
        newLimit = (globalQueue.enqPtr >= localDeq ?                            \
                   (RHT_QUEUE_SIZE - globalQueue.enqPtr) + localDeq :           \
                   localDeq - globalQueue.enqPtr)-1;                            \
    }                                                                           \
    while(newLimit < MIN_PTR_DIST);                                             \
    iterator = 0;                                                               \
    while (newLimit < numIters) {                                               \
        for (; iterator < newLimit; iterator++){                                \
            operation;                                                          \
            globalQueue.content[globalQueue.enqPtr] = value;                    \
            globalQueue.enqPtr = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;     \
        }                                                                       \
        do{                                                                     \
            asm("pause");                                                       \
            localDeq = globalQueue.deqPtr;                                      \
            diff = (globalQueue.enqPtr >= localDeq ?                            \
                   (RHT_QUEUE_SIZE - globalQueue.enqPtr) + localDeq :           \
                   localDeq - globalQueue.enqPtr)-1;                            \
        }                                                                       \
        while(diff < MIN_PTR_DIST);                                             \
        newLimit += diff;                                                       \
    }                                                                           \
    for (; iterator < numIters; iterator++){                                    \
        operation;                                                              \
        globalQueue.content[globalQueue.enqPtr] = value;                        \
        globalQueue.enqPtr = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;         \
    }

#define Macro_AlreadyConsumed_Produce(value)                  \
    nextEnq = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;      \
    while(globalQueue.content[nextEnq] != ALREADY_CONSUMED)   \
        asm("pause");                                         \
    /*producerCount++;*/                                      \
    globalQueue.content[globalQueue.enqPtr] = value;          \
    globalQueue.enqPtr = nextEnq;

#define Macro_Report_Soft_Error(consumerValue, producerValue) \
    printf("\n SOFT ERROR DETECTED, Consumer: %f Producer: %f -- PCount: %ld , CCount: %ld, diff: %ld \n", \
            consumerValue, producerValue, producerCount, consumerCount, producerCount - consumerCount); \
    exit(1);

#define Macro_AlreadyConsumed_Consume_Check(currentValue)                       \
    thisValue = (double) currentValue;                                          \
    otherValue = globalQueue.content[globalQueue.deqPtr];                         \
    if (thisValue != otherValue) {                                              \
        /* des-sync of the queue */                                             \
        if (otherValue == ALREADY_CONSUMED) {                                   \
            /*consumerCount++;*/                                                \
            do asm("pause"); while (globalQueue.content[globalQueue.deqPtr] == ALREADY_CONSUMED); \
            otherValue = globalQueue.content[globalQueue.deqPtr];                 \
                                                                                \
            if (thisValue == otherValue){                                       \
                globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;       \
                globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;   \
            }else{                                                              \
                 Macro_Report_Soft_Error(otherValue, thisValue)                 \
            }                                                                   \
        }else {                                                                 \
            Macro_Report_Soft_Error(otherValue, thisValue)                      \
        }                                                                       \
    }else{                                                                      \
        /*consumerCount++;*/                                                    \
        globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;               \
        globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;           \
    }

#define Macro_AlreadyConsumed_Consume(value)                         \
    value = globalQueue.content[globalQueue.deqPtr];                 \
    if (value == ALREADY_CONSUMED) {                                 \
        do asm("pause"); while (globalQueue.content[globalQueue.deqPtr] == ALREADY_CONSUMED); \
        value = globalQueue.content[globalQueue.deqPtr];             \
    }                                                                \
    /*consumerCount++;*/                                             \
    globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;      \
    globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;

#define RHT_Produce(value) \
    /*Macro_UsingPointers_Produce(value)*/ \
    Macro_AlreadyConsumed_Produce(value)

//#define RHT_Consume_Check(currentValue) \
//   /*Macro_UsingPointers_Consume_Check(value)*/ \
//   Macro_AlreadyConsumed_Consume_Check(currentValue)

#define RHT_Consume(value) \
    /*Macro_UsingPointers_Consume(value)*/ \
    Macro_AlreadyConsumed_Consume(value)

//void RHT_Produce(double value);
void RHT_Consume_Check(double currentValue);
void RHT_Produce_Secure(double value);

//double RHT_Consume();

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

extern RHT_Queue globalQueue;

void RHT_Produce_Volatile(double value);
void RHT_Consume_Volatile(double currentValue);

void RHT_Replication_Init(int numThreads);
void RHT_Replication_Finish();

//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// GLOBAL VARS ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

extern volatile long producerCount;
extern volatile long consumerCount;

extern pthread_t ** consumerThreads;
extern int NUM_CONSUMER_THREADS;

#endif //HPCCG_1_0_SYNCQUEUE_H
