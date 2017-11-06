//
// Created by diego on 30/10/17.
//

#ifndef HPCCG_1_0_SYNCQUEUE_H
#define HPCCG_1_0_SYNCQUEUE_H

#include <cstdlib>
#include <cstdio>
#include <pthread.h>

// -------- Macros ----------

#define RHT_QUEUE_SIZE 1024
#define ALREADY_CONSUMED -2

#define replicate_newLimit(value, operation)                                    \
    do {                                                                        \
        for (i = 0; i < newLimit; i++){                                         \
            value = operation;                                                  \
            RHT_Produce(value);                                                 \
        }                                                                       \
        nextEnq = (globaQueue.enqPtr + 1) % RHT_QUEUE_SIZE;                     \
        while (globaQueue.content[nextEnq] != ALREADY_CONSUMED) asm("pause");   \
        localDeq = globaQueue.deqPtr;                                           \
        newLimit += globaQueue.enqPtr >= localDeq ?                             \
                    (RHT_QUEUE_SIZE - globaQueue.enqPtr) + localDeq :           \
                    localDeq - globaQueue.enqPtr;                               \
    } while (newLimit < n);                                                     \
    for (; i < n; i++){                                                         \
        operation;                                                              \
        RHT_Produce(value);                                                     \
    }


#define Macro_AlreadyConsumed_Produce(value)                    \
    nextEnq = (globaQueue.enqPtr + 1) % RHT_QUEUE_SIZE;         \
    while(globaQueue.content[nextEnq] != ALREADY_CONSUMED){     \
        asm("pause");                                           \
    }                                                           \
    globaQueue.content[globaQueue.enqPtr] = value;              \
    globaQueue.enqPtr = nextEnq;                                \
    /*producerCount++;*/


#define Macro_AlreadyConsumed_Consume(value) {                      \
    value = globaQueue.content[globaQueue.deqPtr];                  \
    if (value == ALREADY_CONSUMED) {                                \
        do asm("pause"); while (globaQueue.content[globaQueue.deqPtr] == ALREADY_CONSUMED); \
        value = globaQueue.content[globaQueue.deqPtr];              \
    }                                                               \
                                                                    \
    globaQueue.content[globaQueue.deqPtr] = ALREADY_CONSUMED;       \
    globaQueue.deqPtr = (globaQueue.deqPtr + 1) % RHT_QUEUE_SIZE;   \
    /*consumerCount++;*/

#define Macro_Report_Soft_Error(consumerValue, producerValue) \
    printf("\n\n SOFT ERROR DETECTED, Consumer: %f vs Producer: %f PCount: %ld -- CCount: %ld, diff: %ld \n", \
            consumerValue, producerValue, producerCount, consumerCount, producerCount - consumerCount); \
    exit(1);

#define Macro_AlreadyConsumed_Consume_Check(currentValue)                       \
    otherValue = globaQueue.content[globaQueue.deqPtr];                         \
    if (currentValue != otherValue) {                                           \
        /* des-sync of the queue */                                             \
        if (otherValue == ALREADY_CONSUMED) {                                   \
            /*consumerCount++;*/                                                \
            do asm("pause");while (globaQueue.content[globaQueue.deqPtr] == ALREADY_CONSUMED); \
            otherValue = globaQueue.content[globaQueue.deqPtr];                 \
                                                                                \
            if (currentValue == otherValue){                                    \
                globaQueue.content[globaQueue.deqPtr] = ALREADY_CONSUMED;       \
                globaQueue.deqPtr = (globaQueue.deqPtr + 1) % RHT_QUEUE_SIZE;   \
            }else                                                               \
                 Macro_Report_Soft_Error(otherValue, currentValue)              \
        }else                                                                   \
            Macro_Report_Soft_Error(otherValue, currentValue)                   \
    }else{                                                                      \
        globaQueue.content[globaQueue.deqPtr] = ALREADY_CONSUMED;               \
        globaQueue.deqPtr = (globaQueue.deqPtr + 1) % RHT_QUEUE_SIZE;           \
        /*consumerCount++;*/                                                    \
    }

#define RHT_Produce(value) \
    Macro_AlreadyConsumed_Produce(value)

//#define RHT_Consume(value) \
//    Macro_AlreadyConsumed_Consume(value)
//
//#define RHT_Consume_Check(currentValue) \
//   Macro_AlreadyConsumed_Consume_Check


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

//void RHT_Produce(double value);
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
extern int nextEnq;

extern volatile long producerCount;
extern volatile long consumerCount;

extern pthread_t ** consumerThreads;
extern int NUM_CONSUMER_THREADS;

#endif //HPCCG_1_0_SYNCQUEUE_H
