    //
// Created by diego on 30/10/17.
//

#ifndef HPCCG_1_0_SYNCQUEUE_H
#define HPCCG_1_0_SYNCQUEUE_H

#include <cstdlib>
#include <cstdio>
#include <pthread.h>

// -------- Macros ----------
// a cache line is 64 bytes, int -> 4 bytes, double -> 8 bytes
// so 8 doubles are 64 bytes, 16 doubles are 2 caches lines (for prefetcher)
#define CACHE_LINE_SIZE 16
#define RHT_QUEUE_SIZE 1024 // 512 is practically the same
#define MIN_PTR_DIST 200 // > 200 makes no real diff
#define ALREADY_CONSUMED -2

typedef struct { ;
    volatile int deqPtr;
    double padding0[CACHE_LINE_SIZE - sizeof(int)];
    volatile int enqPtr;
    double padding1[CACHE_LINE_SIZE - sizeof(int)];
    volatile double *content;
    double padding2[CACHE_LINE_SIZE - sizeof(double *)];
    volatile double volatileValue;
    volatile int checkState;
    double padding3[CACHE_LINE_SIZE - sizeof(int) - sizeof(double)];
    // Producer Local, does not make any diff compared to having them outside of the queue
    int nextEnq, localDeq, newLimit, diff;
    double padding4[CACHE_LINE_SIZE - 4 * sizeof(int)];
    // Consumer Local
    double otherValue, thisValue;
}RHT_QUEUE;


extern RHT_QUEUE globalQueue;

//extern int globalQueue.nextEnq, globalQueue.localDeq, globalQueue.newLimit, globalQueue.diff;
//extern double globalQueue.otherValue, globalQueue.thisValue;

static pthread_t **consumerThreads;
static int consumerThreadCount;

extern long producerCount;
extern long consumerCount;

/// ((RHT_QUEUE_SIZE-globalQueue.enqPtr + globalQueue.localDeq) % RHT_QUEUE_SIZE)-1; this would be faster and
/// with the -1 we make sure that the producer never really catches up to the consumer, but it might still happen
/// the other way around,
/// Also, we may try to implement the globalQueue.diff without spinning using globalQueue.deqPtr to avoid cache trashing, like calculate based
/// solely on the globalQueue.enqPtr the queue entry that should be != than ALREADY_CONSUMED, but then we might me producing less
/// that can actually be produced, because that entry might not be the last one with ALREADY_CONSUMED
    /// TODO the first time it actually needs to sping while is < min(MIN_PTR_DIST, numIters), also the first time we should not yield
    /// the processor without asking, the first time of each do while
//#define replicate_forLoop_no_sync(numIters, iterator, value, operation)             \
//    do{                                                                             \
//        asm("pause");                                                               \
//        globalQueue.localDeq = globalQueue.deqPtr;                                  \
//        globalQueue.newLimit = (globalQueue.enqPtr >= globalQueue.localDeq ?        \
//                   (RHT_QUEUE_SIZE - globalQueue.enqPtr) + globalQueue.localDeq :   \
//                   globalQueue.localDeq - globalQueue.enqPtr)-1;                    \
//    }                                                                               \
//    while(globalQueue.newLimit < MIN_PTR_DIST);                                     \
//    iterator = 0;                                                                   \
//    while (globalQueue.newLimit < numIters) {                                       \
//        for (; iterator < globalQueue.newLimit; iterator++){                        \
//            operation;                                                              \
//            globalQueue.content[globalQueue.enqPtr] = value;                        \
//            globalQueue.enqPtr = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;         \
//        }                                                                           \
//        do{                                                                         \
//            asm("pause");                                                           \
//            globalQueue.localDeq = globalQueue.deqPtr;                              \
//            globalQueue.diff = (globalQueue.enqPtr >= globalQueue.localDeq ?        \
//                   (RHT_QUEUE_SIZE - globalQueue.enqPtr) + globalQueue.localDeq :   \
//                   globalQueue.localDeq - globalQueue.enqPtr)-1;                    \
//        }                                                                           \
//        while(globalQueue.diff < MIN_PTR_DIST);                                     \
//        globalQueue.newLimit += globalQueue.diff;                                   \
//    }                                                                               \
//    for (; iterator < numIters; iterator++){                                        \
//        operation;                                                                  \
//        globalQueue.content[globalQueue.enqPtr] = value;                            \
//        globalQueue.enqPtr = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;             \
//    }

#define replicate_forLoop_no_sync(numIters, iterator, value, operation)             \
    globalQueue.localDeq = globalQueue.deqPtr;                                      \
    globalQueue.newLimit = (globalQueue.enqPtr >= globalQueue.localDeq ?            \
                   (RHT_QUEUE_SIZE - globalQueue.enqPtr) + globalQueue.localDeq :   \
                   globalQueue.localDeq - globalQueue.enqPtr)-1;                    \
    if(globalQueue.newLimit < MIN_PTR_DIST){                                        \
        do{                                                                         \
            asm("pause");                                                           \
            globalQueue.localDeq = globalQueue.deqPtr;                              \
            globalQueue.newLimit = (globalQueue.enqPtr >= globalQueue.localDeq ?    \
                   (RHT_QUEUE_SIZE - globalQueue.enqPtr) + globalQueue.localDeq :   \
                   globalQueue.localDeq - globalQueue.enqPtr) - 1;                  \
        }while(globalQueue.newLimit < MIN_PTR_DIST);                                \
    }                                                                               \
    iterator = 0;                                                                   \
    while (globalQueue.newLimit < numIters) {                                       \
        for (; iterator < globalQueue.newLimit; iterator++){                        \
            operation;                                                              \
            globalQueue.content[globalQueue.enqPtr] = value;                        \
            globalQueue.enqPtr = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;         \
        }                                                                           \
        globalQueue.localDeq = globalQueue.deqPtr;                                  \
            globalQueue.diff = (globalQueue.enqPtr >= globalQueue.localDeq ?        \
                   (RHT_QUEUE_SIZE - globalQueue.enqPtr) + globalQueue.localDeq :   \
                   globalQueue.localDeq - globalQueue.enqPtr)-1;                    \
        if(globalQueue.diff < MIN_PTR_DIST) {                                       \
            do{                                                                         \
                asm("pause");                                                           \
                globalQueue.localDeq = globalQueue.deqPtr;                              \
                globalQueue.diff = (globalQueue.enqPtr >= globalQueue.localDeq ?        \
                       (RHT_QUEUE_SIZE - globalQueue.enqPtr) + globalQueue.localDeq :   \
                    globalQueue.localDeq - globalQueue.enqPtr)-1;                    \
            }while(globalQueue.diff < MIN_PTR_DIST);                                \
        }                                                                           \
        globalQueue.newLimit += globalQueue.diff;                                   \
    }                                                                               \
    for (; iterator < numIters; iterator++){                                        \
        operation;                                                                  \
        globalQueue.content[globalQueue.enqPtr] = value;                            \
        globalQueue.enqPtr = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;             \
    }

#define Report_Soft_Error(consumerValue, producerValue) \
    printf("\n SOFT ERROR DETECTED, Consumer: %f Producer: %f -- PCount: %ld , CCount: %ld\n",  \
            consumerValue, producerValue, producerCount, consumerCount);                        \
    exit(1);

#define RHT_Produce_Volatile(volValue)                          \
    globalQueue.volatileValue = volValue;                       \
    globalQueue.checkState = 0;                                 \
    while (globalQueue.checkState == 0) asm("pause");

#define RHT_Consume_Volatile(volValue)                          \
    while (globalQueue.checkState == 1) asm("pause");           \
    if (volValue != globalQueue.volatileValue){                 \
        Report_Soft_Error(volValue, globalQueue.volatileValue)  \
    }                                                           \
    globalQueue.checkState = 1;

#define Macro_AlreadyConsumed_Produce(value)                  \
    globalQueue.nextEnq = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;      \
    while(globalQueue.content[globalQueue.nextEnq] != ALREADY_CONSUMED)   \
        asm("pause");                                         \
    /*producerCount++;*/                                      \
    globalQueue.content[globalQueue.enqPtr] = value;          \
    globalQueue.enqPtr = globalQueue.nextEnq;

#define Macro_AlreadyConsumed_Consume_Check(currentValue)                       \
    globalQueue.thisValue = (double) currentValue;                                          \
    globalQueue.otherValue = globalQueue.content[globalQueue.deqPtr];                         \
    if (globalQueue.thisValue != globalQueue.otherValue) {                                              \
        /* des-sync of the queue */                                             \
        if (globalQueue.otherValue == ALREADY_CONSUMED) {                                   \
            /*consumerCount++;*/                                                \
            do asm("pause"); while (globalQueue.content[globalQueue.deqPtr] == ALREADY_CONSUMED); \
            globalQueue.otherValue = globalQueue.content[globalQueue.deqPtr];                 \
                                                                                \
            if (globalQueue.thisValue == globalQueue.otherValue){                                       \
                globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;       \
                globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;   \
            }else{                                                              \
                 Report_Soft_Error(globalQueue.otherValue, globalQueue.thisValue)                 \
            }                                                                   \
        }else {                                                                 \
            Report_Soft_Error(globalQueue.otherValue, globalQueue.thisValue)                      \
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

//#define RHT_Produce(value) \
//    /*Macro_UsingPointers_Produce(value)*/ \
//    Macro_AlreadyConsumed_Produce(value)

//#define RHT_Consume_Check(currentValue) \
//   /*Macro_UsingPointers_Consume_Check(value)*/ \
//   Macro_AlreadyConsumed_Consume_Check(currentValue)

//#define RHT_Consume(value) \
//    /*Macro_UsingPointers_Consume(value)*/ \
//    Macro_AlreadyConsumed_Consume(value)


static void SetThreadAffinity(int threadId) {
    cpu_set_t cpuset;

    CPU_ZERO(&cpuset);
    CPU_SET(threadId, &cpuset);

    /* pin the thread to a core */
    if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset)) {
        fprintf(stderr, "Thread pinning failed!\n");
        exit(1);
    }
}

static void Queue_Init() {
    int i = 0;
    globalQueue.content = (double *) (malloc(sizeof(double) * RHT_QUEUE_SIZE));
    for (; i < RHT_QUEUE_SIZE; i++) {
        globalQueue.content[i] = ALREADY_CONSUMED;
    }
    globalQueue.volatileValue = ALREADY_CONSUMED;
    globalQueue.checkState = 1;
    globalQueue.enqPtr = globalQueue.deqPtr = 0;
}

static void createConsumerThreads(int numThreads) {
    int i;

    consumerThreadCount = numThreads;
    consumerThreads = (pthread_t **) malloc(sizeof(pthread_t *) * consumerThreadCount);

    for (i = 0; i < consumerThreadCount; i++)
        consumerThreads[i] = (pthread_t *) malloc(sizeof(pthread_t));
}

static void RHT_Replication_Init(int numThreads) {
    createConsumerThreads(numThreads);
    Queue_Init();
}

static void RHT_Replication_Finish() {
    if (globalQueue.content)
        free((void *) globalQueue.content);
    int i = 0;

    for (i = 0; i < consumerThreadCount; i++)
        free(consumerThreads[i]);
    free(consumerThreads);
}

//////////// INTERNAL QUEUE METHODS BODY //////////////////

// -------- Already Consumed Approach ----------

static inline void AlreadyConsumed_Produce(double value) {
    globalQueue.nextEnq = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;

    while (globalQueue.content[globalQueue.nextEnq] != ALREADY_CONSUMED) {
        asm("pause");
    }

    globalQueue.content[globalQueue.enqPtr] = value;
    globalQueue.enqPtr = globalQueue.nextEnq;
    //producerCount++;
}

static inline double AlreadyConsumed_Consume() {
    double value = globalQueue.content[globalQueue.deqPtr];

    if (value == ALREADY_CONSUMED) {
        //printf("There was a des sync of the queue \n");
        do {
            asm("pause");
        } while (globalQueue.content[globalQueue.deqPtr] == ALREADY_CONSUMED);
        value = globalQueue.content[globalQueue.deqPtr];
    }

    globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;
    globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    //consumerCount++;
    return value;
}

static inline void AlreadyConsumed_Consume_Check(double currentValue) {
    globalQueue.otherValue = globalQueue.content[globalQueue.deqPtr];

    if (currentValue != globalQueue.otherValue) {
        // des-sync of the queue
        if (globalQueue.otherValue == ALREADY_CONSUMED) {
            //consumerCount++;
            //printf("There was a des sync of the queue \n");
            do asm("pause"); while (globalQueue.content[globalQueue.deqPtr] == ALREADY_CONSUMED);

            globalQueue.otherValue = globalQueue.content[globalQueue.deqPtr];

            if (currentValue == globalQueue.otherValue) {
                globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;
                globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
                return;
            }
        }
        Report_Soft_Error(currentValue, globalQueue.otherValue)
    } else {
        globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;
        globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
//        consumerCount++;
    }
}

// -------- Using Pointers Approach ----------

static inline void UsingPointers_Produce(double value) {
    globalQueue.nextEnq = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;

    while (globalQueue.nextEnq == globalQueue.deqPtr) {
        asm("pause");
    }

    globalQueue.content[globalQueue.enqPtr] = value;
    globalQueue.enqPtr = globalQueue.nextEnq;
    //producerCount++;
}

static inline double UsingPointers_Consume() {

    while (globalQueue.deqPtr == globalQueue.enqPtr) {
        asm("pause");
    }

    double value = globalQueue.content[globalQueue.deqPtr];
    globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    //consumerCount++;
    return value;
}

static inline void UsingPointers_Consume_Check(double currentValue) {
    while (globalQueue.deqPtr == globalQueue.enqPtr) {
        asm("pause");
    }

    if (globalQueue.content[globalQueue.deqPtr] != currentValue) {
        Report_Soft_Error(currentValue, globalQueue.content[globalQueue.deqPtr])
    }

    globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    //consumerCount++;
}

// -------- New Limit Approach ----------

static inline void NewLimit_Produce(double value) {
    globalQueue.content[globalQueue.enqPtr] = value;
    globalQueue.enqPtr = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;
}

static inline double NewLimit_Consume() {
//    return UsingPointers_Consume();
    return AlreadyConsumed_Consume();
}

static inline void NewLimit_Consume_Check(double currentValue) {
//    UsingPointers_Consume_Check(currentValue);
    AlreadyConsumed_Consume_Check(currentValue);
}

//////////// PUBLIC QUEUE METHODS //////////////////
extern void RHT_Produce_Secure(double value);
extern void RHT_Produce(double value);
extern void RHT_Consume_Check(double currentValue);
extern double RHT_Consume();

#endif //HPCCG_1_0_SYNCQUEUE_H
