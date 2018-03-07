    //
// Created by diego on 30/10/17.
//

// -D CMAKE_C_COMPILER=/usr/bin/clang-5.0 -D CMAKE_CXX_COMPILER=/usr/bin/clang++-5.0
// -D CMAKE_C_COMPILER=/usr/bin/gcc-7 -D CMAKE_CXX_COMPILER=/usr/bin/g++-7
// gcc -O3 -Q --help=optimizers, to list the optimizations enabled


#ifndef HPCCG_1_0_SYNCQUEUE_H
#define HPCCG_1_0_SYNCQUEUE_H


#include <cstdlib>
#include <cstdio>
#include <pthread.h>
#include <zconf.h>
#include <cmath>
#include "readerwriterqueue.h"

using namespace moodycamel;

// -------- Macros ----------
// a cache line is 64 bytes, int -> 4 bytes, double -> 8 bytes
// so 8 doubles are 64 bytes, 16 doubles are 2 caches lines (for prefetcher)
#define CACHE_LINE_SIZE 16
#define RHT_QUEUE_SIZE 1024 // 512 is practically the same
#define MIN_PTR_DIST 300 // > 200 makes no real diff
#define ALREADY_CONSUMED -25802.89123

#define INLINE inline extern __attribute__((always_inline))

typedef int v4si __attribute__ ((vector_size (4*sizeof(int))));
typedef float v4sf __attribute__ ((vector_size (4*sizeof(float))));

#define INLINE inline extern __attribute__((always_inline))
    INLINE int equal(v4sf v1, v4sf v2) {
#if defined(__ALTIVEC__)
        return vec_all_eq((vector int)v1,(vector int)v2);
#elif defined(__SSE__)
        v4sf compare = __builtin_ia32_cmpeqps(v1,v2);
        return __builtin_ia32_movmskps(compare);
#else
        int * s1 = (int*)&v1;
int * s2 = (int*)&v2;
return (s1[0] == s2[0] && s1[1] == s2[1]
&& s1[2] == s2[2] && s1[3] == s2[3]);
#endif
    }


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
extern ReaderWriterQueue<double> moodyCamelQueue;

//extern int globalQueue.nextEnq, globalQueue.localDeq, globalQueue.newLimit, globalQueue.diff;
//extern double globalQueue.otherValue, globalQueue.thisValue;

static pthread_t **consumerThreads;
static int consumerThreadCount;

extern long producerCount;
extern long consumerCount;

#define EPSILON 0.00001

#define fequal(a,b) fabs(a-b) < EPSILON

/// ((RHT_QUEUE_SIZE-globalQueue.enqPtr + globalQueue.localDeq) % RHT_QUEUE_SIZE)-1; this would be faster and
/// with the -1 we make sure that the producer never really catches up to the consumer, but it might still happen
/// the other way around,

/// Another thing that we may try is inspired on wait-free programming from moodyCamel queue: the code is set up such
/// that each thread can always advance regardless of what the other is doing, when the queue is full the producer
/// allocates more memory

/// Also, we may try to implement the globalQueue.diff without spinning using globalQueue.deqPtr to avoid cache trashing, like calculate based
/// solely on the globalQueue.enqPtr the queue entry that should be != than ALREADY_CONSUMED, but then we might me producing less
/// that can actually be produced, because that entry might not be the last one with ALREADY_CONSUMED

/// TODO the first time it actually needs to sping while is < min(MIN_PTR_DIST, numIters), also the first time we should not yield
/// the processor without asking, the first time of each do while

#define calc_new_distance(waitValue)                                            \
    globalQueue.localDeq = globalQueue.deqPtr;                                  \
    waitValue = (globalQueue.enqPtr >= globalQueue.localDeq ?                   \
                (RHT_QUEUE_SIZE - globalQueue.enqPtr) + globalQueue.localDeq:   \
                globalQueue.localDeq - globalQueue.enqPtr)-1;

#define calc_and_move_normal(operation, value)                          \
    operation;                                                          \
    globalQueue.content[globalQueue.enqPtr] = value;                    \
    globalQueue.enqPtr = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;

#define calc_and_move_inverted(operation, value)                        \
    operation;                                                          \
    globalQueue.nextEnq = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;    \
    globalQueue.content[globalQueue.nextEnq] = ALREADY_CONSUMED;        \
    asm volatile("" ::: "memory");                                      \
    globalQueue.content[globalQueue.enqPtr] = value;                    \
    globalQueue.enqPtr = globalQueue.nextEnq;

#if APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
#define calc_and_move(operation, value) calc_and_move_inverted(operation, value)
#else
#define calc_and_move(operation, value) calc_and_move_normal(operation, value)
#endif

// waits until the distance between pointers is > MIN_PTR_DIST
#define wait_enough_distance(waitValue)         \
    calc_new_distance(waitValue)                \
    if(waitValue < MIN_PTR_DIST){               \
        /*producerCount++;*/                    \
        do{                                     \
            asm("pause");                       \
            asm("pause");                       \
            calc_new_distance(waitValue)        \
        }while(waitValue < MIN_PTR_DIST);       \
    }

// waits until the nextEnqIndex = enq+MIN_PTR_DIST, is ALREADY_CONSUMED
#define wait_next_enqIndex(waitValue)                                                   \
    waitValue = MIN_PTR_DIST;                                                           \
    globalQueue.nextEnq = (globalQueue.enqPtr + MIN_PTR_DIST) % RHT_QUEUE_SIZE;         \
    if(!fequal(ALREADY_CONSUMED, globalQueue.content[globalQueue.nextEnq])){            \
        /*producerCount++;*/                                                            \
        do{                                                                             \
            asm("pause");                                                               \
            asm("pause");                                                               \
            asm("pause");                                                               \
        }while(!fequal(ALREADY_CONSUMED, globalQueue.content[globalQueue.nextEnq]));    \
    }

#if APPROACH_NEW_ENQ_INDEX == 1
#define wait_for_consumer(waitValue) wait_next_enqIndex(waitValue)
#else
#define wait_for_consumer(waitValue) wait_enough_distance(waitValue)
#endif

#define replicate_loop_for(numIters, iterator, value, operation)    \
    wait_for_consumer(globalQueue.newLimit)                         \
    iterator = 0;                                                   \
    while (globalQueue.newLimit < numIters) {                       \
        for (; iterator < globalQueue.newLimit; iterator++){        \
            calc_and_move(operation, value)                         \
        }                                                           \
        wait_for_consumer(globalQueue.diff)                         \
        globalQueue.newLimit += globalQueue.diff;                   \
    }                                                               \
    for (; iterator < numIters; iterator++){                        \
        calc_and_move(operation, value)                             \
    }

#define Report_Soft_Error(consumerValue, producerValue) \
    printf("\n SOFT ERROR DETECTED, Consumer: %f Producer: %f -- PCount: %ld , CCount: %ld\n",  \
            consumerValue, producerValue, producerCount, consumerCount); \
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
    globalQueue.enqPtr = globalQueue.deqPtr = consumerCount = producerCount = 0;
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
}

static inline double AlreadyConsumed_Consume() {
    double value = globalQueue.content[globalQueue.deqPtr];

#if BRANCH_HINT == 1
    if (__builtin_expect(fequal(value, ALREADY_CONSUMED), 0))
#else
    if (fequal(value, ALREADY_CONSUMED))
#endif
    {
#if COUNT_QUEUE_DESYNC == 1
         consumerCount++;
#endif
        do {
            asm("pause");
        } while (globalQueue.content[globalQueue.deqPtr] == ALREADY_CONSUMED);
        value = globalQueue.content[globalQueue.deqPtr];
    }

    globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;
    globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    return value;
}

static inline void AlreadyConsumed_Consume_Check(double currentValue) {
    globalQueue.otherValue = globalQueue.content[globalQueue.deqPtr];

#if BRANCH_HINT == 1
    if(__builtin_expect(fequal(currentValue, globalQueue.otherValue),1))
#else
     if(fequal(currentValue, globalQueue.otherValue))
#endif
    {
        globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;
        globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    } else {
        // des-sync of the queue
#if COUNT_QUEUE_DESYNC == 1
        consumerCount++;
#endif
#if BRANCH_HINT == 1
        if(__builtin_expect(fequal(globalQueue.otherValue, ALREADY_CONSUMED), 1))
#else
        if(fequal(globalQueue.otherValue, ALREADY_CONSUMED))
#endif
        {
            do asm("pause"); while (fequal(globalQueue.content[globalQueue.deqPtr], ALREADY_CONSUMED));
            globalQueue.otherValue = globalQueue.content[globalQueue.deqPtr];

#if BRANCH_HINT == 1
            if (__builtin_expect(fequal(currentValue, globalQueue.otherValue), 1))
#else
            if (fequal(currentValue, globalQueue.otherValue))
#endif
            {
                globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;
                globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
                return;
            }
        }

        Report_Soft_Error(currentValue, globalQueue.otherValue)
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
}

static inline double UsingPointers_Consume() {

    while (globalQueue.deqPtr == globalQueue.enqPtr) {
        asm("pause");
    }

    double value = globalQueue.content[globalQueue.deqPtr];
    globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
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
}

// -------- New Limit Approach ----------

static inline void NewLimit_Produce(double value) {
    globalQueue.content[globalQueue.enqPtr] = value;
    globalQueue.enqPtr = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;
}

static inline double NewLimit_Consume() {
    return AlreadyConsumed_Consume();
}

static inline void NewLimit_Consume_Check(double currentValue) {
    AlreadyConsumed_Consume_Check(currentValue);
}

// -------- Write Inverted New Limit Approach ----------

static inline void WriteInvertedNewLimit_Produce(double value) {
    globalQueue.nextEnq = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;
    globalQueue.content[globalQueue.nextEnq] = ALREADY_CONSUMED;
    asm volatile("" ::: "memory");
    globalQueue.content[globalQueue.enqPtr] = value;
    globalQueue.enqPtr = globalQueue.nextEnq;
}

static inline double WriteInvertedNewLimit_Consume() {
    double value = globalQueue.content[globalQueue.deqPtr];

#if BRANCH_HINT == 1
    if (__builtin_expect(fequal(value, ALREADY_CONSUMED), 0))
#else
    if (fequal(value, ALREADY_CONSUMED))
#endif
    {
#if COUNT_QUEUE_DESYNC == 1
        consumerCount++;
#endif
        do asm("pause"); while (fequal(globalQueue.content[globalQueue.deqPtr], ALREADY_CONSUMED));
        value = globalQueue.content[globalQueue.deqPtr];
    }

    globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    return value;
}

static inline void WriteInvertedNewLimit_Consume_Check(double currentValue) {
    globalQueue.otherValue = globalQueue.content[globalQueue.deqPtr];

#if BRANCH_HINT == 1
    if (__builtin_expect(fequal(currentValue, globalQueue.otherValue), 1))
#else
    if (fequal(currentValue, globalQueue.otherValue))
#endif
    {
        globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    }
    else{
        // des-sync of the queue
#if BRANCH_HINT == 1
        if (__builtin_expect(fequal(globalQueue.otherValue, ALREADY_CONSUMED), 1))
#else
        if (fequal(globalQueue.otherValue, ALREADY_CONSUMED))
#endif
        {
#if COUNT_QUEUE_DESYNC == 1
            consumerCount++;
#endif
            do {
                asm("pause");
            } while (fequal(globalQueue.content[globalQueue.deqPtr], ALREADY_CONSUMED));

#if BRANCH_HINT == 1
            if (__builtin_expect(fequal(currentValue, globalQueue.content[globalQueue.deqPtr]), 1))
#else
            if (fequal(currentValue, globalQueue.content[globalQueue.deqPtr]))
#endif
            {
                globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
                return;
            }
            printf("The second time is also not equal ... \n");
        }
        Report_Soft_Error(currentValue, globalQueue.otherValue)
    }
}

static inline void WriteInverted_Produce_Secure(double value){
    while (globalQueue.content[globalQueue.enqPtr] != ALREADY_CONSUMED) {
        asm("pause");
    }
    globalQueue.nextEnq = (globalQueue.enqPtr + 1) % RHT_QUEUE_SIZE;
    globalQueue.content[globalQueue.nextEnq] = ALREADY_CONSUMED;
    globalQueue.content[globalQueue.enqPtr] = value;
    globalQueue.enqPtr = globalQueue.nextEnq;
}

// -------- Moody Camel Approach ----------

static inline void MoodyCamel_Produce(double value) {
    moodyCamelQueue.enqueue(value);
}

static inline double MoodyCamel_Consume() {
    double result;
    while(!moodyCamelQueue.try_dequeue(result)){
        asm("pause");
    }
    return result;
}

static inline void MoodyCamel_Consume_Check(double currentValue) {
    double otherValue = MoodyCamel_Consume();

    if (otherValue != currentValue) {
        Report_Soft_Error(currentValue, otherValue)
    }

}

//////////// PUBLIC QUEUE METHODS //////////////////
extern void RHT_Produce_Secure(double value);
extern void RHT_Produce(double value);
extern void RHT_Consume_Check(double currentValue);
extern double RHT_Consume();

#endif //HPCCG_1_0_SYNCQUEUE_H
