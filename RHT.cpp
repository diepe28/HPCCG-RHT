//
// Created by diego on 30/10/17.
//

/// Try with different queues as well, so there are no castings necessary... will it be really faster?
/// Produce with no sync and consuming with no sync --> this needs the already consume approach
/// Produce with no sync, consuming with sync --> with just the pointers (to see of cost of writing alreadyConsumed)
/// Make them macros

#include "RHT.h"

RHT_Queue globaQueue;

int nextEnq, localDeq, newLimit;
extern double global_otherValue, global_thisValue;

volatile long producerCount;
volatile long consumerCount;

pthread_t ** consumerThreads;
int NUM_CONSUMER_THREADS;

void SetThreadAffinity(int threadId) {
    cpu_set_t cpuset;

    CPU_ZERO(&cpuset);
    CPU_SET(threadId, &cpuset);

    /* pin the thread to a core */
    if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset)) {
        fprintf(stderr, "Thread pinning failed!\n");
        exit(1);
    }
}

RHT_Queue Queue_Init(){
    RHT_Queue syncQueue;
    int i = 0;
    syncQueue.content = (double*) (malloc(sizeof(double) * RHT_QUEUE_SIZE));
    for(; i < RHT_QUEUE_SIZE; i++){
        syncQueue.content[i] = ALREADY_CONSUMED;
    }
    syncQueue.volatileValue = ALREADY_CONSUMED;
    syncQueue.checkState = 1;
    syncQueue.enqPtr = syncQueue.deqPtr = 0;
    return syncQueue;
}

void createConsumerThreads(int numThreads){
    int i;

    NUM_CONSUMER_THREADS = numThreads;
    consumerThreads = (pthread_t **) malloc(sizeof(pthread_t*) * NUM_CONSUMER_THREADS);

    for (i = 0; i < NUM_CONSUMER_THREADS; i++)
        consumerThreads[i] = (pthread_t*) malloc(sizeof(pthread_t));
}

void RHT_Replication_Init(int numThreads){
    createConsumerThreads(numThreads);
    globaQueue = Queue_Init();
}

void RHT_Replication_Finish(){
    if(globaQueue.content)
        free( (void*) globaQueue.content);
    int i = 0;

    for (i = 0; i < NUM_CONSUMER_THREADS; i++)
        free(consumerThreads[i]);
    free(consumerThreads);
}

//////////// QUEUE METHODS //////////////////
static inline void AlreadyConsumed_Produce(double value);
static inline double AlreadyConsumed_Consume();
static inline void AlreadyConsumed_Consume_Check(double currentValue);

static inline void UsingPointers_Produce(double value);
static inline double UsingPointers_Consume();
static inline void UsingPointers_Consume_Check(double currentValue);

static inline void NewLimit_Produce(double value);

// -------- Public Methods ----------

void RHT_Produce_Volatile(double value){
    globaQueue.volatileValue = value;
    globaQueue.checkState = 0;
    while (globaQueue.checkState == 0) {
        asm("pause");
    }

    // After this, the value has been validated
}

void RHT_Consume_Volatile(double currentValue){
    while (globaQueue.checkState == 1) {
        asm("pause");
    }

    if (currentValue != globaQueue.volatileValue){
        printf("\n\n SOFT ERROR DETECTED AT VOLATILE ACCESS PARAMETER, %f vs %f Producer: %ld -- Consumer: %ld, diff: %ld \n",
               currentValue, globaQueue.volatileValue, producerCount, consumerCount, producerCount - consumerCount);
        exit(1);
    }

    globaQueue.checkState = 1;

    // After this, the value has been validated
}

void RHT_Produce(double value){
//    AlreadyConsumed_Produce(value);
//    UsingPointers_Produce(value);
    NewLimit_Produce(value);
}

void RHT_Consume_Check(double currentValue){
    AlreadyConsumed_Consume_Check(currentValue);
    //UsingPointers_Consume_Check(currentValue);
}

//double RHT_Consume() {
//    AlreadyConsumed_Consume();
////    UsingPointers_Consume();
//}

// -------- Already Consumed Approach ----------

static inline void AlreadyConsumed_Produce(double value){
    int nextEnq = (globaQueue.enqPtr + 1) % RHT_QUEUE_SIZE;

    while(globaQueue.content[nextEnq] != ALREADY_CONSUMED){
        asm("pause");
    }

    globaQueue.content[globaQueue.enqPtr] = value;
    globaQueue.enqPtr = nextEnq;
    //producerCount++;
}

static inline double AlreadyConsumed_Consume() {
    double value = globaQueue.content[globaQueue.deqPtr];

    if (value == ALREADY_CONSUMED) {
        //printf("There was a des sync of the queue \n");
        do {
            asm("pause");
        } while (globaQueue.content[globaQueue.deqPtr] == ALREADY_CONSUMED);
        value = globaQueue.content[globaQueue.deqPtr];
    }

    globaQueue.content[globaQueue.deqPtr] = ALREADY_CONSUMED;
    globaQueue.deqPtr = (globaQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    //consumerCount++;
    return value;
}

static inline void AlreadyConsumed_Consume_Check(double currentValue){
    double otherValue = globaQueue.content[globaQueue.deqPtr];

    if (currentValue != otherValue) {
        // des-sync of the queue
        if (otherValue == ALREADY_CONSUMED) {
            //consumerCount++;
            //printf("There was a des sync of the queue \n");
            do asm("pause"); while (globaQueue.content[globaQueue.deqPtr] == ALREADY_CONSUMED);

            otherValue = globaQueue.content[globaQueue.deqPtr];

            if (currentValue == otherValue){
                globaQueue.content[globaQueue.deqPtr] = ALREADY_CONSUMED;
                globaQueue.deqPtr = (globaQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
                return;
            }
        }
        printf("\n\n SOFT ERROR DETECTED, Consumer: %f vs Producer: %f PCount: %ld -- CCount: %ld, diff: %ld \n",
               currentValue, otherValue, producerCount, consumerCount, producerCount - consumerCount);
        exit(1);
    }else{
        globaQueue.content[globaQueue.deqPtr] = ALREADY_CONSUMED;
        globaQueue.deqPtr = (globaQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
//        consumerCount++;
    }
}

// -------- Using Pointers Approach ----------

static inline void UsingPointers_Produce(double value) {
    int nextEnq = (globaQueue.enqPtr + 1) % RHT_QUEUE_SIZE;

    while (nextEnq == globaQueue.deqPtr){
        asm("pause");
    }

    globaQueue.content[globaQueue.enqPtr] = value;
    globaQueue.enqPtr = nextEnq;
    //producerCount++;
}

static inline double UsingPointers_Consume(){

    while(globaQueue.deqPtr == globaQueue.enqPtr){
        asm("pause");
    }

    double value = globaQueue.content[globaQueue.deqPtr];
    globaQueue.deqPtr = (globaQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    //consumerCount++;
    return value;
}

static inline void UsingPointers_Consume_Check(double currentValue){
    while(globaQueue.deqPtr == globaQueue.enqPtr){
        asm("pause");
    }

    if(globaQueue.content[globaQueue.deqPtr] != currentValue){
        printf("\n\n SOFT ERROR DETECTED, Consumer: %f vs Producer: %f PCount: %ld -- CCount: %ld, diff: %ld \n",
               currentValue, globaQueue.content[globaQueue.deqPtr], producerCount, consumerCount, producerCount - consumerCount);
        exit(1);
    }

    globaQueue.deqPtr = (globaQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
    //consumerCount++;
}

// -------- New Limit Approach ----------

static inline void NewLimit_Produce(double value){
    globaQueue.content[globaQueue.enqPtr] = value;
    globaQueue.enqPtr = (globaQueue.enqPtr + 1) % RHT_QUEUE_SIZE;
}

void RHT_Produce_Secure(double value){
    int nextEnq = (globaQueue.enqPtr + 1) % RHT_QUEUE_SIZE;

    while(globaQueue.content[nextEnq] != ALREADY_CONSUMED){
        asm("pause");
    }

    globaQueue.content[globaQueue.enqPtr] = value;
    globaQueue.enqPtr = nextEnq;
    //producerCount++;
}

static inline void NewLimit_Consume(){
    AlreadyConsumed_Consume();
}

static inline void NewLimit_Consume_Check(double currentValue){
    AlreadyConsumed_Consume_Check(currentValue);
}