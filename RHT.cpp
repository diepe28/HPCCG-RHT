//
// Created by diego on 30/10/17.
//

#include "RHT.h"

SyncQueue globaQueue;

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

SyncQueue SyncQueue_Init(){
    SyncQueue syncQueue;
    int i = 0;
    syncQueue.content = (double*) (malloc(sizeof(double) * SYNC_QUEUE_SIZE));
    for(; i < SYNC_QUEUE_SIZE; i++){
        syncQueue.content[i] = ALREADY_CONSUMED;
    }
    syncQueue.volatileValue = ALREADY_CONSUMED;
    syncQueue.checkState = 1;
    syncQueue.enqPtr = syncQueue.deqPtr = 0;
    return syncQueue;
}

void SyncQueue_Produce_Simple(double value){
    int nextEnq = (globaQueue.enqPtr + 1) % SYNC_QUEUE_SIZE;

    while(globaQueue.content[nextEnq] != ALREADY_CONSUMED){
        asm("pause");
    }

    globaQueue.content[globaQueue.enqPtr] = value;
    globaQueue.enqPtr = nextEnq;
    producerCount++;
}

void SyncQueue_Produce_Volatile(double value){

    globaQueue.volatileValue = value;
    globaQueue.checkState = 0;
    while (globaQueue.checkState == 0) {
        asm("pause");
    }

    // After this, the value has been validated
}

double SyncQueue_Consume() {
    double value = globaQueue.content[globaQueue.deqPtr];

    if (value == ALREADY_CONSUMED) {
        //printf("There was a des sync of the queue \n");
        do {
            asm("pause");
        } while (globaQueue.content[globaQueue.deqPtr] == ALREADY_CONSUMED);
    }

    value = globaQueue.content[globaQueue.deqPtr];
    globaQueue.content[globaQueue.deqPtr] = ALREADY_CONSUMED;
    globaQueue.deqPtr = (globaQueue.deqPtr + 1) % SYNC_QUEUE_SIZE;
    consumerCount++;
    return value;
}

void SyncQueue_Consume_Check(double currentValue){
    double otherValue = globaQueue.content[globaQueue.deqPtr];

    if (currentValue != otherValue) {
        // des-sync of the queue
        if (otherValue == ALREADY_CONSUMED) {
            //consumerCount++;
            //printf("There was a des sync of the queue \n");
            do {
                asm("pause");
            } while (globaQueue.content[globaQueue.deqPtr] == ALREADY_CONSUMED);

            otherValue = globaQueue.content[globaQueue.deqPtr];

            if (currentValue == otherValue){
                globaQueue.content[globaQueue.deqPtr] = ALREADY_CONSUMED;
                globaQueue.deqPtr = (globaQueue.deqPtr + 1) % SYNC_QUEUE_SIZE;
                return;
            }
        }
        printf("\n\n SOFT ERROR DETECTED, Consumer: %f vs Producer: %f PCount: %ld -- CCount: %ld, diff: %ld \n",
               currentValue, otherValue, producerCount, consumerCount, producerCount - consumerCount);
        exit(1);
    }else{
        globaQueue.content[globaQueue.deqPtr] = ALREADY_CONSUMED;
        globaQueue.deqPtr = (globaQueue.deqPtr + 1) % SYNC_QUEUE_SIZE;
        consumerCount++;
    }
}

void SyncQueue_Consume_Volatile(double currentValue) {

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

void createConsumerThreads(int numThreads){
    int i;

    NUM_CONSUMER_THREADS = numThreads;
    consumerThreads = (pthread_t **) malloc(sizeof(pthread_t*) * NUM_CONSUMER_THREADS);

    for (i = 0; i < NUM_CONSUMER_THREADS; i++)
        consumerThreads[i] = (pthread_t*) malloc(sizeof(pthread_t));
}

void Replication_Init(int numThreads){
    createConsumerThreads(numThreads);
    globaQueue = SyncQueue_Init();
}

void Replication_Finish(){
    if(globaQueue.content)
        free( (void*) globaQueue.content);
    int i = 0;

    for (i = 0; i < NUM_CONSUMER_THREADS; i++)
        free(consumerThreads[i]);
    free(consumerThreads);
}