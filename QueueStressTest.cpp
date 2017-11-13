//
// Created by diego on 13/11/17.
//

#include <thread>
#include "QueueStressTest.h"

using namespace moodycamel;


static void consumer_RHT(void * args) {
    SetThreadAffinity(2);

    int i, j;
    long result;

    for(i = 0; i < MAX_ROWS; i++){
        for(j = 0; j < MAX_COLS; j++){
            result = i + j;
            RHT_Consume_Check(result);
        }
    }
}

static void producer_RHT() {
    SetThreadAffinity(0);

    int i, j;
    long result;

    for(i = 0; i < MAX_ROWS; i++){

        replicate_forLoop_no_sync(MAX_COLS, j, result, result = i + j)

        //for(j = 0; j < MAX_COLS; j++){
            //result = i + j;
            //RHT_Produce(result);
        //}
    }
}

static void test_rht_queue(){
    RHT_Replication_Init(1);
    int err = pthread_create(consumerThreads[0], NULL, (void *(*)(void *)) consumer_RHT, (void *) NULL);
    if (err) {
        fprintf(stderr, "Failed to create thread %d\n", 1);
        exit(1);
    }

    producer_RHT();
    pthread_join(*consumerThreads[0], NULL);
    RHT_Replication_Finish();
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

BlockingReaderWriterQueue<double> q;
ReaderWriterQueue<double> q2(100);

static void consumer_cameron(void * args) {
    SetThreadAffinity(0);

    int i, j;
    double result, otherValue;

    for(i = 0; i < MAX_ROWS; i++) {
        for (j = 0; j < MAX_COLS; j++) {
            result = i + j;

          tryDequeue:
            // Fully-blocking:
            //q.wait_dequeue(otherValue);
            if(q2.try_dequeue(otherValue)) {
                if (result != otherValue) {
                    printf("Soft error encountered at iteration %d , %d\n", i, j);
                    exit(1);

                }
            }else{
                goto tryDequeue;
            }
        }
    }
}

static void producer_cameron() {
    SetThreadAffinity(2);

    int i, j;
    double result;

    for(i = 0; i < MAX_ROWS; i++) {
        for (j = 0; j < MAX_COLS; j++) {
            // Fully-blocking:
            result = i + j;
            //q.enqueue(result);
            q2.enqueue(result);
        }
    }
}


static void test_cameron_queue(){
    int i;
    consumerThreads = (pthread_t **) malloc(sizeof(pthread_t *) * 2);

    for (i = 0; i < consumerThreadCount; i++)
        consumerThreads[i] = (pthread_t *) malloc(sizeof(pthread_t));


    int err = pthread_create(consumerThreads[0], NULL, (void *(*)(void *)) consumer_cameron, (void *) NULL);
    if (err) {
        fprintf(stderr, "Failed to create thread %d\n", 1);
        exit(1);
    }

    producer_cameron();

    pthread_join(*consumerThreads[0], NULL);

    for (i = 0; i < consumerThreadCount; i++)
        free(consumerThreads[i]);
    free(consumerThreads);
}

void TestQueues(){

    double milliseconds_elapsed;
    struct timespec start, finish;
    int numRuns = 5, i;

    for(i = 0; i < numRuns; i++) {
        clock_gettime(CLOCK_MONOTONIC, &start);

        test_rht_queue();

        clock_gettime(CLOCK_MONOTONIC, &finish);
        milliseconds_elapsed = (finish.tv_sec - start.tv_sec + (finish.tv_nsec - start.tv_nsec) / 1000000000.0) * 1000;

        printf("RHT -- Miliseconds elapsed %f \n", milliseconds_elapsed);
    }

    printf("\n --------------------------------------\n\n");

    for(i = 0; i < numRuns; i++) {
        clock_gettime(CLOCK_MONOTONIC, &start);

        test_cameron_queue();

        clock_gettime(CLOCK_MONOTONIC, &finish);
        milliseconds_elapsed = (finish.tv_sec - start.tv_sec + (finish.tv_nsec - start.tv_nsec) / 1000000000.0) * 1000;

        printf("Cameron -- Miliseconds elapsed %f \n", milliseconds_elapsed);
    }
}