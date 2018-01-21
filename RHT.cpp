//
// Created by diego on 10/11/17.
//

#include "RHT.h"

RHT_QUEUE globalQueue;
ReaderWriterQueue<double> moodyCamelQueue(256);

long producerCount;
long consumerCount;

//int nextEnq, localDeq, newLimit, diff;
//// Consumer Local
//double otherValue, thisValue;


void RHT_Produce_Secure(double value) {
#if APPROACH_USING_POINTERS == 1
    UsingPointers_Produce(value);
#elif APPROACH_ALREADY_CONSUMED == 1
    AlreadyConsumed_Produce(value);
#elif APPROACH_NEW_LIMIT == 1
    AlreadyConsumed_Produce(value);
#elif APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
    WriteInverted_Produce_Secure(value);
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif
}

void RHT_Produce(double value) {
#if APPROACH_USING_POINTERS == 1
    UsingPointers_Produce(value);
#elif APPROACH_ALREADY_CONSUMED == 1
    AlreadyConsumed_Produce(value);
#elif APPROACH_NEW_LIMIT == 1
    NewLimit_Produce(value);
#elif APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
    WriteInvertedNewLimit_Produce(value);
#elif APPROACH_MOODY_CAMEL == 1
    MoodyCamel_Produce(value);
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif
}

void RHT_Consume_Check(double currentValue) {
#if APPROACH_USING_POINTERS == 1
    UsingPointers_Consume_Check(currentValue);
#elif APPROACH_ALREADY_CONSUMED == 1
    AlreadyConsumed_Consume_Check(currentValue);
#elif APPROACH_NEW_LIMIT == 1
    NewLimit_Consume_Check(currentValue);
#elif APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
    WriteInvertedNewLimit_Consume_Check(currentValue);
#elif APPROACH_MOODY_CAMEL == 1
    MoodyCamel_Consume_Check(currentValue);
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif
}

void RHT_Consume_CheckSpecial(double currentValue, int rank) {

    globalQueue.otherValue = globalQueue.content[globalQueue.deqPtr];

    if(fequal(currentValue, globalQueue.otherValue)){
        globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;
        globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
        consumerCount++;
    } else {
        // des-sync of the queue
        if(fequal(globalQueue.otherValue, ALREADY_CONSUMED)){
            //consumerCount++;
            do asm("pause"); while (fequal(globalQueue.content[globalQueue.deqPtr], ALREADY_CONSUMED));

            globalQueue.otherValue = globalQueue.content[globalQueue.deqPtr];

            if(fequal(currentValue, globalQueue.otherValue)){
                globalQueue.content[globalQueue.deqPtr] = ALREADY_CONSUMED;
                globalQueue.deqPtr = (globalQueue.deqPtr + 1) % RHT_QUEUE_SIZE;
                return;
            }
        }
        printf("The rank %d is the one with the problem\n",rank);
        Report_Soft_Error(currentValue, globalQueue.otherValue)
    }

}

double RHT_Consume() {
#if APPROACH_USING_POINTERS == 1
    return UsingPointers_Consume();
#elif APPROACH_ALREADY_CONSUMED == 1
    return AlreadyConsumed_Consume();
#elif APPROACH_NEW_LIMIT == 1
    return NewLimit_Consume();
#elif APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
    return WriteInvertedNewLimit_Consume();
#elif APPROACH_MOODY_CAMEL == 1
    return MoodyCamel_Consume();
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif

}
