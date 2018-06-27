//
// Created by diego on 10/11/17.
//

#include "RHT.h"

RHT_QUEUE globalQueue;

#ifdef APPROACH_WANG
WANG_QUEUE wangQueue;
#endif

long producerCount;

int wait_var;
double wait_calc;

double groupVarProducer;
double groupVarConsumer;
int groupIncompleteConsumer;
int groupIncompleteProducer;
__thread long iterCountProducer;
__thread long iterCountConsumer;

long consumerCount;
int printValues = 0;

void RHT_Produce(double value) {
#if APPROACH_USING_POINTERS == 1
    UsingPointers_Produce(value);
#elif APPROACH_ALREADY_CONSUMED == 1 || APPROACH_CONSUMER_NO_SYNC == 1 // consumer no sync uses alreadyConsume produce method
    AlreadyConsumed_Produce(value);
#elif APPROACH_NEW_LIMIT == 1
    AlreadyConsumed_Produce(value);
#elif APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
    WriteInverted_Produce_Secure(value);
#elif VAR_GROUPING == 1 && (APPROACH_WANG == 1 || APPROACH_MIX_WANG == 1) // var grouping for wang and mix
    VG_Produce(value);
#elif APPROACH_WANG == 1
    Wang_Produce(value);
#elif APPROACH_MIX_WANG == 1
    Mix_Produce(value);
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif
}

// directly pushes a new value in the queue (regardless of var grouping)
void RHT_Produce_NoCheck(double value) {
#if APPROACH_MIX_WANG == 1
    Mix_Produce(value);
#elif APPROACH_WANG == 1
    Wang_Produce(value);
#else
    RHT_Produce(value);
#endif
}

void RHT_Consume_Check(double currentValue) {
#if APPROACH_USING_POINTERS == 1
    UsingPointers_Consume_Check(currentValue);
#elif APPROACH_ALREADY_CONSUMED == 1
    AlreadyConsumed_Consume_Check(currentValue);
#elif APPROACH_CONSUMER_NO_SYNC == 1 || APPROACH_NEW_LIMIT == 1 || APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
    NoSyncConsumer_Consume_Check(currentValue); // they all use the no sync consumer
#elif VAR_GROUPING == 1 && (APPROACH_WANG == 1 || APPROACH_MIX_WANG == 1) // var grouping for wang and mix
    VG_Consume_Check(currentValue);
#elif APPROACH_WANG == 1
    Wang_Consume_Check(currentValue);
#elif APPROACH_MIX_WANG == 1
    Mix_Consume_Check(currentValue);
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif
}

double RHT_Consume() {
#if APPROACH_USING_POINTERS == 1
    return UsingPointers_Consume();
#elif APPROACH_ALREADY_CONSUMED == 1
    return AlreadyConsumed_Consume();
#elif APPROACH_CONSUMER_NO_SYNC == 1 || APPROACH_NEW_LIMIT == 1 || APPROACH_WRITE_INVERTED_NEW_LIMIT == 1
    return NoSyncConsumer_Consume(); // they all use the no sync consumer
#elif APPROACH_WANG == 1
    return Wang_Consume();
#elif APPROACH_MIX_WANG == 1
    return Mix_Consume();
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif

}
