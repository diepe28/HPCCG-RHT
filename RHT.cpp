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
#ifdef APPROACH_USING_POINTERS
    UsingPointers_Produce(value);
#endif
#ifdef APPROACH_ALREADY_CONSUMED
    AlreadyConsumed_Produce(value);
#endif
#ifdef APPROACH_NEW_LIMIT
    AlreadyConsumed_Produce(value);
#endif
#ifdef APPROACH_WRITE_INVERTED_NEW_LIMIT
    WriteLineInverted_Produce_Secure(value);
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif
}

void RHT_Produce(double value) {
#ifdef APPROACH_USING_POINTERS
    UsingPointers_Produce(value);
#endif
#ifdef APPROACH_ALREADY_CONSUMED
    AlreadyConsumed_Produce(value);
#endif
#ifdef APPROACH_NEW_LIMIT
    NewLimit_Produce(value);
#endif
#ifdef APPROACH_WRITE_INVERTED_NEW_LIMIT
    WriteInvertedNewLimit_Produce(value);
#endif
#ifdef APPROACH_MOODY_CAMEL
    MoodyCamel_Produce(value);
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif
}

void RHT_Consume_Check(double currentValue) {
#ifdef APPROACH_USING_POINTERS
    UsingPointers_Consume_Check(currentValue);
#endif
#ifdef APPROACH_ALREADY_CONSUMED
    AlreadyConsumed_Consume_Check(currentValue);
#endif
#ifdef APPROACH_NEW_LIMIT
    NewLimit_Consume_Check(currentValue);
#endif
#ifdef APPROACH_WRITE_INVERTED_NEW_LIMIT
    WriteInvertedNewLimit_Consume_Check(currentValue);
#endif
#ifdef APPROACH_MOODY_CAMEL
    MoodyCamel_Consume_Check(currentValue);
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif
}

double RHT_Consume() {
#ifdef APPROACH_USING_POINTERS
    return UsingPointers_Consume();
#endif
#ifdef APPROACH_ALREADY_CONSUMED
    return AlreadyConsumed_Consume();
#endif
#ifdef APPROACH_NEW_LIMIT
    return NewLimit_Consume();
#endif
#ifdef APPROACH_WRITE_INVERTED_NEW_LIMIT
    return WriteInvertedNewLimit_Consume();
#endif
#ifdef APPROACH_MOODY_CAMEL
    return MoodyCamel_Consume();
#else
    printf("NO APPROACH SPECIFIED\n");
    exit(1);
#endif

}
