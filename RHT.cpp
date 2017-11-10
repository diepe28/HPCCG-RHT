//
// Created by diego on 10/11/17.
//

#include "RHT.h"

RHT_QUEUE globalQueue;

long producerCount;
long consumerCount;

int nextEnq, localDeq, newLimit, diff;
// Consumer Local
double otherValue, thisValue;


void RHT_Produce_Secure(double value) {
//    UsingPointers_Produce(value);
    AlreadyConsumed_Produce(value);
}

void RHT_Produce(double value) {
//    AlreadyConsumed_Produce(value);
//    UsingPointers_Produce(value);
    NewLimit_Produce(value);
}

void RHT_Consume_Check(double currentValue) {
//    AlreadyConsumed_Consume_Check(currentValue);
//    UsingPointers_Consume_Check(currentValue);
    NewLimit_Consume_Check(currentValue);
}

double RHT_Consume() {
//    return AlreadyConsumed_Consume();
//    return UsingPointers_Consume();
    return NewLimit_Consume();
}