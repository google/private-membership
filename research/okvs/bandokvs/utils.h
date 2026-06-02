#ifndef BANDOKVS_UTILS_H_
#define BANDOKVS_UTILS_H_

#include <chrono>
#include <iostream>

namespace band_okvs {

void StartTimer();

void EndTimer();

std::chrono::duration<double> GetElapsedTime();

void PrintElapsedTime();

}

#endif //BANDOKVS_UTILS_H_
