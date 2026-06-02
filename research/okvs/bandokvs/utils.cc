#ifndef BANDOKVS_UTILS_H_
#define BANDOKVS_UTILS_H_

#include <chrono>
#include <iostream>

#include "utils.h"

namespace band_okvs {

auto start = std::chrono::high_resolution_clock::now();
auto end = std::chrono::high_resolution_clock::now();

void StartTimer() {
  start = std::chrono::high_resolution_clock::now();
}

void EndTimer() {
  end = std::chrono::high_resolution_clock::now();
}

std::chrono::duration<double> GetElapsedTime() {
  return end - start;
}

void PrintElapsedTime() {
  std::cout << GetElapsedTime().count() << " seconds" << std::endl;
}

}

#endif //BANDOKVS_UTILS_H_