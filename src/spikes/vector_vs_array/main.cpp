#include <chrono>
#include <iostream>
#include "functions.hpp"

using namespace std;

int main() {

    // std::vector performance

    int n = 1;
    for (int i = 0; i < 25; i++, n *= 2) {
        auto start = chrono::steady_clock::now();
        float sum = vector_summer(n);
        auto end = chrono::steady_clock::now();
        long time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
        cout
                << "i=" << i
                << ",n=" << n
                << ",kind=std::vector"
                << ",time=" << time_us << "µs"
                << ",time/n=" << time_us / (float)n
                << endl;
    }
}