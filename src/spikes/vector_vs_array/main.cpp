#include <chrono>
#include <iostream>
#include "functions.hpp"

#define NUMBER_OF_LOOPS 22

using namespace std;

void vector_test() {
    int n = 1;
    for (int i = 0; i < NUMBER_OF_LOOPS; i++, n *= 2) {
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

void array_test() {
    int n = 1;
    for (int i = 0; i < NUMBER_OF_LOOPS; i++, n *= 2) {
        auto start = chrono::steady_clock::now();
        float sum = array_summer(n);
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

int main() {

    auto start = chrono::steady_clock::now();
    vector_test();
    auto end = chrono::steady_clock::now();
    long time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "vector total time: " << time_us << "\n\n";

    start = chrono::steady_clock::now();
    array_test();
    end = chrono::steady_clock::now();
    time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "array total time: " << time_us << "\n\n";
}

