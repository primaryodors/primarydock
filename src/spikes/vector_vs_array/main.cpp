#include <chrono>
#include <iostream>
#include "functions.hpp"

#define NUMBER_OF_LOOPS 22

using namespace std;

void vector_traditional_forloop_test() {
    int n = 1;
    for (int i = 0; i < NUMBER_OF_LOOPS; i++, n *= 2) {
        auto start = chrono::steady_clock::now();
        float sum = vector_traditional_forloop(n);
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

void vector_ranged_forloop_test() {
    int n = 1;
    for (int i = 0; i < NUMBER_OF_LOOPS; i++, n *= 2) {
        auto start = chrono::steady_clock::now();
        float sum = vector_ranged_forloop(n);
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
                << ",kind=array"
                << ",time=" << time_us << "µs"
                << ",time/n=" << time_us / (float)n
                << endl;
    }
}

int main() {

    auto start = chrono::steady_clock::now();
    vector_traditional_forloop_test();
    auto end = chrono::steady_clock::now();
    long time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "vector traditional forloop total time: " << time_us << "\n\n";

    start = chrono::steady_clock::now();
    vector_ranged_forloop_test();
    end = chrono::steady_clock::now();
    time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "vector ranged forloop total time: " << time_us << "\n\n";

    start = chrono::steady_clock::now();
    array_test();
    end = chrono::steady_clock::now();
    time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "array total time: " << time_us << "\n\n";
}

