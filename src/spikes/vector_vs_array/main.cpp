#include <chrono>
#include <iostream>
#include "functions.hpp"

#define NUMBER_OF_LOOPS 20

using namespace std;

void vector_traditional_forloop_test() {
    int n = 1;
    for (int i = 0; i < NUMBER_OF_LOOPS; i++, n *= 2) {
        auto start = chrono::steady_clock::now();
        float sum = vector_traditional_forloop(n);
        auto end = chrono::steady_clock::now();
        long time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    }
}

void vector_ranged_forloop_test() {
    int n = 1;
    for (int i = 0; i < NUMBER_OF_LOOPS; i++, n *= 2) {
        auto start = chrono::steady_clock::now();
        float sum = vector_ranged_forloop(n);
        auto end = chrono::steady_clock::now();
        long time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    }
}

void vector_preallocated_ranged_forloop_test() {
    int n = 1;
    for (int i = 0; i < NUMBER_OF_LOOPS; i++, n *= 2) {
        auto start = chrono::steady_clock::now();
        float sum = vector_preallocated_ranged_forloop(n);
        auto end = chrono::steady_clock::now();
        long time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    }
}

void array_stack_test() {
    int n = 1;
    for (int i = 0; i < NUMBER_OF_LOOPS; i++, n *= 2) {
        auto start = chrono::steady_clock::now();
        float sum = array_stack_allocation(n);
        auto end = chrono::steady_clock::now();
        long time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    }
}

void array_heap_test() {
    int n = 1;
    for (int i = 0; i < NUMBER_OF_LOOPS; i++, n *= 2) {
        auto start = chrono::steady_clock::now();
        float sum = array_heap_allocation(n);
        auto end = chrono::steady_clock::now();
        long time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
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
    vector_preallocated_ranged_forloop_test();
    end = chrono::steady_clock::now();
    time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "vector preallocated ranged forloop total time: " << time_us << "\n\n";

    start = chrono::steady_clock::now();
    array_stack_test();
    end = chrono::steady_clock::now();
    time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "array stack allocation total time: " << time_us << "\n\n";

    start = chrono::steady_clock::now();
    array_heap_test();
    end = chrono::steady_clock::now();
    time_us = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "array heap allocation total time: " << time_us << "\n\n";
}

