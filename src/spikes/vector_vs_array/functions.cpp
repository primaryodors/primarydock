#include "functions.hpp"
#include <iostream>
#include <vector>

// Made up problem: Calculate sum of integers 1..n
// in a really stupid way: create a vector/array
// with the integers, then loop through them accumulating
// into a float variable!

float vector_summer(int n) {
    std::vector<int> integers;

    for (int i = 1; i <= n; i++)
        integers.push_back(i);

    float sum = 0;
    for (int i = 1; i <= n; i++)
        sum += (float)integers[i];

    return sum;
}

float array_summer(int n) {
    int* integers = new int[n];

    for (int i = 0; i < n; i++)
        integers[i] = i;

    float sum = 0;
    for (int i = 0; i < n; i++)
        sum += (float)integers[i];

	delete integers;
    return sum;
}
