#include <chrono>
#include <iostream>
#include <unistd.h>
#include "functions.hpp"

using namespace std;

int main() {
    auto start = chrono::steady_clock::now();

    sleep(1);

    auto end = chrono::steady_clock::now();

    cout << "Elapsed time in microseconds: "
         << chrono::duration_cast<chrono::microseconds>(end - start).count()
         << " µs" << endl;

    hello();
}