#include "Iteration.hpp"

// CSWM model;
CSWM model;
int main(void) {
    clock_t start, stop;
    start = clock();
    Init::Init2d(model);
    Iteration::Leapfrog(model);
    stop = clock();
    std::cout << double(stop - start) / CLOCKS_PER_SEC << std::endl;
    return 0;
}