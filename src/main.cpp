#include "Iteration.hpp"

// CSWM model;
CSWM model;
int main(void) {
    Init::Init2d(model);
    Iteration::Leapfrog(model);
    return 0;
}