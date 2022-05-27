#include "Iteration.hpp"

// CSWM model;
CSWM model;
int main(void) {
    Init::Init2d(model);
    Iteration::Leapfrog(model);
    // for (int idx = 0; idx < NY; idx++) {
    //     double a = model.cswm[5].lon[1][idx];
    //     double b = model.cswm[3].lon[idx][0];
    //     if (a < 0) a += 2 * M_PI;
    //     if (b < 0) b += 2 * M_PI;
    //     std::cout << a << " " << b << " ";
    //     std::cout << std::endl;
    // }
    return 0;
}