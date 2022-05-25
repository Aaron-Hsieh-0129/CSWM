#include "Declare.hpp"

// CSWM model;
CSWM model;
int main(void) {
    model.BoundaryProcess(model);

    // for (int i = 1; i < NX-1; i++) {
    //     std::cout << model.cswm[5].x[i][1] << " ";
    // }
    // std::cout << std::endl;
    // std::cout << std::endl;

    for (int j = 1; j < NY-1; j++) {
        std::cout << model.cswm[0].lat[1][j] << " ";
    }
    // std::cout << std::endl;
    // for (int i = 0; i < NX; i++) {
    //     std::cout << model.gUpper[i][1][0] << " " << model.gUpper[i][1][1] << " " << model.gUpper[i][1][2] << " " << model.gUpper[i][1][3];
    //     std::cout << std::endl;
    // }
    return 0;
}