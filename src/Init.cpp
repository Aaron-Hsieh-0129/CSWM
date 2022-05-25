#include "Init.hpp"

void Init::Init2d(CSWM & model) {
    // Init h
    // for (int i = 0; i < (NX-2)*4; i++) {
    //     for (int j = 0; j < (NY-2)*3; j++) {
    //         model.cswm[0].h[i][j] = 
    //     }
    // }
}


double Init::GetH(int phyIdxLon, int phyIdxLat) {
    double h0 = 1000;
    double lonC = 3. * M_PI / 2., latC = 0;
    double rd = radius * acos(sin(lonC) * sin((phyIdxLon * DX) * (180. / M_PI)) + cos(latC) * sin((phyIdxLat * DY) * (180. / M_PI)) * cos(DX * (180. / M_PI)));
    double r0 = radius / 3.;
    if (rd < r0) return h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return 0;
}

double Init::GetU(int phyIdxLon, int phyIdxLat) {
    double u0 = 2 * M_PI * radius / (6. * 86400);
    double alpha0 = 0.;
    double u = u0 * (cos(alpha0) * cos(phyIdxLat * DY * (180. / M_PI)) + sin(alpha0) * sin(phyIdxLat * DX * (180. / M_PI)));
    return u;
}

double Init::GetV(int phyIdxLon) {
    double alpha0 = 0.;
    double u0 = 2 * M_PI * radius / (6. * 86400);
    double v = u0 * sin(alpha0) * cos(phyIdxLon * DX * (180. / M_PI));
    return v;
}