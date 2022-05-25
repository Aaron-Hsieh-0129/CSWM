#include "Init.hpp"

void Init::Init2d(CSWM & model) {
    // Init h, u, v
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                model.cswm[p].h[i][j] = GetH(model.cswm[p].lon[i][j], model.cswm[p].lat[i][j]);
                model.cswm[p].u[i][j] = GetU(0.5*(model.cswm[p].lon[i-1][j]+model.cswm[p].lon[i][j]), model.cswm[p].lat[i][j]);
                model.cswm[p].v[i][j] = GetV(model.cswm[p].lon[i][j]);

                model.cswm[p].hm[i][j] = model.cswm[p].h[i][j]; 
                model.cswm[p].um[i][j] = model.cswm[p].u[i][j]; 
                model.cswm[p].vm[i][j] = model.cswm[p].v[i][j];   

                model.cswm[p].hp[i][j] = model.cswm[p].h[i][j]; 
                model.cswm[p].up[i][j] = model.cswm[p].u[i][j]; 
                model.cswm[p].vp[i][j] = model.cswm[p].v[i][j];             
            }
        }
    }
}


double Init::GetH(double lon, double lat) {
    double h0 = 1000;
    double lonC = 3. * M_PI / 2., latC = 0;
    double rd = radius * acos(sin(latC) * sin(lat) + cos(latC) * cos(lat) * cos(lon));
    double r0 = radius / 3.;
    if (rd < r0) return h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return 0;
}

double Init::GetU(double lon, double lat) {
    double u0 = 2 * M_PI * radius / (12. * 86400);
    double alpha0 = 0.;
    double u = u0 * (cos(alpha0) * cos(lat) + sin(alpha0) * cos(lon) * sin(lat));
    return u;
}

double Init::GetV(double lon) {
    double alpha0 = 0.;
    double u0 = 2 * M_PI * radius / (6. * 86400);
    double v = u0 * sin(alpha0) * cos(lon);
    return v;
}