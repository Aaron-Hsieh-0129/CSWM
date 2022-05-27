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

    /*
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[0].hp[0][idx] = model.cswm[3].hp[NX-2][idx];
        // right
        model.cswm[0].hp[NX-1][idx] = model.cswm[1].hp[1][idx];
        // up
        model.cswm[0].hp[idx][NY-1] = model.cswm[4].hp[idx][1];
        // bottom
        model.cswm[0].hp[idx][0] = model.cswm[5].hp[idx][NY-2];
    }

    // patch 2
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[1].hp[0][idx] = model.cswm[0].hp[NX-2][idx];
        // right
        model.cswm[1].hp[NX-1][idx] = model.cswm[2].hp[1][idx];
        // up
        model.cswm[1].hp[idx][NY-1] = model.cswm[4].hp[NX-2][idx];
        // bottom
        model.cswm[1].hp[idx][0] = model.cswm[5].hp[NX-2][(NY-1)-idx];
    }

    // patch 3
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[2].hp[0][idx] = model.cswm[1].hp[NX-2][idx];
        // right
        model.cswm[2].hp[NX-1][idx] = model.cswm[3].hp[1][idx];
        // up
        model.cswm[2].hp[idx][NY-1] = model.cswm[4].hp[(NX-1)-idx][NY-2];
        // bottom
        model.cswm[2].hp[idx][0] = model.cswm[5].hp[(NX-1)-idx][1];
    }

    // patch 4
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[3].hp[0][idx] = model.cswm[2].hp[NX-2][idx];
        // right
        model.cswm[3].hp[NX-1][idx] = model.cswm[0].hp[1][idx];
        // up
        model.cswm[3].hp[idx][NY-1] = model.cswm[4].hp[1][(NY-1)-idx];
        // bottom
        model.cswm[3].hp[idx][0] = model.cswm[5].hp[1][idx];
    }

    // patch 5
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[4].hp[0][idx] = model.cswm[3].hp[(NX-1)-idx][NY-2];
        // right
        model.cswm[4].hp[NX-1][idx] = model.cswm[1].hp[idx][NY-2];
        // up
        model.cswm[4].hp[idx][NY-1] = model.cswm[2].hp[(NX-1)-idx][NY-2];
        // bottom
        model.cswm[4].hp[idx][0] = model.cswm[0].hp[idx][NY-2];
    }

    // patch 6
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[5].hp[0][idx] = model.cswm[3].hp[idx][1];
        // right
        model.cswm[5].hp[NX-1][idx] = model.cswm[1].hp[(NX-1)-idx][1];
        // up
        model.cswm[5].hp[idx][NY-1] = model.cswm[0].hp[idx][1];
        // bottom
        model.cswm[5].hp[idx][0] = model.cswm[2].hp[(NX-1)-idx][1];
    }

    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                model.cswm[p].u[i][j] = model.cswm[p].AInverse[i][j][0] * model.cswm[p].u[i][j] + model.cswm[p].AInverse[i][j][1] * model.cswm[p].v[i][j];
                model.cswm[p].v[i][j] = model.cswm[p].AInverse[i][j][2] * model.cswm[p].u[i][j] + model.cswm[p].AInverse[i][j][3] * model.cswm[p].v[i][j];
                model.cswm[p].um[i][j] = model.cswm[p].u[i][j];
                model.cswm[p].vm[i][j] = model.cswm[p].v[i][j];
                model.cswm[p].up[i][j] = model.cswm[p].u[i][j];
                model.cswm[p].vp[i][j] = model.cswm[p].v[i][j];
            }
        }
    }
    */



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
    double alpha0 = M_PI / 2.;
    double u = u0 * (cos(alpha0) * cos(lat) + sin(alpha0) * sin(lon) * sin(lat));
    return u;
}

double Init::GetV(double lon) {
    double alpha0 = M_PI / 2.;
    double u0 = 2 * M_PI * radius / (12. * 86400);
    double v = u0 * sin(alpha0) * cos(lon);
    return v;
}