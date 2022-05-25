#include "Declare.hpp"


CSWM::patch::patch() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            hp[i][j] = h[i][j] = hm[i][j] = 0.;
            up[i][j] = u[i][j] = um[i][j] = 0.;
            vp[i][j] = v[i][j] = vm[i][j] = 0.;
            lon[i][j] = lat[i][j] = 0.;
            x[i][j] = y[i][j] = 0.;
        }
    }
}

CSWM::CSWM() {
    double alpha[NX], beta[NY];
    double alpha2D[NX][NY], beta2D[NX][NY];
    // TODO: deal with the odd NX/NY and even NX/NY
    for (int i = 0; i < NX; i++) alpha[i] = -M_PI/4 + (i-0.5) * ((M_PI / 2.) / (NX-2));
    for (int j = 0; j < NY; j++) beta[j] = -M_PI/4 + (j-0.5) * ((M_PI / 2.) / (NY-2));

    // init alpha/beta
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                alpha2D[i][j] = alpha[i];
                beta2D[i][j] = beta[j];
            }
        }
    }


    // calculate gamma, sqrtG, gUpper
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            gamma[i][j] = sqrt(1 + pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D[i][j]), 2));
            sqrtG[i][j] = 1. / (pow(gamma[i][j], 3) * pow(cos(alpha2D[i][j]), 2) * pow(cos(beta2D[i][j]), 2));
            gUpper[i][j][0] = pow(gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2) * (1 + pow(tan(beta2D[i][j]), 2));
            gUpper[i][j][1] = pow(gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2) * (tan(alpha2D[i][j]) * tan(beta2D[i][j]));
            gUpper[i][j][2] = gUpper[i][j][1];
            gUpper[i][j][3] = pow(gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2) * (1 + pow(tan(alpha2D[i][j]), 2));
        }
    }


     // calculate lon/lat, x/y, A
    for (int p = 0; p < 6; p++) {
        if (p == 0 || p == 1 || p == 2 || p == 3) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    // lon/lat
                    cswm[p].lon[i][j] = alpha2D[i][j] + p * M_PI / 2.;
                    cswm[p].lat[i][j] = atan(tan(beta2D[i][j]) * cos(alpha2D[i][j]));

                    // x/y
                    cswm[p].x[i][j] = radius * (cswm[p].lon[i][j] - p * M_PI / 2.);
                    cswm[p].y[i][j] = radius * atan(tan(cswm[p].lat[i][j]) / cos(cswm[p].lon[i][j] - p * M_PI / 2.));

                    // A/AInverse
                    cswm[p].A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * gamma[i][j] * cos(beta2D[i][j]);
                    cswm[p].A[i][j][1] = 0;
                    cswm[p].A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-tan(alpha2D[i][j]) * sin(beta2D[i][j]));
                    cswm[p].A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) / cos(beta2D[i][j]);
                    
                    cswm[p].AInverse[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) / cos(beta2D[i][j]);
                    cswm[p].AInverse[i][j][1] = 0;
                    cswm[p].AInverse[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (tan(alpha2D[i][j]) * sin(beta2D[i][j]));
                    cswm[p].AInverse[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]));
                }
            }
        }

        if (p == 4) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    // lon/lat
                    cswm[p].lon[i][j] = atan2(tan(alpha2D[i][j]), -tan(beta2D[i][j]));
                    cswm[p].lat[i][j] = atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D[i][j]), 2)));

                    // x/y
                    cswm[p].x[i][j] = radius * atan(sin(cswm[p].lon[i][j]) / tan(cswm[p].lat[i][j]));
                    cswm[p].y[i][j] = radius * atan(-cos(cswm[p].lon[i][j]) / tan(cswm[p].lat[i][j]));

                    // A/AInverse
                    cswm[p].A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][1] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon[i][j]));
                    
                    cswm[p].AInverse[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon[i][j]));
                }
            }
        }

        if (p == 5) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    // lon/lat
                    cswm[p].lon[i][j] = atan2(tan(alpha2D[i][j]), tan(beta2D[i][j]));
                    cswm[p].lat[i][j] = -atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D[i][j]), 2)));

                    // x/y
                    cswm[p].x[i][j] = radius * atan(-sin(cswm[p].lon[i][j]) / tan(cswm[p].lat[i][j]));
                    cswm[p].y[i][j] = radius * atan(-cos(cswm[p].lon[i][j]) / tan(cswm[p].lat[i][j]));

                    // A/AInverse
                    cswm[p].A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][1] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon[i][j]));

                    cswm[p].AInverse[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon[i][j]));
                }
            }
        }
    }   
}

// TODO: Interpolation
void CSWM::BoundaryProcess(CSWM &model) {
    // patch 1
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[0].hp[0][idx] = model.cswm[3].hp[NX-2][idx];
        model.cswm[0].up[0][idx] = model.cswm[3].up[NX-2][idx];
        model.cswm[0].vp[0][idx] = model.cswm[3].vp[NX-2][idx];

        // right
        model.cswm[0].hp[NX-1][idx] = model.cswm[1].hp[1][idx];
        model.cswm[0].up[NX-1][idx] = model.cswm[1].up[1][idx];
        model.cswm[0].vp[NX-1][idx] = model.cswm[1].vp[1][idx];

        // up
        model.cswm[0].hp[idx][NY-1] = model.cswm[4].hp[idx][1];
        model.cswm[0].up[idx][NY-1] = model.cswm[4].up[idx][1];
        model.cswm[0].vp[idx][NY-1] = model.cswm[4].vp[idx][1];

        // bottom
        model.cswm[0].hp[idx][0] = model.cswm[5].hp[idx][NY-2];
        model.cswm[0].up[idx][0] = model.cswm[5].up[idx][NY-2];
        model.cswm[0].vp[idx][0] = model.cswm[5].vp[idx][NY-2];
    }

    // patch 2
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[1].hp[0][idx] = model.cswm[0].hp[NX-2][idx];
        model.cswm[1].up[0][idx] = model.cswm[0].up[NX-2][idx];
        model.cswm[1].vp[0][idx] = model.cswm[0].vp[NX-2][idx];

        // right
        model.cswm[1].hp[NX-1][idx] = model.cswm[2].hp[1][idx];
        model.cswm[1].up[NX-1][idx] = model.cswm[2].up[1][idx];
        model.cswm[1].vp[NX-1][idx] = model.cswm[2].vp[1][idx];

        // up
        model.cswm[1].hp[idx][NY-1] = model.cswm[4].hp[NX-2][idx];
        model.cswm[1].up[idx][NY-1] = model.cswm[4].vp[NX-2][idx];
        model.cswm[1].vp[idx][NY-1] = -model.cswm[4].up[NX-2][idx];

        // bottom
        model.cswm[1].hp[idx][0] = model.cswm[5].hp[NX-2][(NY-1)-idx];
        model.cswm[1].up[idx][0] = -model.cswm[5].vp[NX-2][(NY-1)-idx];
        model.cswm[1].vp[idx][0] = model.cswm[5].up[NX-2][(NY-1)-idx];
    }

    // patch 3
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[2].hp[0][idx] = model.cswm[1].hp[NX-2][idx];
        model.cswm[2].up[0][idx] = model.cswm[1].up[NX-2][idx];
        model.cswm[2].vp[0][idx] = model.cswm[1].vp[NX-2][idx];

        // right
        model.cswm[2].hp[NX-1][idx] = model.cswm[3].hp[1][idx];
        model.cswm[2].up[NX-1][idx] = model.cswm[3].up[1][idx];
        model.cswm[2].vp[NX-1][idx] = model.cswm[3].vp[1][idx];

        // up
        model.cswm[2].hp[idx][NY-1] = model.cswm[4].hp[(NX-1)-idx][NY-2];
        model.cswm[2].up[idx][NY-1] = -model.cswm[4].up[(NX-1)-idx][NY-2];
        model.cswm[2].vp[idx][NY-1] = -model.cswm[4].vp[(NX-1)-idx][NY-2];

        // bottom
        model.cswm[2].hp[idx][0] = model.cswm[5].hp[(NX-1)-idx][1];
        model.cswm[2].up[idx][0] = -model.cswm[5].up[(NX-1)-idx][1];
        model.cswm[2].vp[idx][0] = -model.cswm[5].vp[(NX-1)-idx][1];
    }

    // patch 4
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[3].hp[0][idx] = model.cswm[2].hp[NX-2][idx];
        model.cswm[3].up[0][idx] = model.cswm[2].up[NX-2][idx];
        model.cswm[3].vp[0][idx] = model.cswm[2].vp[NX-2][idx];

        // right
        model.cswm[3].hp[NX-1][idx] = model.cswm[0].hp[1][idx];
        model.cswm[3].up[NX-1][idx] = model.cswm[0].up[1][idx];
        model.cswm[3].vp[NX-1][idx] = model.cswm[0].vp[1][idx];

        // up
        model.cswm[3].hp[idx][NY-1] = model.cswm[4].hp[1][(NY-1)-idx];
        model.cswm[3].up[idx][NY-1] = -model.cswm[4].vp[1][(NY-1)-idx];
        model.cswm[3].vp[idx][NY-1] = model.cswm[4].up[1][(NY-1)-idx];

        // bottom
        model.cswm[3].hp[idx][0] = model.cswm[5].hp[1][idx];
        model.cswm[3].up[idx][0] = model.cswm[5].vp[1][idx];
        model.cswm[3].vp[idx][0] = -model.cswm[5].up[1][idx];
    }

    // patch 5
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[4].hp[0][idx] = model.cswm[3].hp[(NX-1)-idx][NY-2];
        model.cswm[4].up[0][idx] = model.cswm[3].vp[(NX-1)-idx][NY-2];
        model.cswm[4].vp[0][idx] = -model.cswm[3].up[(NX-1)-idx][NY-2];

        // right
        model.cswm[4].hp[NX-1][idx] = model.cswm[1].hp[idx][NY-2];
        model.cswm[4].up[NX-1][idx] = -model.cswm[1].vp[idx][NY-2];
        model.cswm[4].vp[NX-1][idx] = model.cswm[1].up[idx][NY-2];

        // up
        model.cswm[4].hp[idx][NY-1] = model.cswm[2].hp[(NX-1)-idx][NY-2];
        model.cswm[4].up[idx][NY-1] = -model.cswm[2].up[(NX-1)-idx][NY-2];
        model.cswm[4].vp[idx][NY-1] = -model.cswm[2].vp[(NX-1)-idx][NY-2];

        // bottom
        model.cswm[4].hp[idx][0] = model.cswm[0].hp[idx][1];
        model.cswm[4].up[idx][0] = model.cswm[0].up[idx][1];
        model.cswm[4].vp[idx][0] = model.cswm[0].vp[idx][1];
    }

    // patch 6
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[5].hp[0][idx] = model.cswm[3].hp[idx][1];
        model.cswm[5].up[0][idx] = -model.cswm[3].vp[idx][1];
        model.cswm[5].vp[0][idx] = model.cswm[3].up[idx][1];

        // right
        model.cswm[5].hp[NX-1][idx] = model.cswm[1].hp[(NX-1)-idx][1];
        model.cswm[5].up[NX-1][idx] = model.cswm[1].vp[(NX-1)-idx][1];
        model.cswm[5].vp[NX-1][idx] = -model.cswm[1].up[(NX-1)-idx][1];

        // up
        model.cswm[5].hp[idx][NY-1] = model.cswm[0].hp[idx][1];
        model.cswm[5].up[idx][NY-1] = model.cswm[0].up[idx][1];
        model.cswm[5].vp[idx][NY-1] = model.cswm[0].vp[idx][1];

        // bottom
        model.cswm[5].hp[idx][0] = model.cswm[2].hp[(NX-1)-idx][1];
        model.cswm[5].up[idx][0] = -model.cswm[2].up[(NX-1)-idx][1];
        model.cswm[5].vp[idx][0] = -model.cswm[2].vp[(NX-1)-idx][1];
    }
}
