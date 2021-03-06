#include "Declare.hpp"


CSWM::patch::patch() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            hp[i][j] = h[i][j] = hm[i][j] = FILLVALUE;
            up[i][j] = u[i][j] = um[i][j] = FILLVALUE;
            vp[i][j] = v[i][j] = vm[i][j] = FILLVALUE;
            lon[i][j] = lat[i][j] = FILLVALUE;
            lon_u[i][j] = lat_v[i][j] = FILLVALUE;
            x[i][j] = y[i][j] = FILLVALUE;
            x_u[i][j] = y_u[i][j] = FILLVALUE;
            x_v[i][j] = y_v[i][j] = FILLVALUE;
        }
    }
}

CSWM::CSWM() {
    double alpha[NX], beta[NY], alpha_u[NX], beta_v[NY];
    for (int i = 0; i < NX; i++) {
        alpha[i] = -M_PI/4 + (i-0.5) * ((M_PI / 2.) / (NX-2));
        alpha_u[i] = -M_PI/4 + (i-1) * ((M_PI / 2.) / (NX-2));
    }
    for (int j = 0; j < NY; j++) {
        beta[j] = -M_PI/4 + (j-0.5) * ((M_PI / 2.) / (NY-2));
        beta_v[j] = -M_PI/4 + (j-1) * ((M_PI / 2.) / (NY-2));
    }

    // init alpha/beta
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            alpha2D[i][j] = alpha[i];
            beta2D[i][j] = beta[j];
            alpha2D_u[i][j] = alpha_u[i];
            beta2D_v[i][j] = beta_v[j];
        }
    }
    // bug test
    // for (int idx = 0; idx < NX; idx++) {
    //     std::cout << alpha2D[idx][NX/2] << " ";
    // }
    // std::cout << std::endl;
    // for (int idx = 0; idx < NX; idx++) {
    //     std::cout << beta2D[NX/2][idx] << " ";
    // }
    // std::cout << std::endl;

    // calculate gamma, sqrtG, gUpper
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            gamma[i][j] = sqrt(1 + pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D[i][j]), 2));
            sqrtG[i][j] = 1. / (pow(gamma[i][j], 3) * pow(cos(alpha2D[i][j]), 2) * pow(cos(beta2D[i][j]), 2));
            gUpper[i][j][0] = pow(gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2) * (1 + pow(tan(beta2D[i][j]), 2));
            gUpper[i][j][1] = pow(gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2) * (tan(alpha2D[i][j]) * tan(beta2D[i][j]));
            gUpper[i][j][2] = gUpper[i][j][1];
            gUpper[i][j][3] = pow(gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2) * (1 + pow(tan(alpha2D[i][j]), 2));

            gLower[i][j][0] = 1. / (pow(gamma[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2)) * (1 + pow(tan(alpha2D[i][j]), 2));
            gLower[i][j][1] = 1. / (pow(gamma[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2)) * (-tan(alpha2D[i][j]) * tan(beta2D[i][j]));
            gLower[i][j][2] = gLower[i][j][1];
            gLower[i][j][3] = 1. / (pow(gamma[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2)) * (1 + pow(tan(beta2D[i][j]), 2));
        }
    }

    // calculate gamma_u, sqrtG_u, gUpper_u
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            gamma_u[i][j] = sqrt(1 + pow(tan(alpha2D_u[i][j]), 2) + pow(tan(beta2D[i][j]), 2));
            sqrtG_u[i][j] = 1. / (pow(gamma_u[i][j], 3) * pow(cos(alpha2D_u[i][j]), 2) * pow(cos(beta2D[i][j]), 2));
            gUpper_u[i][j][0] = pow(gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]), 2) * (1 + pow(tan(beta2D[i][j]), 2));
            gUpper_u[i][j][1] = pow(gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]), 2) * (tan(alpha2D_u[i][j]) * tan(beta2D[i][j]));
            gUpper_u[i][j][2] = gUpper_u[i][j][1];
            gUpper_u[i][j][3] = pow(gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]), 2) * (1 + pow(tan(alpha2D_u[i][j]), 2));

            gLower_u[i][j][0] = 1. / (pow(gamma_u[i][j], 4) * pow(cos(alpha2D_u[i][j]) * cos(beta2D[i][j]), 2)) * (1 + pow(tan(alpha2D_u[i][j]), 2));
            gLower_u[i][j][1] = 1. / (pow(gamma_u[i][j], 4) * pow(cos(alpha2D_u[i][j]) * cos(beta2D[i][j]), 2)) * (-tan(alpha2D_u[i][j]) * tan(beta2D[i][j]));
            gLower_u[i][j][2] = gLower_u[i][j][1];
            gLower_u[i][j][3] = 1. / (pow(gamma_u[i][j], 4) * pow(cos(alpha2D_u[i][j]) * cos(beta2D[i][j]), 2)) * (1 + pow(tan(beta2D[i][j]), 2));
        }
    }

    // calculate gamma_v, sqrtG_v, gUpper_v
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            gamma_v[i][j] = sqrt(1 + pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D_v[i][j]), 2));
            sqrtG_v[i][j] = 1. / (pow(gamma_v[i][j], 3) * pow(cos(alpha2D[i][j]), 2) * pow(cos(beta2D_v[i][j]), 2));
            gUpper_v[i][j][0] = pow(gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]), 2) * (1 + pow(tan(beta2D_v[i][j]), 2));
            gUpper_v[i][j][1] = pow(gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]), 2) * (tan(alpha2D[i][j]) * tan(beta2D_v[i][j]));
            gUpper_v[i][j][2] = gUpper_v[i][j][1];
            gUpper_v[i][j][3] = pow(gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]), 2) * (1 + pow(tan(alpha2D[i][j]), 2));

            gLower_v[i][j][0] = 1. / (pow(gamma_v[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D_v[i][j]), 2)) * (1 + pow(tan(alpha2D[i][j]), 2));
            gLower_v[i][j][1] = 1. / (pow(gamma_v[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D_v[i][j]), 2)) * (-tan(alpha2D[i][j]) * tan(beta2D_v[i][j]));
            gLower_v[i][j][2] = gLower_v[i][j][1];
            gLower_v[i][j][3] = 1. / (pow(gamma_v[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D_v[i][j]), 2)) * (1 + pow(tan(beta2D_v[i][j]), 2));
        }
    }

     // calculate lon/lat, x/y, A/AInverse
    for (int p = 0; p < 6; p++) {
        if (p == 0 || p == 1 || p == 2 || p == 3) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    // lon/lat
                    cswm[p].lon[i][j] = alpha2D[i][j] + p * M_PI / 2.;
                    cswm[p].lat[i][j] = atan(tan(beta2D[i][j]) * cos(alpha2D[i][j]));

                    // lon_u/lat_u
                    cswm[p].lon_u[i][j] = alpha2D_u[i][j] + p * M_PI / 2.;
                    cswm[p].lat_u[i][j] = atan(tan(beta2D[i][j]) * cos(alpha2D_u[i][j]));

                    // lon_v/lat_v
                    cswm[p].lon_v[i][j] = alpha2D[i][j] + p * M_PI / 2.;
                    cswm[p].lat_v[i][j] = atan(tan(beta2D_v[i][j]) * cos(alpha2D[i][j]));

                    // x/y
                    cswm[p].x[i][j] = radius * (cswm[p].lon[i][j] - p * M_PI / 2.);
                    cswm[p].y[i][j] = radius * atan(tan(cswm[p].lat[i][j]) / cos(cswm[p].lon[i][j] - p * M_PI / 2.));

                    // x_u/y_u
                    cswm[p].x_u[i][j] = radius * (cswm[p].lon_u[i][j] - p * M_PI / 2.);
                    cswm[p].y_u[i][j] = radius * atan(tan(cswm[p].lat_u[i][j]) / cos(cswm[p].lon_u[i][j] - p * M_PI / 2.));

                    // x_v/y_v
                    cswm[p].x_v[i][j] = radius * (cswm[p].lon_v[i][j] - p * M_PI / 2.);
                    cswm[p].y_v[i][j] = radius * atan(tan(cswm[p].lat_v[i][j]) / cos(cswm[p].lon_v[i][j] - p * M_PI / 2.));

                    // A/AInverse
                    cswm[p].A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * gamma[i][j] * cos(beta2D[i][j]);
                    cswm[p].A[i][j][1] = 0;
                    cswm[p].A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-tan(alpha2D[i][j]) * sin(beta2D[i][j]));
                    cswm[p].A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) / cos(beta2D[i][j]);
                    
                    cswm[p].AInverse[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) / cos(beta2D[i][j]);
                    cswm[p].AInverse[i][j][1] = 0;
                    cswm[p].AInverse[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (tan(alpha2D[i][j]) * sin(beta2D[i][j]));
                    cswm[p].AInverse[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]));

                    // A_u/AInverse_u
                    cswm[p].A_u[i][j][0] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * gamma_u[i][j] * cos(beta2D[i][j]);
                    cswm[p].A_u[i][j][1] = 0;
                    cswm[p].A_u[i][j][2] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * (-tan(alpha2D_u[i][j]) * sin(beta2D[i][j]));
                    cswm[p].A_u[i][j][3] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) / cos(beta2D[i][j]);
                    
                    cswm[p].AInverse_u[i][j][0] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) / cos(beta2D[i][j]);
                    cswm[p].AInverse_u[i][j][1] = 0;
                    cswm[p].AInverse_u[i][j][2] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (tan(alpha2D_u[i][j]) * sin(beta2D[i][j]));
                    cswm[p].AInverse_u[i][j][3] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (gamma_u[i][j] * cos(beta2D[i][j]));

                    // A_v/AInverse_v
                    cswm[p].A_v[i][j][0] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * gamma_v[i][j] * cos(beta2D_v[i][j]);
                    cswm[p].A_v[i][j][1] = 0;
                    cswm[p].A_v[i][j][2] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * (-tan(alpha2D[i][j]) * sin(beta2D_v[i][j]));
                    cswm[p].A_v[i][j][3] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) / cos(beta2D_v[i][j]);
                    
                    cswm[p].AInverse_v[i][j][0] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) / cos(beta2D_v[i][j]);
                    cswm[p].AInverse_v[i][j][1] = 0;
                    cswm[p].AInverse_v[i][j][2] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (tan(alpha2D[i][j]) * sin(beta2D_v[i][j]));
                    cswm[p].AInverse_v[i][j][3] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (gamma_v[i][j] * cos(beta2D_v[i][j]));

                }
            }
        }

        if (p == 4) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    // lon/lat
                    cswm[p].lon[i][j] = atan2(tan(alpha2D[i][j]), -tan(beta2D[i][j]));
                    cswm[p].lat[i][j] = atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D[i][j]), 2)));

                    // lon_u/lat_u
                    cswm[p].lon_u[i][j] = atan2(tan(alpha2D_u[i][j]), -tan(beta2D[i][j]));
                    cswm[p].lat_u[i][j] = atan(1 / sqrt(pow(tan(alpha2D_u[i][j]), 2) + pow(tan(beta2D[i][j]), 2)));

                    // lon_v/lat_v
                    cswm[p].lon_v[i][j] = atan2(tan(alpha2D[i][j]), -tan(beta2D_v[i][j]));
                    cswm[p].lat_v[i][j] = atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D_v[i][j]), 2)));

                    // x/y
                    cswm[p].x[i][j] = radius * atan(sin(cswm[p].lon[i][j]) / tan(cswm[p].lat[i][j]));
                    cswm[p].y[i][j] = radius * atan(-cos(cswm[p].lon[i][j]) / tan(cswm[p].lat[i][j]));

                    // x_u/y_u
                    cswm[p].x_u[i][j] = radius * atan(sin(cswm[p].lon_u[i][j]) / tan(cswm[p].lat_u[i][j]));
                    cswm[p].y_u[i][j] = radius * atan(-cos(cswm[p].lon_u[i][j]) / tan(cswm[p].lat_u[i][j]));

                    // x_v/y_v
                    cswm[p].x_v[i][j] = radius * atan(sin(cswm[p].lon_v[i][j]) / tan(cswm[p].lat_v[i][j]));
                    cswm[p].y_v[i][j] = radius * atan(-cos(cswm[p].lon_v[i][j]) / tan(cswm[p].lat_v[i][j]));

                    // A/AInverse
                    cswm[p].A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][1] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon[i][j]));
                    
                    cswm[p].AInverse[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon[i][j]));

                    // A_u/AInverse_u
                    cswm[p].A_u[i][j][0] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * (gamma_u[i][j] * cos(beta2D[i][j]) / cos(alpha2D_u[i][j]) * cos(cswm[p].lon_u[i][j]));
                    cswm[p].A_u[i][j][1] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * (gamma_u[i][j] * cos(alpha2D_u[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon_u[i][j]));
                    cswm[p].A_u[i][j][2] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * (-cos(beta2D[i][j]) / cos(alpha2D_u[i][j]) * sin(cswm[p].lon_u[i][j]));
                    cswm[p].A_u[i][j][3] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D_u[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon_u[i][j]));
                    
                    cswm[p].AInverse_u[i][j][0] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D_u[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon_u[i][j]));
                    cswm[p].AInverse_u[i][j][1] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (-gamma_u[i][j] * cos(alpha2D_u[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon_u[i][j]));
                    cswm[p].AInverse_u[i][j][2] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (cos(beta2D[i][j]) / cos(alpha2D_u[i][j]) * sin(cswm[p].lon_u[i][j]));
                    cswm[p].AInverse_u[i][j][3] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D_u[i][j]) * cos(cswm[p].lon_u[i][j]));

                    // A_v/AInverse_v
                    cswm[p].A_v[i][j][0] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * (gamma_v[i][j] * cos(beta2D_v[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon_v[i][j]));
                    cswm[p].A_v[i][j][1] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * (gamma_v[i][j] * cos(alpha2D[i][j]) / cos(beta2D_v[i][j]) * sin(cswm[p].lon_v[i][j]));
                    cswm[p].A_v[i][j][2] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * (-cos(beta2D_v[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon_v[i][j]));
                    cswm[p].A_v[i][j][3] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D_v[i][j]) * cos(cswm[p].lon_v[i][j]));
                    
                    cswm[p].AInverse_v[i][j][0] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D_v[i][j]) * cos(cswm[p].lon_v[i][j]));
                    cswm[p].AInverse_v[i][j][1] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (-gamma_v[i][j] * cos(alpha2D[i][j]) / cos(beta2D_v[i][j]) * sin(cswm[p].lon_v[i][j]));
                    cswm[p].AInverse_v[i][j][2] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (cos(beta2D_v[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon_v[i][j]));
                    cswm[p].AInverse_v[i][j][3] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (gamma_v[i][j] * cos(beta2D_v[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon_v[i][j]));
                }
            }
        }

        if (p == 5) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    // lon/lat
                    cswm[p].lon[i][j] = atan2(tan(alpha2D[i][j]), tan(beta2D[i][j]));
                    cswm[p].lat[i][j] = -atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D[i][j]), 2)));

                    // lon_u/lat_u
                    cswm[p].lon_u[i][j] = atan2(tan(alpha2D_u[i][j]), tan(beta2D[i][j]));
                    cswm[p].lat_u[i][j] = -atan(1 / sqrt(pow(tan(alpha2D_u[i][j]), 2) + pow(tan(beta2D[i][j]), 2)));

                    // lon_v/lat_v
                    cswm[p].lon_v[i][j] = atan2(tan(alpha2D[i][j]), tan(beta2D_v[i][j]));
                    cswm[p].lat_v[i][j] = -atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D_v[i][j]), 2)));

                    // x/y
                    cswm[p].x[i][j] = radius * atan(-sin(cswm[p].lon[i][j]) / tan(cswm[p].lat[i][j]));
                    cswm[p].y[i][j] = radius * atan(-cos(cswm[p].lon[i][j]) / tan(cswm[p].lat[i][j]));

                    // x_u/y_u
                    cswm[p].x_u[i][j] = radius * atan(-sin(cswm[p].lon_u[i][j]) / tan(cswm[p].lat_u[i][j]));
                    cswm[p].y_u[i][j] = radius * atan(-cos(cswm[p].lon_u[i][j]) / tan(cswm[p].lat_u[i][j]));

                    // x_v/y_v
                    cswm[p].x_v[i][j] = radius * atan(-sin(cswm[p].lon_v[i][j]) / tan(cswm[p].lat_v[i][j]));
                    cswm[p].y_v[i][j] = radius * atan(-cos(cswm[p].lon_v[i][j]) / tan(cswm[p].lat_v[i][j]));

                    // A/AInverse
                    cswm[p].A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][1] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon[i][j]));

                    cswm[p].AInverse[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon[i][j]));
                    cswm[p].AInverse[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon[i][j]));

                    // A_u/AInverse_u
                    cswm[p].A_u[i][j][0] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * (gamma_u[i][j] * cos(beta2D[i][j]) / cos(alpha2D_u[i][j]) * cos(cswm[p].lon_u[i][j]));
                    cswm[p].A_u[i][j][1] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * (-gamma_u[i][j] * cos(alpha2D_u[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon_u[i][j]));
                    cswm[p].A_u[i][j][2] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * (cos(beta2D[i][j]) / cos(alpha2D_u[i][j]) * sin(cswm[p].lon_u[i][j]));
                    cswm[p].A_u[i][j][3] = 1 / (pow(gamma_u[i][j], 2) * cos(alpha2D_u[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D_u[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon_u[i][j]));

                    cswm[p].AInverse_u[i][j][0] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D_u[i][j]) / cos(beta2D[i][j]) * cos(cswm[p].lon_u[i][j]));
                    cswm[p].AInverse_u[i][j][1] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (gamma_u[i][j] * cos(alpha2D_u[i][j]) / cos(beta2D[i][j]) * sin(cswm[p].lon_u[i][j]));
                    cswm[p].AInverse_u[i][j][2] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (-cos(beta2D[i][j]) / cos(alpha2D_u[i][j]) * sin(cswm[p].lon_u[i][j]));
                    cswm[p].AInverse_u[i][j][3] = gamma_u[i][j] * cos(alpha2D_u[i][j]) * cos(beta2D[i][j]) * (gamma_u[i][j] * cos(beta2D[i][j]) / cos(alpha2D_u[i][j]) * cos(cswm[p].lon_u[i][j]));

                    // A_v/AInverse_v
                    cswm[p].A_v[i][j][0] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * (gamma_v[i][j] * cos(beta2D_v[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon_v[i][j]));
                    cswm[p].A_v[i][j][1] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * (-gamma_v[i][j] * cos(alpha2D[i][j]) / cos(beta2D_v[i][j]) * sin(cswm[p].lon_v[i][j]));
                    cswm[p].A_v[i][j][2] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * (cos(beta2D_v[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon_v[i][j]));
                    cswm[p].A_v[i][j][3] = 1 / (pow(gamma_v[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D_v[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D_v[i][j]) * cos(cswm[p].lon_v[i][j]));

                    cswm[p].AInverse_v[i][j][0] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D_v[i][j]) * cos(cswm[p].lon_v[i][j]));
                    cswm[p].AInverse_v[i][j][1] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (gamma_v[i][j] * cos(alpha2D[i][j]) / cos(beta2D_v[i][j]) * sin(cswm[p].lon_v[i][j]));
                    cswm[p].AInverse_v[i][j][2] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (-cos(beta2D_v[i][j]) / cos(alpha2D[i][j]) * sin(cswm[p].lon_v[i][j]));
                    cswm[p].AInverse_v[i][j][3] = gamma_v[i][j] * cos(alpha2D[i][j]) * cos(beta2D_v[i][j]) * (gamma_v[i][j] * cos(beta2D_v[i][j]) / cos(alpha2D[i][j]) * cos(cswm[p].lon_v[i][j]));
                }
            }
        }
    }


    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                cswm[p].CF[i][j] = 2 * omega * (-cos(cswm[p].lon[i][j]) * cos(cswm[p].lat[i][j]) * sin(ALPHA0) + sin(cswm[p].lat[i][j]) * cos(ALPHA0));
                cswm[p].CF_u[i][j] = 2 * omega * (-cos(cswm[p].lon_u[i][j]) * cos(cswm[p].lat_u[i][j]) * sin(ALPHA0) + sin(cswm[p].lat_u[i][j]) * cos(ALPHA0));
                cswm[p].CF_v[i][j] = 2 * omega * (-cos(cswm[p].lon_v[i][j]) * cos(cswm[p].lat_v[i][j]) * sin(ALPHA0) + sin(cswm[p].lat_v[i][j]) * cos(ALPHA0));
            }
        }
    }

    for (int idx = 0; idx < NX; idx++) {
        std::cout << cswm[0].lat_u[NX-2][idx] << " " << cswm[3].lat_u[0][idx];
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void CSWM::BoundaryProcess(CSWM &model) {
    // patch 1
    for (int idx = 0; idx < NX; idx++) {
        if (idx < NX / 2) {
            model.cswm[0].hp[0][idx] = Interpolate4lat(model.cswm[3].lat[NX-2][idx], model.cswm[3].lat[NX-2][idx+1], model.cswm[3].hp[NX-2][idx], model.cswm[3].hp[NX-2][idx+1], model.cswm[0].lat[0][idx]);
            model.cswm[0].hp[NX-1][idx] = Interpolate4lat(model.cswm[1].lat[1][idx], model.cswm[1].lat[1][idx+1], model.cswm[1].hp[1][idx], model.cswm[1].hp[1][idx+1], model.cswm[0].lat[NX-1][idx]);
            model.cswm[0].hp[idx][NY-1] = Interpolate4lat(model.cswm[4].lon[idx][1], model.cswm[4].lon[idx+1][1], model.cswm[4].hp[idx][1], model.cswm[4].hp[idx+1][1], model.cswm[0].lon[idx][NY-1]);
            model.cswm[0].hp[idx][0] = Interpolate4lat(model.cswm[5].lon[idx][NY-2], model.cswm[5].lon[idx+1][NY-2], model.cswm[5].hp[idx][NY-2], model.cswm[5].hp[idx+1][NY-2], model.cswm[0].lon[idx][0]);
        }
        else if (idx == NX / 2) {
            model.cswm[0].hp[0][idx] = model.cswm[3].hp[NX-2][idx];
            model.cswm[0].hp[NX-1][idx] = model.cswm[1].hp[1][idx];
            model.cswm[0].hp[idx][NY-1] = model.cswm[4].hp[idx][1];
            model.cswm[0].hp[idx][0] = model.cswm[5].hp[idx][NY-2];
        }
        else {
            model.cswm[0].hp[0][idx] = Interpolate4lat(model.cswm[3].lat[NX-2][idx-1], model.cswm[3].lat[NX-2][idx], model.cswm[3].hp[NX-2][idx-1], model.cswm[3].hp[NX-2][idx], model.cswm[0].lat[0][idx]);
            model.cswm[0].hp[NX-1][idx] = Interpolate4lat(model.cswm[1].lat[1][idx-1], model.cswm[1].lat[1][idx], model.cswm[1].hp[1][idx-1], model.cswm[1].hp[1][idx], model.cswm[0].lat[NX-1][idx]);
            model.cswm[0].hp[idx][NY-1] = Interpolate4lat(model.cswm[4].lon[idx-1][1], model.cswm[4].lon[idx][1], model.cswm[4].hp[idx-1][1], model.cswm[4].hp[idx][1], model.cswm[0].lon[idx][NY-1]);
            model.cswm[0].hp[idx][0] = Interpolate4lat(model.cswm[5].lon[idx-1][NY-2], model.cswm[5].lon[idx][NY-2], model.cswm[5].hp[idx-1][NY-2], model.cswm[5].hp[idx][NY-2], model.cswm[0].lon[idx][0]);
        }
    }

    // patch 2
    for (int idx = 0; idx < NX; idx++) {
        if (idx < NX / 2) {
            model.cswm[1].hp[0][idx] = Interpolate4lat(model.cswm[0].lat[NX-2][idx], model.cswm[0].lat[NX-2][idx+1], model.cswm[0].hp[NX-2][idx], model.cswm[0].hp[NX-2][idx+1], model.cswm[1].lat[0][idx]);
            model.cswm[1].hp[NX-1][idx] = Interpolate4lat(model.cswm[2].lat[1][idx], model.cswm[2].lat[1][idx+1], model.cswm[2].hp[1][idx], model.cswm[2].hp[1][idx+1], model.cswm[1].lat[NX-1][idx]);
            model.cswm[1].hp[idx][NY-1] = Interpolate4lat(model.cswm[4].lon[NX-2][idx], model.cswm[4].lon[NX-2][idx+1], model.cswm[4].hp[NX-2][idx], model.cswm[4].hp[NX-2][idx+1], model.cswm[1].lon[idx][NY-1]);
            model.cswm[1].hp[idx][0] = Interpolate4lat(model.cswm[5].lon[NX-2][(NY-1)-idx], model.cswm[5].lon[NX-2][(NY-1)-idx-1], model.cswm[5].hp[NX-2][(NY-1)-idx], model.cswm[5].hp[NX-2][(NY-1)-idx-1], model.cswm[1].lon[idx][0]);
        }
        else if (idx == NX / 2) {
            model.cswm[1].hp[0][idx] = model.cswm[0].hp[NX-2][idx];
            model.cswm[1].hp[NX-1][idx] = model.cswm[2].hp[1][idx];
            model.cswm[1].hp[idx][NY-1] = model.cswm[4].hp[NX-2][idx];
            model.cswm[1].hp[idx][0] = model.cswm[5].hp[NX-2][(NY-1)-idx];
        }
        else {
            model.cswm[1].hp[0][idx] = Interpolate4lat(model.cswm[0].lat[NX-2][idx-1], model.cswm[0].lat[NX-2][idx], model.cswm[0].hp[NX-2][idx-1], model.cswm[0].hp[NX-2][idx], model.cswm[1].lat[0][idx]);
            model.cswm[1].hp[NX-1][idx] = Interpolate4lat(model.cswm[2].lat[1][idx-1], model.cswm[2].lat[1][idx], model.cswm[2].hp[1][idx-1], model.cswm[2].hp[1][idx], model.cswm[1].lat[NX-1][idx]);
            model.cswm[1].hp[idx][NY-1] = Interpolate4lat(model.cswm[4].lon[NX-2][idx-1], model.cswm[4].lon[NX-2][idx], model.cswm[4].hp[NX-2][idx-1], model.cswm[4].hp[NX-2][idx], model.cswm[1].lon[idx][NY-1]);
            model.cswm[1].hp[idx][0] = Interpolate4lat(model.cswm[5].lon[NX-2][(NY-1)-idx+1], model.cswm[5].lon[NX-2][(NY-1)-idx], model.cswm[5].hp[NX-2][(NY-1)-idx+1], model.cswm[5].hp[NX-2][(NY-1)-idx], model.cswm[1].lon[idx][0]);
        }
        
    }

    // patch 3
    for (int idx = 0; idx < NX; idx++) {
        // left
        if (idx < NX / 2) {
            model.cswm[2].hp[0][idx] = Interpolate4lat(model.cswm[1].lat[NX-2][idx], model.cswm[1].lat[NX-2][idx+1], model.cswm[1].hp[NX-2][idx], model.cswm[1].hp[NX-2][idx+1], model.cswm[2].lat[0][idx]);
            model.cswm[2].hp[NX-1][idx] = Interpolate4lat(model.cswm[3].lat[1][idx], model.cswm[3].lat[1][idx+1], model.cswm[3].hp[1][idx], model.cswm[3].hp[1][idx+1], model.cswm[2].lat[NX-1][idx]);
            model.cswm[2].hp[idx][NY-1] = Interpolate(model.cswm[4].lon[(NX-1)-idx][NY-2], model.cswm[4].lon[(NX-1)-idx-1][NY-2], model.cswm[4].hp[(NX-1)-idx][NY-2], model.cswm[4].hp[(NX-1)-idx-1][NY-2], model.cswm[2].lon[idx][NY-1]);
            model.cswm[2].hp[idx][0] = Interpolate(model.cswm[5].lon[(NX-1)-idx][1], model.cswm[5].lon[(NX-1)-idx-1][1], model.cswm[5].hp[(NX-1)-idx][1], model.cswm[5].hp[(NX-1)-idx-1][1], model.cswm[2].lon[idx][0]);
        }
        else if (idx == NX / 2) {
            model.cswm[2].hp[0][idx] = model.cswm[1].hp[NX-2][idx];
            model.cswm[2].hp[NX-1][idx] = model.cswm[3].hp[1][idx];
            model.cswm[2].hp[idx][NY-1] = model.cswm[4].hp[(NX-1)-idx][NY-2];
            model.cswm[2].hp[idx][0] = model.cswm[5].hp[(NX-1)-idx][1];
        }
        else {
            model.cswm[2].hp[0][idx] = Interpolate4lat(model.cswm[1].lat[NX-2][idx-1], model.cswm[1].lat[NX-2][idx], model.cswm[1].hp[NX-2][idx-1], model.cswm[1].hp[NX-2][idx], model.cswm[2].lat[0][idx]);
            model.cswm[2].hp[NX-1][idx] = Interpolate4lat(model.cswm[3].lat[1][idx-1], model.cswm[3].lat[1][idx], model.cswm[3].hp[1][idx-1], model.cswm[3].hp[1][idx], model.cswm[2].lat[NX-1][idx]);
            model.cswm[2].hp[idx][NY-1] = Interpolate(model.cswm[4].lon[(NX-1)-idx+1][NY-2], model.cswm[4].lon[(NX-1)-idx][NY-2], model.cswm[4].hp[(NX-1)-idx+1][NY-2], model.cswm[4].hp[(NX-1)-idx][NY-2], model.cswm[2].lon[idx][NY-1]);
            model.cswm[2].hp[idx][0] = Interpolate(model.cswm[5].lon[(NX-1)-idx+1][1], model.cswm[5].lon[(NX-1)-idx][1], model.cswm[5].hp[(NX-1)-idx+1][1], model.cswm[5].hp[(NX-1)-idx][1], model.cswm[2].lon[idx][0]);
        }
        
    }

    // patch 4
    for (int idx = 0; idx < NX; idx++) {
        // left
        if (idx < NX / 2) {
            model.cswm[3].hp[0][idx] = Interpolate4lat(model.cswm[2].lat[NX-2][idx], model.cswm[2].lat[NX-2][idx+1], model.cswm[2].hp[NX-2][idx], model.cswm[2].hp[NX-2][idx+1], model.cswm[3].lat[0][idx]);
            model.cswm[3].hp[NX-1][idx] = Interpolate4lat(model.cswm[0].lat[1][idx], model.cswm[0].lat[1][idx+1], model.cswm[0].hp[1][idx], model.cswm[0].hp[1][idx+1], model.cswm[3].lat[NX-1][idx]);
            model.cswm[3].hp[idx][NY-1] = Interpolate(model.cswm[4].lon[1][(NY-1)-idx], model.cswm[4].lon[1][(NY-1)-idx-1], model.cswm[4].hp[1][(NY-1)-idx], model.cswm[4].hp[1][(NY-1)-idx-1], model.cswm[3].lon[idx][NY-1]);
            model.cswm[3].hp[idx][0] = Interpolate(model.cswm[5].lon[1][idx], model.cswm[5].lon[1][idx+1], model.cswm[5].hp[1][idx], model.cswm[5].hp[1][idx+1], model.cswm[3].lon[idx][0]);
        }
        else if (idx == NX / 2) {
            model.cswm[3].hp[0][idx] = model.cswm[2].hp[NX-2][idx];
            model.cswm[3].hp[NX-1][idx] = model.cswm[0].hp[1][idx];
            model.cswm[3].hp[idx][NY-1] = model.cswm[4].hp[1][(NY-1)-idx];
            model.cswm[3].hp[idx][0] = model.cswm[5].hp[1][idx];
        }
        else {
            model.cswm[3].hp[0][idx] = Interpolate4lat(model.cswm[2].lat[NX-2][idx-1], model.cswm[2].lat[NX-2][idx], model.cswm[2].hp[NX-2][idx-1], model.cswm[2].hp[NX-2][idx], model.cswm[3].lat[0][idx]);
            model.cswm[3].hp[NX-1][idx] = Interpolate4lat(model.cswm[0].lat[1][idx-1], model.cswm[0].lat[1][idx], model.cswm[0].hp[1][idx-1], model.cswm[0].hp[1][idx], model.cswm[3].lat[NX-1][idx]);
            model.cswm[3].hp[idx][NY-1] = Interpolate(model.cswm[4].lon[1][(NY-1)-idx+1], model.cswm[4].lon[1][(NY-1)-idx], model.cswm[4].hp[1][(NY-1)-idx+1], model.cswm[4].hp[1][(NY-1)-idx], model.cswm[3].lon[idx][NY-1]);
            model.cswm[3].hp[idx][0] = Interpolate(model.cswm[5].lon[1][idx-1], model.cswm[5].lon[1][idx], model.cswm[5].hp[1][idx-1], model.cswm[5].hp[1][idx], model.cswm[3].lon[idx][0]);
        }
    }

    // patch 5
    for (int idx = 0; idx < NX; idx++) {
        // left
        if (idx < NX / 2) {
            model.cswm[4].hp[0][idx] = Interpolate(model.cswm[3].lon[(NX-1)-idx][NY-2], model.cswm[3].lon[(NX-1)-idx-1][NY-2], model.cswm[3].hp[(NX-1)-idx][NY-2], model.cswm[3].hp[(NX-1)-idx-1][NY-2], model.cswm[4].lon[0][idx]);
            model.cswm[4].hp[NX-1][idx] = Interpolate4lat(model.cswm[1].lon[idx][NY-2], model.cswm[1].lon[idx+1][NY-2], model.cswm[1].hp[idx][NY-2], model.cswm[1].hp[idx+1][NY-2], model.cswm[4].lon[NX-1][idx]);
            model.cswm[4].hp[idx][NY-1] = Interpolate(model.cswm[2].lon[(NX-1)-idx][NY-2], model.cswm[2].lon[(NX-1)-idx-1][NY-2], model.cswm[2].hp[(NX-1)-idx][NY-2], model.cswm[2].hp[(NX-1)-idx-1][NY-2], model.cswm[4].lon[idx][NY-1]);
            model.cswm[4].hp[idx][0] = Interpolate4lat(model.cswm[0].lon[idx][NY-2], model.cswm[0].lon[idx+1][NY-2], model.cswm[0].hp[idx][NY-2], model.cswm[0].hp[idx+1][NY-2], model.cswm[4].lon[idx][0]);
        }
        else if (idx == NX / 2) {
            model.cswm[4].hp[0][idx] = model.cswm[3].hp[(NX-1)-idx][NY-2];
            model.cswm[4].hp[NX-1][idx] = model.cswm[1].hp[idx][NY-2];
            model.cswm[4].hp[idx][NY-1] = model.cswm[2].hp[(NX-1)-idx][NY-2];
            model.cswm[4].hp[idx][0] = model.cswm[0].hp[idx][NY-2];
        }
        else {
            model.cswm[4].hp[0][idx] = Interpolate(model.cswm[3].lon[(NX-1)-idx+1][NY-2], model.cswm[3].lon[(NX-1)-idx][NY-2], model.cswm[3].hp[(NX-1)-idx+1][NY-2], model.cswm[3].hp[(NX-1)-idx][NY-2], model.cswm[4].lon[0][idx]);
            model.cswm[4].hp[NX-1][idx] = Interpolate4lat(model.cswm[1].lon[idx-1][NY-2], model.cswm[1].lon[idx][NY-2], model.cswm[1].hp[idx-1][NY-2], model.cswm[1].hp[idx][NY-2], model.cswm[4].lon[NX-1][idx]);
            model.cswm[4].hp[idx][NY-1] = Interpolate(model.cswm[2].lon[(NX-1)-idx+1][NY-2], model.cswm[2].lon[(NX-1)-idx][NY-2], model.cswm[2].hp[(NX-1)-idx+1][NY-2], model.cswm[2].hp[(NX-1)-idx][NY-2], model.cswm[4].lon[idx][NY-1]);
            model.cswm[4].hp[idx][0] = Interpolate4lat(model.cswm[0].lon[idx-1][NY-2], model.cswm[0].lon[idx][NY-2], model.cswm[0].hp[idx-1][NY-2], model.cswm[0].hp[idx][NY-2], model.cswm[4].lon[idx][0]);
        }
    }
    
    // patch 6
    for (int idx = 0; idx < NX; idx++) {
        // left
        if (idx < NX / 2) {
            model.cswm[5].hp[0][idx] = Interpolate(model.cswm[3].lon[idx][1], model.cswm[3].lon[idx+1][1], model.cswm[3].hp[idx][1], model.cswm[3].hp[idx+1][1], model.cswm[5].lon[0][idx]);
            model.cswm[5].hp[NX-1][idx] = Interpolate4lat(model.cswm[1].lon[(NX-1)-idx][1], model.cswm[1].lon[(NX-1)-idx-1][1], model.cswm[1].hp[(NX-1)-idx][1], model.cswm[1].hp[(NX-1)-idx-1][1], model.cswm[5].lon[NX-1][idx]);
            model.cswm[5].hp[idx][NY-1] = Interpolate4lat(model.cswm[0].lon[idx][1], model.cswm[0].lon[idx+1][1], model.cswm[0].hp[idx][1], model.cswm[0].hp[idx+1][1], model.cswm[5].lon[idx][NY-1]);
            model.cswm[5].hp[idx][0] = Interpolate(model.cswm[2].lon[(NX-1)-idx][1], model.cswm[2].lon[(NX-1)-idx+1][1], model.cswm[2].hp[(NX-1)-idx][1], model.cswm[2].hp[(NX-1)-idx+1][1], model.cswm[5].lon[idx][0]);
        }
        else if (idx == NX / 2) {
            model.cswm[5].hp[0][idx] = model.cswm[3].hp[idx][1];
            model.cswm[5].hp[NX-1][idx] = model.cswm[1].hp[(NX-1)-idx][1];
            model.cswm[5].hp[idx][NY-1] = model.cswm[0].hp[idx][1];
            model.cswm[5].hp[idx][0] = model.cswm[2].hp[(NX-1)-idx][1];
        }
        else {
            model.cswm[5].hp[0][idx] = Interpolate(model.cswm[3].lon[idx-1][1], model.cswm[3].lon[idx][1], model.cswm[3].hp[idx-1][1], model.cswm[3].hp[idx][1], model.cswm[5].lon[0][idx]);
            model.cswm[5].hp[NX-1][idx] = Interpolate4lat(model.cswm[1].lon[(NX-1)-idx+1][1], model.cswm[1].lon[(NX-1)-idx][1], model.cswm[1].hp[(NX-1)-idx+1][1], model.cswm[1].hp[(NX-1)-idx][1], model.cswm[5].lon[NX-1][idx]);
            model.cswm[5].hp[idx][NY-1] = Interpolate4lat(model.cswm[0].lon[idx-1][1], model.cswm[0].lon[idx][1], model.cswm[0].hp[idx-1][1], model.cswm[0].hp[idx][1], model.cswm[5].lon[idx][NY-1]);
            model.cswm[5].hp[idx][0] = Interpolate(model.cswm[2].lon[(NX-1)-idx+1][1], model.cswm[2].lon[(NX-1)-idx][1], model.cswm[2].hp[(NX-1)-idx+1][1], model.cswm[2].hp[(NX-1)-idx][1], model.cswm[5].lon[idx][0]);
        }
    }
}

double CSWM::Interpolate(double A1, double A2, double y1, double y2, double B) {
    if (A1 <0) A1 += 2 * M_PI;
    if (A2 < 0) A1 += 2 * M_PI;
    if (B < 0) B += 2 * M_PI;
    return y1 + (y2-y1) * (B-A1) / (A2 - A1);
}

double CSWM::Interpolate4lat(double A1, double A2, double y1, double y2, double B) {
    return y1 + (y2-y1) * (B-A1) / (A2 - A1);
}

/*
void CSWM::BoundaryTransform(CSWM &model) {
    // patch 1
    for (int idx = 0; idx < NX; idx++) {
        // left
        // model.cswm[0].up[0][idx] = ConvertUPatch2Patch(model, 0, 3, 0, NX-2, idx, idx);
        // model.cswm[0].vp[0][idx] = ConvertVPatch2Patch(model, 0, 3, 0, NX-2, idx, idx);
        model.cswm[0].up[0][idx] = model.cswm[3].up[NX-2][idx];
        model.cswm[0].vp[0][idx] = model.cswm[3].vp[NX-2][idx];

        // right
        // model.cswm[0].up[NX-1][idx] = ConvertUPatch2Patch(model, 0, 1, NX-1, 1, idx, idx);
        // model.cswm[0].vp[NX-1][idx] = ConvertVPatch2Patch(model, 0, 1, NX-1, 1, idx, idx);
        model.cswm[0].up[NX-1][idx] = model.cswm[1].up[1][idx];
        model.cswm[0].vp[NX-1][idx] = model.cswm[1].vp[1][idx];

        // up
        // model.cswm[0].up[idx][NY-1] = ConvertUPatch2Patch(model, 0, 4, idx, idx, NY-1, 1);
        // model.cswm[0].vp[idx][NY-1] = ConvertVPatch2Patch(model, 0, 4, idx, idx, NY-1, 1);
        model.cswm[0].up[idx][NY-1] = model.cswm[4].up[idx][1];
        model.cswm[0].vp[idx][NY-1] = model.cswm[4].vp[idx][1];

        // bottom
        // model.cswm[0].up[idx][0] = ConvertUPatch2Patch(model, 0, 5, idx, idx, 0, NY-2);
        // model.cswm[0].vp[idx][0] = ConvertVPatch2Patch(model, 0, 5, idx, idx, 0, NY-2);
        model.cswm[0].up[idx][0] = model.cswm[5].up[idx][NY-2];
        model.cswm[0].vp[idx][0] = model.cswm[5].vp[idx][NY-2];
    }

    // patch 2
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        // model.cswm[1].up[0][idx] = ConvertUPatch2Patch(model, 1, 0, 0, NX-2, idx, idx);
        // model.cswm[1].vp[0][idx] = ConvertVPatch2Patch(model, 1, 0, 0, NX-2, idx, idx);
        model.cswm[1].up[0][idx] = model.cswm[0].up[NX-2][idx];
        model.cswm[1].vp[0][idx] = model.cswm[0].vp[NX-2][idx];

        // right
        // model.cswm[1].up[NX-1][idx] = ConvertUPatch2Patch(model, 1, 2, NX-1, 1, idx, idx);
        // model.cswm[1].vp[NX-1][idx] = ConvertVPatch2Patch(model, 1, 2, NX-1, 1, idx, idx);
        model.cswm[1].up[NX-1][idx] = model.cswm[2].up[1][idx];
        model.cswm[1].vp[NX-1][idx] = model.cswm[2].vp[1][idx];

        // up
        // model.cswm[1].up[idx][NY-1] = ConvertUPatch2Patch(model, 1, 4, idx, NX-2, NY-1, idx);
        // model.cswm[1].vp[idx][NY-1] = ConvertVPatch2Patch(model, 1, 4, idx, NX-2, NY-1, idx);
        // model.cswm[1].up[idx][NY-1] = ConvertBV2AUPatch2Patch(model, 1, 4, idx, NX-2, NY-1, idx);
        // model.cswm[1].vp[idx][NY-1] = ConvertBU2AVPatch2Patch(model, 1, 4, idx, NX-2, NY-1, idx);
        model.cswm[1].up[idx][NY-1] = model.cswm[4].vp[NX-2][idx];
        model.cswm[1].vp[idx][NY-1] = -model.cswm[4].up[NX-2][idx];

        // bottom
        // model.cswm[1].up[idx][0] = ConvertUPatch2Patch(model, 1, 5, idx, NX-2, 0, (NY-1)-idx);
        // model.cswm[1].vp[idx][0] = ConvertVPatch2Patch(model, 1, 5, idx, NX-2, 0, (NY-1)-idx);
        // model.cswm[1].up[idx][0] = ConvertBV2AUPatch2Patch(model, 1, 5, idx, NX-2, 0, (NY-1)-idx);
        // model.cswm[1].vp[idx][0] = ConvertBU2AVPatch2Patch(model, 1, 5, idx, NX-2, 0, (NY-1)-idx);
        model.cswm[1].up[idx][0] = -model.cswm[5].vp[NX-2][(NY-1)-idx];
        model.cswm[1].vp[idx][0] = model.cswm[5].up[NX-2][(NY-1)-idx];
    }

    // patch 3
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        // model.cswm[2].up[0][idx] = ConvertUPatch2Patch(model, 2, 1, 0, NX-2, idx, idx);
        // model.cswm[2].vp[0][idx] = ConvertVPatch2Patch(model, 2, 1, 0, NX-2, idx, idx);
        model.cswm[2].up[0][idx] = model.cswm[1].up[NX-2][idx];
        model.cswm[2].vp[0][idx] = model.cswm[1].vp[NX-2][idx];

        // right
        // model.cswm[2].up[NX-1][idx] = ConvertUPatch2Patch(model, 2, 3, NX-1, 1, idx, idx);
        // model.cswm[2].vp[NX-1][idx] = ConvertVPatch2Patch(model, 2, 3, NX-1, 1, idx, idx);
        model.cswm[2].up[NX-1][idx] = model.cswm[3].up[1][idx];
        model.cswm[2].vp[NX-1][idx] = model.cswm[3].vp[1][idx];

        // up
        // model.cswm[2].up[idx][NY-1] = ConvertUPatch2Patch(model, 2, 4, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[2].vp[idx][NY-1] = ConvertVPatch2Patch(model, 2, 4, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[2].up[idx][NY-1] = ConvertUPatch2Patch(model, 2, 4, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[2].vp[idx][NY-1] = ConvertVPatch2Patch(model, 2, 4, idx, (NX-1)-idx, NY-1, NY-2);
        model.cswm[2].up[idx][NY-1] = -model.cswm[4].up[(NX-1)-idx][NY-2];
        model.cswm[2].vp[idx][NY-1] = -model.cswm[4].vp[(NX-1)-idx][NY-2];

        // bottom
        // model.cswm[2].up[idx][0] = ConvertUPatch2Patch(model, 2, 5, idx, (NX-1)-idx, 0, 1);
        // model.cswm[2].vp[idx][0] = ConvertVPatch2Patch(model, 2, 5, idx, (NX-1)-idx, 0, 1);
        // model.cswm[2].up[idx][0] = ConvertUPatch2Patch(model, 2, 5, idx, (NX-1)-idx, 0, 1);
        // model.cswm[2].vp[idx][0] = ConvertVPatch2Patch(model, 2, 5, idx, (NX-1)-idx, 0, 1);
        model.cswm[2].up[idx][0] = -model.cswm[5].up[(NX-1)-idx][1];
        model.cswm[2].vp[idx][0] = -model.cswm[5].vp[(NX-1)-idx][1];
    }

    // patch 4
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        // model.cswm[3].up[0][idx] = ConvertUPatch2Patch(model, 3, 2, 0, NX-2, idx, idx);
        // model.cswm[3].vp[0][idx] = ConvertVPatch2Patch(model, 3, 2, 0, NX-2, idx, idx);
        model.cswm[3].up[0][idx] = model.cswm[2].up[NX-2][idx];
        model.cswm[3].vp[0][idx] = model.cswm[2].vp[NX-2][idx];

        // right
        // model.cswm[3].up[NX-1][idx] = ConvertUPatch2Patch(model, 3, 0, NX-1, 1, idx, idx);
        // model.cswm[3].vp[NX-1][idx] = ConvertVPatch2Patch(model, 3, 0, NX-1, 1, idx, idx);
        model.cswm[3].up[NX-1][idx] = model.cswm[0].up[1][idx];
        model.cswm[3].vp[NX-1][idx] = model.cswm[0].vp[1][idx];

        // up
        // model.cswm[3].up[idx][NY-1] = ConvertUPatch2Patch(model, 3, 4, idx, 1, NY-1, (NY-1)-idx);
        // model.cswm[3].vp[idx][NY-1] = ConvertVPatch2Patch(model, 3, 4, idx, 1, NY-1, (NY-1)-idx);
        // model.cswm[3].up[idx][NY-1] = ConvertBV2AUPatch2Patch(model, 3, 4, idx, 1, NY-1, (NY-1)-idx);
        // model.cswm[3].vp[idx][NY-1] = ConvertBU2AVPatch2Patch(model, 3, 4, idx, 1, NY-1, (NY-1)-idx);
        model.cswm[3].up[idx][NY-1] = -model.cswm[4].vp[1][(NY-1)-idx];
        model.cswm[3].vp[idx][NY-1] = model.cswm[4].up[1][(NY-1)-idx];

        // bottom
        // model.cswm[3].up[idx][0] = ConvertUPatch2Patch(model, 3, 5, idx, 1, 0, idx);
        // model.cswm[3].vp[idx][0] = ConvertVPatch2Patch(model, 3, 5, idx, 1, 0, idx);
        // model.cswm[3].up[idx][0] = ConvertBV2AUPatch2Patch(model, 3, 5, idx, 1, 0, idx);
        // model.cswm[3].vp[idx][0] = ConvertBU2AVPatch2Patch(model, 3, 5, idx, 1, 0, idx);
        model.cswm[3].up[idx][0] = model.cswm[5].vp[1][idx];
        model.cswm[3].vp[idx][0] = -model.cswm[5].up[1][idx];
    }

    // patch 5
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        // model.cswm[4].up[0][idx] = ConvertUPatch2Patch(model, 4, 3, 0, (NX-1)-idx, idx, NY-2);
        // model.cswm[4].vp[0][idx] = ConvertVPatch2Patch(model, 4, 3, 0, (NX-1)-idx, idx, NY-2);
        // model.cswm[4].up[0][idx] = ConvertBV2AUPatch2Patch(model, 4, 3, 0, (NX-1)-idx, idx, NY-2);
        // model.cswm[4].vp[0][idx] = ConvertBU2AVPatch2Patch(model, 4, 3, 0, (NX-1)-idx, idx, NY-2);
        model.cswm[4].up[0][idx] = model.cswm[3].vp[(NX-1)-idx][NY-2];
        model.cswm[4].vp[0][idx] = -model.cswm[3].up[(NX-1)-idx][NY-2];

        // right
        // model.cswm[4].up[NX-1][idx] = ConvertUPatch2Patch(model, 4, 1, NX-1, idx, idx, NY-2);
        // model.cswm[4].vp[NX-1][idx] = ConvertVPatch2Patch(model, 4, 1, NX-1, idx, idx, NY-2);
        // model.cswm[4].up[NX-1][idx] = ConvertBV2AUPatch2Patch(model, 4, 1, NX-1, idx, idx, NY-2);
        // model.cswm[4].vp[NX-1][idx] = ConvertBU2AVPatch2Patch(model, 4, 1, NX-1, idx, idx, NY-2);
        model.cswm[4].up[NX-1][idx] = -model.cswm[1].vp[idx][NY-2];
        model.cswm[4].vp[NX-1][idx] = model.cswm[1].up[idx][NY-2];

        // up
        // model.cswm[4].up[idx][NY-1] = ConvertUPatch2Patch(model, 4, 2, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[4].vp[idx][NY-1] = ConvertVPatch2Patch(model, 4, 2, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[4].up[idx][NY-1] = ConvertUPatch2Patch(model, 4, 2, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[4].vp[idx][NY-1] = ConvertVPatch2Patch(model, 4, 2, idx, (NX-1)-idx, NY-1, NY-2);
        model.cswm[4].up[idx][NY-1] = -model.cswm[2].up[(NX-1)-idx][NY-2];
        model.cswm[4].vp[idx][NY-1] = -model.cswm[2].vp[(NX-1)-idx][NY-2];

        // bottom
        // model.cswm[4].up[idx][0] = ConvertUPatch2Patch(model, 4, 0, idx, idx, 0, NY-2);
        // model.cswm[4].vp[idx][0] = ConvertVPatch2Patch(model, 4, 0, idx, idx, 0, NY-2);
        model.cswm[4].up[idx][0] = model.cswm[0].up[idx][NY-2];
        model.cswm[4].vp[idx][0] = model.cswm[0].vp[idx][NY-2];
    }

    // patch 6
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        // model.cswm[5].up[0][idx] = ConvertUPatch2Patch(model, 5, 3, 0, idx, idx, 1);
        // model.cswm[5].vp[0][idx] = ConvertVPatch2Patch(model, 5, 3, 0, idx, idx, 1);
        // model.cswm[5].up[0][idx] = ConvertBV2AUPatch2Patch(model, 5, 3, 0, idx, idx, 1);
        // model.cswm[5].vp[0][idx] = ConvertBU2AVPatch2Patch(model, 5, 3, 0, idx, idx, 1);
        model.cswm[5].up[0][idx] = -model.cswm[3].vp[idx][1];
        model.cswm[5].vp[0][idx] = model.cswm[3].up[idx][1];

        // right
        // model.cswm[5].up[NX-1][idx] = ConvertUPatch2Patch(model, 5, 1, NX-1, (NX-1)-idx, idx, 1);
        // model.cswm[5].vp[NX-1][idx] = ConvertVPatch2Patch(model, 5, 1, NX-1, (NX-1)-idx, idx, 1);
        // model.cswm[5].up[NX-1][idx] = ConvertBV2AUPatch2Patch(model, 5, 1, NX-1, (NX-1)-idx, idx, 1);
        // model.cswm[5].vp[NX-1][idx] = ConvertBU2AVPatch2Patch(model, 5, 1, NX-1, (NX-1)-idx, idx, 1);
        model.cswm[5].up[NX-1][idx] = model.cswm[1].vp[(NX-1)-idx][1];
        model.cswm[5].vp[NX-1][idx] = -model.cswm[1].up[(NX-1)-idx][1];

        // up
        // model.cswm[5].up[idx][NY-1] = ConvertUPatch2Patch(model, 5, 0, idx, idx, NY-1, 1);
        // model.cswm[5].vp[idx][NY-1] = ConvertVPatch2Patch(model, 5, 0, idx, idx, NY-1, 1);
        model.cswm[5].up[idx][NY-1] = model.cswm[0].up[idx][1];
        model.cswm[5].vp[idx][NY-1] = model.cswm[0].vp[idx][1];

        // bottom
        // model.cswm[5].up[idx][0] = ConvertUPatch2Patch(model, 5, 2, idx, (NX-1)-idx, 0, 1);
        // model.cswm[5].vp[idx][0] = ConvertVPatch2Patch(model, 5, 2, idx, (NX-1)-idx, 0, 1);
        // model.cswm[5].up[idx][0] = ConvertUPatch2Patch(model, 5, 2, idx, (NX-1)-idx, 0, 1);
        // model.cswm[5].vp[idx][0] = ConvertVPatch2Patch(model, 5, 2, idx, (NX-1)-idx, 0, 1);
        model.cswm[5].up[idx][0] = -model.cswm[2].up[(NX-1)-idx][1];
        model.cswm[5].vp[idx][0] = -model.cswm[2].vp[(NX-1)-idx][1];
    }
    return;
}
*/

void CSWM::BoundaryTransform(CSWM &model) {
    // patch 1
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[0].up[0][idx] = ConvertUPatch2Patch(model, 0, 3, 0, NX-2, idx, idx);
        model.cswm[0].vp[0][idx] = ConvertVPatch2Patch(model, 0, 3, 0, NX-2, idx, idx);
        // model.cswm[0].up[0][idx] = model.cswm[3].up[NX-2][idx];
        // model.cswm[0].vp[0][idx] = model.cswm[3].vp[NX-2][idx];

        // right
        model.cswm[0].up[NX-1][idx] = ConvertUPatch2Patch(model, 0, 1, NX-1, 1, idx, idx);
        model.cswm[0].vp[NX-1][idx] = ConvertVPatch2Patch(model, 0, 1, NX-1, 1, idx, idx);
        // model.cswm[0].up[NX-1][idx] = model.cswm[1].up[1][idx];
        // model.cswm[0].vp[NX-1][idx] = model.cswm[1].vp[1][idx];

        // up
        model.cswm[0].up[idx][NY-1] = ConvertUPatch2Patch(model, 0, 4, idx, idx, NY-1, 1);
        model.cswm[0].vp[idx][NY-1] = ConvertVPatch2Patch(model, 0, 4, idx, idx, NY-1, 1);
        // model.cswm[0].up[idx][NY-1] = model.cswm[4].up[idx][1];
        // model.cswm[0].vp[idx][NY-1] = model.cswm[4].vp[idx][1];

        // bottom
        model.cswm[0].up[idx][0] = ConvertUPatch2Patch(model, 0, 5, idx, idx, 0, NY-2);
        model.cswm[0].vp[idx][0] = ConvertVPatch2Patch(model, 0, 5, idx, idx, 0, NY-2);
        // model.cswm[0].up[idx][0] = model.cswm[5].up[idx][NY-2];
        // model.cswm[0].vp[idx][0] = model.cswm[5].vp[idx][NY-2];
    }

    // patch 2
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[1].up[0][idx] = ConvertUPatch2Patch(model, 1, 0, 0, NX-2, idx, idx);
        model.cswm[1].vp[0][idx] = ConvertVPatch2Patch(model, 1, 0, 0, NX-2, idx, idx);
        // model.cswm[1].up[0][idx] = model.cswm[0].up[NX-2][idx];
        // model.cswm[1].vp[0][idx] = model.cswm[0].vp[NX-2][idx];

        // right
        model.cswm[1].up[NX-1][idx] = ConvertUPatch2Patch(model, 1, 2, NX-1, 1, idx, idx);
        model.cswm[1].vp[NX-1][idx] = ConvertVPatch2Patch(model, 1, 2, NX-1, 1, idx, idx);
        // model.cswm[1].up[NX-1][idx] = model.cswm[2].up[1][idx];
        // model.cswm[1].vp[NX-1][idx] = model.cswm[2].vp[1][idx];

        // up
        // model.cswm[1].up[idx][NY-1] = ConvertUPatch2Patch(model, 1, 4, idx, NX-2, NY-1, idx);
        // model.cswm[1].vp[idx][NY-1] = ConvertVPatch2Patch(model, 1, 4, idx, NX-2, NY-1, idx);
        model.cswm[1].up[idx][NY-1] = ConvertBV2AUPatch2Patch(model, 1, 4, idx, NX-2, NY-1, idx);
        model.cswm[1].vp[idx][NY-1] = ConvertBU2AVPatch2Patch(model, 1, 4, idx, NX-2, NY-1, idx);
        // model.cswm[1].up[idx][NY-1] = model.cswm[4].vp[NX-2][idx];
        // model.cswm[1].vp[idx][NY-1] = -model.cswm[4].up[NX-2][idx];

        // bottom
        // model.cswm[1].up[idx][0] = ConvertUPatch2Patch(model, 1, 5, idx, NX-2, 0, (NY-1)-idx);
        // model.cswm[1].vp[idx][0] = ConvertVPatch2Patch(model, 1, 5, idx, NX-2, 0, (NY-1)-idx);
        model.cswm[1].up[idx][0] = ConvertBV2AUPatch2Patch(model, 1, 5, idx, NX-2, 0, (NY-1)-idx);
        model.cswm[1].vp[idx][0] = ConvertBU2AVPatch2Patch(model, 1, 5, idx, NX-2, 0, (NY-1)-idx);
        // model.cswm[1].up[idx][0] = -model.cswm[5].vp[NX-2][(NY-1)-idx];
        // model.cswm[1].vp[idx][0] = model.cswm[5].up[NX-2][(NY-1)-idx];
    }

    // patch 3
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[2].up[0][idx] = ConvertUPatch2Patch(model, 2, 1, 0, NX-2, idx, idx);
        model.cswm[2].vp[0][idx] = ConvertVPatch2Patch(model, 2, 1, 0, NX-2, idx, idx);
        // model.cswm[2].up[0][idx] = model.cswm[1].up[NX-2][idx];
        // model.cswm[2].vp[0][idx] = model.cswm[1].vp[NX-2][idx];

        // right
        model.cswm[2].up[NX-1][idx] = ConvertUPatch2Patch(model, 2, 3, NX-1, 1, idx, idx);
        model.cswm[2].vp[NX-1][idx] = ConvertVPatch2Patch(model, 2, 3, NX-1, 1, idx, idx);
        // model.cswm[2].up[NX-1][idx] = model.cswm[3].up[1][idx];
        // model.cswm[2].vp[NX-1][idx] = model.cswm[3].vp[1][idx];

        // up
        // model.cswm[2].up[idx][NY-1] = ConvertUPatch2Patch(model, 2, 4, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[2].vp[idx][NY-1] = ConvertVPatch2Patch(model, 2, 4, idx, (NX-1)-idx, NY-1, NY-2);
        model.cswm[2].up[idx][NY-1] = ConvertUPatch2Patch(model, 2, 4, idx, (NX-1)-idx, NY-1, NY-2);
        model.cswm[2].vp[idx][NY-1] = ConvertVPatch2Patch(model, 2, 4, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[2].up[idx][NY-1] = -model.cswm[4].up[(NX-1)-idx][NY-2];
        // model.cswm[2].vp[idx][NY-1] = -model.cswm[4].vp[(NX-1)-idx][NY-2];

        // bottom
        // model.cswm[2].up[idx][0] = ConvertUPatch2Patch(model, 2, 5, idx, (NX-1)-idx, 0, 1);
        // model.cswm[2].vp[idx][0] = ConvertVPatch2Patch(model, 2, 5, idx, (NX-1)-idx, 0, 1);
        model.cswm[2].up[idx][0] = ConvertUPatch2Patch(model, 2, 5, idx, (NX-1)-idx, 0, 1);
        model.cswm[2].vp[idx][0] = ConvertVPatch2Patch(model, 2, 5, idx, (NX-1)-idx, 0, 1);
        // model.cswm[2].up[idx][0] = -model.cswm[5].up[(NX-1)-idx][1];
        // model.cswm[2].vp[idx][0] = -model.cswm[5].vp[(NX-1)-idx][1];
    }

    // patch 4
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        model.cswm[3].up[0][idx] = ConvertUPatch2Patch(model, 3, 2, 0, NX-2, idx, idx);
        model.cswm[3].vp[0][idx] = ConvertVPatch2Patch(model, 3, 2, 0, NX-2, idx, idx);
        // model.cswm[3].up[0][idx] = model.cswm[2].up[NX-2][idx];
        // model.cswm[3].vp[0][idx] = model.cswm[2].vp[NX-2][idx];

        // right
        model.cswm[3].up[NX-1][idx] = ConvertUPatch2Patch(model, 3, 0, NX-1, 1, idx, idx);
        model.cswm[3].vp[NX-1][idx] = ConvertVPatch2Patch(model, 3, 0, NX-1, 1, idx, idx);
        // model.cswm[3].up[NX-1][idx] = model.cswm[0].up[1][idx];
        // model.cswm[3].vp[NX-1][idx] = model.cswm[0].vp[1][idx];

        // up
        // model.cswm[3].up[idx][NY-1] = ConvertUPatch2Patch(model, 3, 4, idx, 1, NY-1, (NY-1)-idx);
        // model.cswm[3].vp[idx][NY-1] = ConvertVPatch2Patch(model, 3, 4, idx, 1, NY-1, (NY-1)-idx);
        model.cswm[3].up[idx][NY-1] = ConvertBV2AUPatch2Patch(model, 3, 4, idx, 1, NY-1, (NY-1)-idx);
        model.cswm[3].vp[idx][NY-1] = ConvertBU2AVPatch2Patch(model, 3, 4, idx, 1, NY-1, (NY-1)-idx);
        // model.cswm[3].up[idx][NY-1] = -model.cswm[4].vp[1][(NY-1)-idx];
        // model.cswm[3].vp[idx][NY-1] = model.cswm[4].up[1][(NY-1)-idx];

        // bottom
        // model.cswm[3].up[idx][0] = ConvertUPatch2Patch(model, 3, 5, idx, 1, 0, idx);
        // model.cswm[3].vp[idx][0] = ConvertVPatch2Patch(model, 3, 5, idx, 1, 0, idx);
        model.cswm[3].up[idx][0] = ConvertBV2AUPatch2Patch(model, 3, 5, idx, 1, 0, idx);
        model.cswm[3].vp[idx][0] = ConvertBU2AVPatch2Patch(model, 3, 5, idx, 1, 0, idx);
        // model.cswm[3].up[idx][0] = model.cswm[5].vp[1][idx];
        // model.cswm[3].vp[idx][0] = -model.cswm[5].up[1][idx];
    }

    // patch 5
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        // model.cswm[4].up[0][idx] = ConvertUPatch2Patch(model, 4, 3, 0, (NX-1)-idx, idx, NY-2);
        // model.cswm[4].vp[0][idx] = ConvertVPatch2Patch(model, 4, 3, 0, (NX-1)-idx, idx, NY-2);
        model.cswm[4].up[0][idx] = ConvertBV2AUPatch2Patch(model, 4, 3, 0, (NX-1)-idx, idx, NY-2);
        model.cswm[4].vp[0][idx] = ConvertBU2AVPatch2Patch(model, 4, 3, 0, (NX-1)-idx, idx, NY-2);
        // model.cswm[4].up[0][idx] = model.cswm[3].vp[(NX-1)-idx][NY-2];
        // model.cswm[4].vp[0][idx] = -model.cswm[3].up[(NX-1)-idx][NY-2];

        // right
        // model.cswm[4].up[NX-1][idx] = ConvertUPatch2Patch(model, 4, 1, NX-1, idx, idx, NY-2);
        // model.cswm[4].vp[NX-1][idx] = ConvertVPatch2Patch(model, 4, 1, NX-1, idx, idx, NY-2);
        model.cswm[4].up[NX-1][idx] = ConvertBV2AUPatch2Patch(model, 4, 1, NX-1, idx, idx, NY-2);
        model.cswm[4].vp[NX-1][idx] = ConvertBU2AVPatch2Patch(model, 4, 1, NX-1, idx, idx, NY-2);
        // model.cswm[4].up[NX-1][idx] = -model.cswm[1].vp[idx][NY-2];
        // model.cswm[4].vp[NX-1][idx] = model.cswm[1].up[idx][NY-2];

        // up
        // model.cswm[4].up[idx][NY-1] = ConvertUPatch2Patch(model, 4, 2, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[4].vp[idx][NY-1] = ConvertVPatch2Patch(model, 4, 2, idx, (NX-1)-idx, NY-1, NY-2);
        model.cswm[4].up[idx][NY-1] = ConvertUPatch2Patch(model, 4, 2, idx, (NX-1)-idx, NY-1, NY-2);
        model.cswm[4].vp[idx][NY-1] = ConvertVPatch2Patch(model, 4, 2, idx, (NX-1)-idx, NY-1, NY-2);
        // model.cswm[4].up[idx][NY-1] = -model.cswm[2].up[(NX-1)-idx][NY-2];
        // model.cswm[4].vp[idx][NY-1] = -model.cswm[2].vp[(NX-1)-idx][NY-2];

        // bottom
        model.cswm[4].up[idx][0] = ConvertUPatch2Patch(model, 4, 0, idx, idx, 0, NY-2);
        model.cswm[4].vp[idx][0] = ConvertVPatch2Patch(model, 4, 0, idx, idx, 0, NY-2);
        // model.cswm[4].up[idx][0] = model.cswm[0].up[idx][NY-2];
        // model.cswm[4].vp[idx][0] = model.cswm[0].vp[idx][NY-2];
    }

    // patch 6
    for (int idx = 1; idx < NX-1; idx++) {
        // left
        // model.cswm[5].up[0][idx] = ConvertUPatch2Patch(model, 5, 3, 0, idx, idx, 1);
        // model.cswm[5].vp[0][idx] = ConvertVPatch2Patch(model, 5, 3, 0, idx, idx, 1);
        model.cswm[5].up[0][idx] = ConvertBV2AUPatch2Patch(model, 5, 3, 0, idx, idx, 1);
        model.cswm[5].vp[0][idx] = ConvertBU2AVPatch2Patch(model, 5, 3, 0, idx, idx, 1);
        // model.cswm[5].up[0][idx] = -model.cswm[3].vp[idx][1];
        // model.cswm[5].vp[0][idx] = model.cswm[3].up[idx][1];

        // right
        // model.cswm[5].up[NX-1][idx] = ConvertUPatch2Patch(model, 5, 1, NX-1, (NX-1)-idx, idx, 1);
        // model.cswm[5].vp[NX-1][idx] = ConvertVPatch2Patch(model, 5, 1, NX-1, (NX-1)-idx, idx, 1);
        model.cswm[5].up[NX-1][idx] = ConvertBV2AUPatch2Patch(model, 5, 1, NX-1, (NX-1)-idx, idx, 1);
        model.cswm[5].vp[NX-1][idx] = ConvertBU2AVPatch2Patch(model, 5, 1, NX-1, (NX-1)-idx, idx, 1);
        // model.cswm[5].up[NX-1][idx] = model.cswm[1].vp[(NX-1)-idx][1];
        // model.cswm[5].vp[NX-1][idx] = -model.cswm[1].up[(NX-1)-idx][1];

        // up
        model.cswm[5].up[idx][NY-1] = ConvertUPatch2Patch(model, 5, 0, idx, idx, NY-1, 1);
        model.cswm[5].vp[idx][NY-1] = ConvertVPatch2Patch(model, 5, 0, idx, idx, NY-1, 1);
        // model.cswm[5].up[idx][NY-1] = model.cswm[0].up[idx][1];
        // model.cswm[5].vp[idx][NY-1] = model.cswm[0].vp[idx][1];

        // bottom
        // model.cswm[5].up[idx][0] = ConvertUPatch2Patch(model, 5, 2, idx, (NX-1)-idx, 0, 1);
        // model.cswm[5].vp[idx][0] = ConvertVPatch2Patch(model, 5, 2, idx, (NX-1)-idx, 0, 1);
        model.cswm[5].up[idx][0] = ConvertUPatch2Patch(model, 5, 2, idx, (NX-1)-idx, 0, 1);
        model.cswm[5].vp[idx][0] = ConvertVPatch2Patch(model, 5, 2, idx, (NX-1)-idx, 0, 1);
        // model.cswm[5].up[idx][0] = -model.cswm[2].up[(NX-1)-idx][1];
        // model.cswm[5].vp[idx][0] = -model.cswm[2].vp[(NX-1)-idx][1];
    }
    return;
}

void CSWM::BoundaryWindInterpolation(CSWM &model) {
    return;
}

double CSWM::ConvertUPatch2Sphere(CSWM &model, int p, int i, int j) {
    return (model.gUpper_u[i][j][0] * model.cswm[p].A_u[i][j][0] + model.gUpper_u[i][j][1] * model.cswm[p].A_u[i][j][2]) * model.cswm[p].u[i][j] + 
           (model.gUpper_u[i][j][0] * model.cswm[p].A_u[i][j][1] + model.gUpper_u[i][j][1] * model.cswm[p].A_u[i][j][3]) * model.cswm[p].v[i][j];
}

double CSWM::ConvertVPatch2Sphere(CSWM &model, int p, int i, int j) {
    return (model.gUpper_v[i][j][2] * model.cswm[p].A_v[i][j][0] + model.gUpper_v[i][j][3] * model.cswm[p].A_v[i][j][2]) * model.cswm[p].u[i][j] + 
           (model.gUpper_v[i][j][2] * model.cswm[p].A_v[i][j][1] + model.gUpper_v[i][j][3] * model.cswm[p].A_v[i][j][3]) * model.cswm[p].v[i][j];
}

double CSWM::ConvertUSphere2Patch(CSWM &model, int p, int i, int j) {
    return (model.cswm[p].AInverse_u[i][j][0] * model.cswm[p].u[i][j] + 
            model.cswm[p].AInverse_u[i][j][1] * 0.25 * (model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j] + model.cswm[p].v[i-1][j+1] + model.cswm[p].v[i-1][j]));
}

double CSWM::ConvertVSphere2Patch(CSWM &model, int p, int i, int j) {
    return (model.cswm[p].AInverse_v[i][j][0] * 0.25 * (model.cswm[p].u[i+1][j] + model.cswm[p].u[i][j] + model.cswm[p].u[i+1][j-1] + model.cswm[p].u[i][j-1]) + 
            model.cswm[p].AInverse_v[i][j][1] * model.cswm[p].v[i][j+1]);
}

double CSWM::ConvertUPatch2Patch(CSWM &model, int p1, int p2, int i1, int i2, int j1, int j2) {
    double C[2][2];
    double D[2][2];
    double final[2][2];
    double gLowerA[2][2], gUpperB[2][2], AInverseA[2][2], AB[2][2];
    int count = 0;

    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            C[i][j] = D[i][j] = final[i][j] = gLowerA[i][j] = gUpperB[i][j] = AInverseA[i][j] = AB[i][j] = 0.;
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            gLowerA[i][j] = model.gLower_u[i1][j1][count];
            gUpperB[i][j] = model.gUpper_u[i2][j2][count];
            AInverseA[i][j] = model.cswm[p1].AInverse_u[i1][j1][count];
            AB[i][j] = model.cswm[p2].A_u[i2][j2][count];
            count++;
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {   
            for (int k = 0; k < 2; k++) {
                C[i][j] += gLowerA[i][k] * AInverseA[k][j];
                D[i][j] += AB[i][k] * gUpperB[k][j];
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {  
            for (int k = 0; k < 2; k++) {
                final[i][j] += C[i][k] * D[k][j];
            }
        }
    }
    return final[0][0] * model.cswm[p2].up[i2][j2] + final[0][1] * 0.25*(model.cswm[p2].vp[i2][j2+1] + model.cswm[p2].vp[i2][j2] + model.cswm[p2].vp[i2-1][j2+1] + model.cswm[p2].vp[i2-1][j2]);
}

double CSWM::ConvertVPatch2Patch(CSWM &model, int p1, int p2, int i1, int i2, int j1, int j2) {
    double C[2][2];
    double D[2][2];
    double final[2][2];
    double gLowerA[2][2], gUpperB[2][2], AInverseA[2][2], AB[2][2];
    int count = 0;

    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            C[i][j] = D[i][j] = final[i][j] = gLowerA[i][j] = gUpperB[i][j] = AInverseA[i][j] = AB[i][j] = 0.;
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            gLowerA[i][j] = model.gLower_v[i1][j1][count];
            gUpperB[i][j] = model.gUpper_v[i2][j2][count];
            AInverseA[i][j] = model.cswm[p1].AInverse_v[i1][j1][count];
            AB[i][j] = model.cswm[p2].A_v[i2][j2][count];
            count++;
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {   
            for (int k = 0; k < 2; k++) {
                C[i][j] += gLowerA[i][k] * AInverseA[k][j];
                D[i][j] += AB[i][k] * gUpperB[k][j];
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {  
            for (int k = 0; k < 2; k++) {
                final[i][j] += C[i][k] * D[k][j];
            }
        }
    }
    return final[1][0] * 0.25*(model.cswm[p2].up[i2+1][j2] + model.cswm[p2].up[i2][j2] + model.cswm[p2].up[i2+1][j2-1] + model.cswm[p2].up[i2][j2-1]) + final[1][1] * model.cswm[p2].vp[i2][j2];
}

double CSWM::ConvertBV2AUPatch2Patch(CSWM &model, int p1, int p2, int i1, int i2, int j1, int j2) {
    double C[2][2];
    double D[2][2];
    double final[2][2];
    double gLowerA[2][2], gUpperB[2][2], AInverseA[2][2], AB[2][2];
    int count = 0;

    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            C[i][j] = D[i][j] = final[i][j] = gLowerA[i][j] = gUpperB[i][j] = AInverseA[i][j] = AB[i][j] = 0.;
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            gLowerA[i][j] = model.gLower_u[i1][j1][count];
            gUpperB[i][j] = model.gUpper_v[i2][j2][count];
            AInverseA[i][j] = model.cswm[p1].AInverse_u[i1][j1][count];
            AB[i][j] = model.cswm[p2].A_v[i2][j2][count];
            count++;
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {   
            for (int k = 0; k < 2; k++) {
                C[i][j] += gLowerA[i][k] * AInverseA[k][j];
                D[i][j] += AB[i][k] * gUpperB[k][j];
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {  
            for (int k = 0; k < 2; k++) {
                final[i][j] += C[i][k] * D[k][j];
            }
        }
    }
    return final[0][0] * model.cswm[p2].up[i2][j2] + final[0][1] * 0.25*(model.cswm[p2].vp[i2][j2+1] + model.cswm[p2].vp[i2][j2] + model.cswm[p2].vp[i2-1][j2+1] + model.cswm[p2].vp[i2-1][j2]);
}

double CSWM::ConvertBU2AVPatch2Patch(CSWM &model, int p1, int p2, int i1, int i2, int j1, int j2) {
    double C[2][2];
    double D[2][2];
    double final[2][2];
    double gLowerA[2][2], gUpperB[2][2], AInverseA[2][2], AB[2][2];
    int count = 0;

    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            C[i][j] = D[i][j] = final[i][j] = gLowerA[i][j] = gUpperB[i][j] = AInverseA[i][j] = AB[i][j] = 0.;
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            gLowerA[i][j] = model.gLower_v[i1][j1][count];
            gUpperB[i][j] = model.gUpper_u[i2][j2][count];
            AInverseA[i][j] = model.cswm[p1].AInverse_v[i1][j1][count];
            AB[i][j] = model.cswm[p2].A_u[i2][j2][count];
            count++;
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {   
            for (int k = 0; k < 2; k++) {
                C[i][j] += gLowerA[i][k] * AInverseA[k][j];
                D[i][j] += AB[i][k] * gUpperB[k][j];
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {  
            for (int k = 0; k < 2; k++) {
                final[i][j] += C[i][k] * D[k][j];
            }
        }
    }
    return final[1][0] * 0.25*(model.cswm[p2].up[i2+1][j2] + model.cswm[p2].up[i2][j2] + model.cswm[p2].up[i2+1][j2-1] + model.cswm[p2].up[i2][j2-1]) + final[1][1] * model.cswm[p2].vp[i2][j2];
}


void CSWM::ExtrapolationBoundary(CSWM &model) {
    for (int p = 0; p < 6; p++) {
        for (int idx = 0; idx < NX; idx++) {
            // left
            model.cswm[p].up[0][idx] = ExtrapolateRight(model.cswm[p].x_u[1][idx], model.cswm[p].x_u[2][idx], model.cswm[p].up[1][idx], model.cswm[p].up[2][idx], model.cswm[p].x_u[0][idx]);
            model.cswm[p].vp[0][idx] = ExtrapolateRight(model.cswm[p].x_v[1][idx], model.cswm[p].x_v[2][idx], model.cswm[p].vp[1][idx], model.cswm[p].vp[2][idx], model.cswm[p].x_v[0][idx]);

            // right
            model.cswm[p].up[NX-1][idx] = ExtrapolateRight(model.cswm[p].x_u[NX-3][idx], model.cswm[p].x_u[NX-2][idx], model.cswm[p].up[NX-3][idx], model.cswm[p].up[NX-2][idx], model.cswm[p].x_u[NX-1][idx]);
            model.cswm[p].vp[NX-1][idx] = ExtrapolateRight(model.cswm[p].x_v[NX-3][idx], model.cswm[p].x_v[NX-2][idx], model.cswm[p].vp[NX-3][idx], model.cswm[p].vp[NX-2][idx], model.cswm[p].x_v[NX-1][idx]);


            // up
            model.cswm[p].up[idx][NX-1] = ExtrapolateRight(model.cswm[p].y_u[idx][NX-3], model.cswm[p].y_u[idx][NX-2], model.cswm[p].up[idx][NX-3], model.cswm[p].up[idx][NX-2], model.cswm[p].y_u[idx][NX-1]);
            model.cswm[p].vp[idx][NX-1] = ExtrapolateRight(model.cswm[p].y_v[idx][NX-3], model.cswm[p].y_v[idx][NX-2], model.cswm[p].vp[idx][NX-3], model.cswm[p].vp[idx][NX-2], model.cswm[p].y_v[idx][NX-1]);

            // bottom
            model.cswm[p].up[idx][0] = ExtrapolateRight(model.cswm[p].y_u[idx][1], model.cswm[p].y_u[idx][2], model.cswm[p].up[idx][1], model.cswm[p].up[idx][2], model.cswm[p].y_u[idx][0]);
            model.cswm[p].vp[idx][0] = ExtrapolateRight(model.cswm[p].y_v[idx][1], model.cswm[p].y_v[idx][2], model.cswm[p].vp[idx][1], model.cswm[p].vp[idx][2], model.cswm[p].y_v[idx][0]);
        }
    }
}

double CSWM::ExtrapolateLeft(double x1, double x2, double y1, double y2, double x) { 
    return y2 - (y2-y1) * (x2-x) / (x2-x1);
}

double CSWM::ExtrapolateRight(double x1, double x2, double y1, double y2, double x) { 
    return y1 + (y2-y1) * (x-x1) / (x2-x1);
}