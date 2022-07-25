#include "Init.hpp"

void Init::Init2d(CSWM & model) {
    // Init h, u, v
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                // Jung
                #ifdef Jung
                    model.cswm[p].hp[i][j] = GetH(model.cswm[p].lon[i][j], model.cswm[p].lat[i][j]);
                    model.cswm[p].up[i][j] = (model.gLower_u[i][j][0] * model.cswm[p].AInverse_u[i][j][0] + model.gLower_u[i][j][1] * model.cswm[p].AInverse_u[i][j][2]) * GetU(model.cswm[p].lon_u[i][j], model.cswm[p].lat_u[i][j]) + 
                                             (model.gLower_u[i][j][0] * model.cswm[p].AInverse_u[i][j][1] + model.gLower_u[i][j][1] * model.cswm[p].AInverse_u[i][j][3]) * GetV(model.cswm[p].lon_u[i][j]);
                    model.cswm[p].vp[i][j] = (model.gLower_v[i][j][2] * model.cswm[p].AInverse_v[i][j][0] + model.gLower_v[i][j][3] * model.cswm[p].AInverse_v[i][j][2]) * GetU(model.cswm[p].lon_v[i][j], model.cswm[p].lat_v[i][j]) + 
                                             (model.gLower_v[i][j][2] * model.cswm[p].AInverse_v[i][j][1] + model.gLower_v[i][j][3] * model.cswm[p].AInverse_v[i][j][3]) * GetV(model.cswm[p].lon_v[i][j]);
                    model.cswm[p].up[i][j] = -model.cswm[p].up[i][j];
                    model.cswm[p].vp[i][j] = -model.cswm[p].vp[i][j];
                #endif

                // Williamson
                #ifdef Williamson
                    model.cswm[p].hp[i][j] = GetHH(model.cswm[p].lon[i][j], model.cswm[p].lat[i][j]);
                    model.cswm[p].up[i][j] = (model.gLower_u[i][j][0] * model.cswm[p].AInverse_u[i][j][0] + model.gLower_u[i][j][1] * model.cswm[p].AInverse_u[i][j][2]) * GetUU(model.cswm[p].lon_u[i][j], model.cswm[p].lat_u[i][j]) + 
                                             (model.gLower_u[i][j][0] * model.cswm[p].AInverse_u[i][j][1] + model.gLower_u[i][j][1] * model.cswm[p].AInverse_u[i][j][3]) * GetVV(model.cswm[p].lon_u[i][j]);
                    model.cswm[p].vp[i][j] = (model.gLower_v[i][j][2] * model.cswm[p].AInverse_v[i][j][0] + model.gLower_v[i][j][3] * model.cswm[p].AInverse_v[i][j][2]) * GetUU(model.cswm[p].lon_v[i][j], model.cswm[p].lat_v[i][j]) + 
                                             (model.gLower_v[i][j][2] * model.cswm[p].AInverse_v[i][j][1] + model.gLower_v[i][j][3] * model.cswm[p].AInverse_v[i][j][3]) * GetVV(model.cswm[p].lon_v[i][j]);
                #endif


                // Vertical
                #ifdef VerticalAdvection
                    model.cswm[p].up[i][j] = 0; 
                    if (p == 0 || p == 4 || p == 5) {
                        model.cswm[p].vp[i][j] = 40;
                    }
                    else if (p == 2) {
                        model.cswm[p].vp[i][j] = -40; 
                    }    
                    else {
                        model.cswm[p].vp[i][j] = 0.; 
                    } 
                #endif

                // Horizontal
                #ifdef HorizontalAdvection
                    model.cswm[p].vp[i][j] = 0; 
                    if (p == 0 || p == 1 || p == 2 || p == 3) {
                        model.cswm[p].up[i][j] = 40;
                    }
                    else {
                        model.cswm[p].up[i][j] = 0; 
                    }    
                #endif

                #ifdef Geostrophy
                    // model.cswm[p].hp[i][j] = GeostrophyH(model.cswm[p].lon[i][j], model.cswm[p].lat[i][j]);
                    model.cswm[p].hp[i][j] = GetH(model.cswm[p].lon[i][j], model.cswm[p].lat[i][j]);
                    model.cswm[p].up[i][j] = 0.;
                    model.cswm[p].vp[i][j] = 0.;
                #endif

                #ifdef ConvergenceRate
                    model.cswm[p].hp[i][j] = ConvergenceRateH(model.cswm[p].lon[i][j], model.cswm[p].lat[i][j]);
                    model.cswm[p].up[i][j] = (model.gLower_u[i][j][0] * model.cswm[p].AInverse_u[i][j][0] + model.gLower_u[i][j][1] * model.cswm[p].AInverse_u[i][j][2]) * GetUU(model.cswm[p].lon_u[i][j], model.cswm[p].lat_u[i][j]) + 
                                             (model.gLower_u[i][j][0] * model.cswm[p].AInverse_u[i][j][1] + model.gLower_u[i][j][1] * model.cswm[p].AInverse_u[i][j][3]) * GetVV(model.cswm[p].lon_u[i][j]);
                    model.cswm[p].vp[i][j] = (model.gLower_v[i][j][2] * model.cswm[p].AInverse_v[i][j][0] + model.gLower_v[i][j][3] * model.cswm[p].AInverse_v[i][j][2]) * GetUU(model.cswm[p].lon_v[i][j], model.cswm[p].lat_v[i][j]) + 
                                             (model.gLower_v[i][j][2] * model.cswm[p].AInverse_v[i][j][1] + model.gLower_v[i][j][3] * model.cswm[p].AInverse_v[i][j][3]) * GetVV(model.cswm[p].lon_v[i][j]);
                #endif

                #ifdef SteadyGeostrophy
                    model.cswm[p].hp[i][j] = GeostrophyH(model.cswm[p].lon[i][j], model.cswm[p].lat[i][j]);
                    model.cswm[p].up[i][j] = (model.gLower_u[i][j][0] * model.cswm[p].AInverse_u[i][j][0] + model.gLower_u[i][j][1] * model.cswm[p].AInverse_u[i][j][2]) * GetUU(model.cswm[p].lon_u[i][j], model.cswm[p].lat_u[i][j]) + 
                                             (model.gLower_u[i][j][0] * model.cswm[p].AInverse_u[i][j][1] + model.gLower_u[i][j][1] * model.cswm[p].AInverse_u[i][j][3]) * GetVV(model.cswm[p].lon_u[i][j]);
                    model.cswm[p].vp[i][j] = (model.gLower_v[i][j][2] * model.cswm[p].AInverse_v[i][j][0] + model.gLower_v[i][j][3] * model.cswm[p].AInverse_v[i][j][2]) * GetUU(model.cswm[p].lon_v[i][j], model.cswm[p].lat_v[i][j]) + 
                                             (model.gLower_v[i][j][2] * model.cswm[p].AInverse_v[i][j][1] + model.gLower_v[i][j][3] * model.cswm[p].AInverse_v[i][j][3]) * GetVV(model.cswm[p].lon_v[i][j]);
                #endif

            }
        }
    }
    // model.cswm[4].vp[1][NX/2] = 10;
    model.BoundaryProcess(model);
    // std::cout << "v: " << model.cswm[3].vp[NX/2][NX-1] << std::endl;
    

    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                model.cswm[p].hm[i][j] = model.cswm[p].hp[i][j]; 
                model.cswm[p].um[i][j] = model.cswm[p].up[i][j]; 
                model.cswm[p].vm[i][j] = model.cswm[p].vp[i][j];   

                model.cswm[p].h[i][j] = model.cswm[p].hp[i][j]; 
                model.cswm[p].u[i][j] = model.cswm[p].up[i][j]; 
                model.cswm[p].v[i][j] = model.cswm[p].vp[i][j];  
            }
        }
    }

    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                if (model.cswm[p].hm[i][j] == FILLVALUE || model.cswm[p].h[i][j] == FILLVALUE || model.cswm[p].hp[i][j] == FILLVALUE) {
                    std::cout << "WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                }
                if (model.cswm[p].um[i][j] == FILLVALUE || model.cswm[p].u[i][j] == FILLVALUE || model.cswm[p].up[i][j] == FILLVALUE) {
                    std::cout << "WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                }
                if (model.cswm[p].vm[i][j] == FILLVALUE || model.cswm[p].v[i][j] == FILLVALUE || model.cswm[p].vp[i][j] == FILLVALUE) {
                    std::cout << "WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                }

                if (model.cswm[p].hm[i][j] != model.cswm[p].h[i][j] || model.cswm[p].hp[i][j] != model.cswm[p].h[i][j]) {
                    std::cout << "WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                }

                if (model.cswm[p].um[i][j] != model.cswm[p].u[i][j] || model.cswm[p].up[i][j] != model.cswm[p].u[i][j]) {
                    std::cout << "WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                }

                if (model.cswm[p].vm[i][j] != model.cswm[p].v[i][j] || model.cswm[p].vp[i][j] != model.cswm[p].v[i][j]) {
                    std::cout << "WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                }
            }
        }
    }    

}


double Init::GetH(double lon, double lat) {
    double h0 = 1000;
    double lonC = 0., latC = 0.;
    double rd = radius * acos(sin(latC) * sin(lat) + cos(latC) * cos(lat) * cos(lon-lonC));
    double r0 = radius / 3.;
    if (rd < r0) return h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return 0.;
}

double Init::GetU(double lon, double lat) {
    double u0 = 2 * M_PI * radius / (12. * 86400);
    double u = u0 * (cos(ALPHA0) * cos(lat) + sin(ALPHA0) * sin(lon) * sin(lat));
    return u;
}

double Init::GetV(double lon) {
    double u0 = 2 * M_PI * radius / (12. * 86400);
    double v = u0 * sin(ALPHA0) * cos(lon);
    return v;
}

double Init::GetHH(double lon, double lat) {
    double h0 = 1000;
    double lonC = 3*M_PI/2., latC = 0.;
    double rd = radius * acos(sin(latC) * sin(lat) + cos(latC) * cos(lat) * cos(lon-lonC));
    double r0 = radius / 3.;
    if (rd < r0) return h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return 0.;
}

double Init::GetUU(double lon, double lat) {
    double u0 = 2 * M_PI * radius / (12. * 86400);
    double u = u0 * (cos(ALPHA0) * cos(lat) + sin(ALPHA0) * cos(lon) * sin(lat));
    return u;
}

double Init::GetVV(double lon) {
    double u0 = 2 * M_PI * radius / (12. * 86400);
    double v = -u0 * sin(ALPHA0) * sin(lon);
    return v;
}

double Init::GeostrophyH(double lon, double lat) {
    double h0 = 2.94E4 / gravity;
    double u0 = 2 * M_PI * radius / (12. * 86400);
    return h0 - (radius * omega * u0 + u0 * u0 / 2) * pow(-cos(lon) * cos(lat) * sin(ALPHA0) + sin(lat) * cos(ALPHA0), 2) / gravity;
}

double Init::ConvergenceRateH(double lon, double lat) {
    double lonP = atan(cos(lat) * sin(lon) / (cos(lat) * cos(lon) * cos(ALPHA0) + sin(lat)  *sin(ALPHA0)));
    double latP = asin(sin(lat) * cos(ALPHA0) - cos(lat) * cos(lon) * sin(ALPHA0));
    int k = 4;
    int p = 4;
    return pow(cos(latP), p) * sin(k * lonP);
}