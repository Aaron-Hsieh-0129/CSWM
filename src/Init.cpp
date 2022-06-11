#include "Init.hpp"

void Init::Init2d(CSWM & model) {
    // Init h, u, v
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                // Jung
                // model.cswm[p].hp[i][j] = GetH(model.cswm[p].lon[i][j], model.cswm[p].lat[i][j]);
                // model.cswm[p].up[i][j] = GetU(model.cswm[p].lon_u[i][j], model.cswm[p].lat_u[i][j]);
                // model.cswm[p].vp[i][j] = GetV(model.cswm[p].lon_v[i][j]);

                // model.cswm[p].up[i][j] = model.cswm[p].AInverse_u[i][j][0] * model.cswm[p].up[i][j] + model.cswm[p].AInverse_u[i][j][1] * GetV(model.cswm[p].lon_u[i][j]);
                // model.cswm[p].vp[i][j] = model.cswm[p].AInverse_v[i][j][2] * GetU(model.cswm[p].lon_v[i][j], model.cswm[p].lat_v[i][j]) + model.cswm[p].AInverse_v[i][j][3] * model.cswm[p].vp[i][j];
                
                // Williamson
                model.cswm[p].hp[i][j] = GetHH(model.cswm[p].lon[i][j], model.cswm[p].lat[i][j]);
                model.cswm[p].up[i][j] = GetUU(model.cswm[p].lon_u[i][j], model.cswm[p].lat_u[i][j]);
                model.cswm[p].vp[i][j] = GetVV(model.cswm[p].lon_v[i][j]);  

                model.cswm[p].up[i][j] = model.cswm[p].AInverse_u[i][j][0] * model.cswm[p].up[i][j] + model.cswm[p].AInverse_u[i][j][1] * GetVV(model.cswm[p].lon_u[i][j]);
                model.cswm[p].vp[i][j] = model.cswm[p].AInverse_v[i][j][2] * GetUU(model.cswm[p].lon_v[i][j], model.cswm[p].lat_v[i][j]) + model.cswm[p].AInverse_v[i][j][3] * model.cswm[p].vp[i][j];

                // Jung Verical test
                // if (p == 1 || p == 3) {
                //     model.cswm[p].up[i][j] = 0.;
                //     model.cswm[p].vp[i][j] = 0.;
                // }

                // William Verical test
                // if (p == 0 || p == 2) {
                //     model.cswm[p].up[i][j] = 0.;
                //     model.cswm[p].vp[i][j] = 0.;
                // }

                // William Horizontal
                // if (p == 4 || p == 5) {
                //     model.cswm[p].up[i][j] = 0.;
                //     model.cswm[p].vp[i][j] = 0.;
                // }

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
            }
        }
    }
    model.BoundaryProcess(model);

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

    // for (int p = 0; p < 6; p++) {
    //     for (int idx = 0; idx < NX; idx++) {
    //         model.cswm[p].u[0][idx] = GetU(model.cswm[p].lon[0][idx] - 0.5 * (model.cswm[p].lon[1][idx] - model.cswm[p].lon[0][idx]), model.cswm[p].lat[0][idx]);
    //         model.cswm[p].u[NX-1][idx] = GetU(model.cswm[p].lon[NX-1][idx] + 0.5 * (model.cswm[p].lon[NX-1][idx] - model.cswm[p].lon[NX-2][idx]), model.cswm[p].lat[NX-1][idx]);
    //         model.cswm[p].u[idx][0] = GetU(0.5*(model.cswm[p].lon[idx][0] + model.cswm[p].lon[idx-1][0]), model.cswm[p].lat[idx][0]);
    //         model.cswm[p].u[idx][NY-1] = GetU(0.5*(model.cswm[p].lon[idx][NY-1] + model.cswm[p].lon[idx-1][NY-1]), model.cswm[p].lat[idx][NY-1]);


    //         model.cswm[p].v[0][idx] = GetV(model.cswm[p].lon[0][idx]);
    //         model.cswm[p].v[NX-1][idx] = GetV(model.cswm[p].lon[NX-1][idx]);
    //         model.cswm[p].v[idx][0] = GetV(model.cswm[p].lon[idx][0]);
    //         model.cswm[p].v[idx][NY-1] = GetV(model.cswm[p].lon[idx][NY-1]);

    //         model.cswm[p].u[0][j] = GetUU(model.cswm[p].lon[0][j] - 0.5 * (model.cswm[p].lon[1][j] - model.cswm[p].lon[0][j]), model.cswm[p].lat[0][j]);
    //         model.cswm[p].u[NX-1][j] = GetUU(model.cswm[p].lon[NX-1][j] + 0.5 * (model.cswm[p].lon[NX-1][j] - model.cswm[p].lon[NX-2][j]), model.cswm[p].lat[NX-1][j]);
    //     }
    // }


    // for (int p = 0; p < 6; p++) {
    //     for (int i = 0; i < NX; i++) {
    //         for (int j = 0; j < NY; j++) {
    //             model.cswm[p].hm[i][j] = model.cswm[p].h[i][j]; 
    //             model.cswm[p].um[i][j] = model.cswm[p].u[i][j]; 
    //             model.cswm[p].vm[i][j] = model.cswm[p].v[i][j];   

    //             model.cswm[p].hp[i][j] = model.cswm[p].h[i][j]; 
    //             model.cswm[p].up[i][j] = model.cswm[p].u[i][j]; 
    //             model.cswm[p].vp[i][j] = model.cswm[p].v[i][j];  
    //         }
    //     }
    // }

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

    // for (int idx = 0; idx < NX; idx++) {
    //     model.cswm[2].v[idx][0] = model.cswm[5].v[(NX-1)-idx][1];
    //     model.cswm[2].u[idx][0] = model.cswm[5].u[(NX-1)-idx][1];

    //     model.cswm[5].v[idx][0] = model.cswm[2].v[(NX-1)-idx][1];
    //     model.cswm[5].u[idx][0] = model.cswm[2].u[(NX-1)-idx][1];
    // }
    // for (int idx = 0; idx < NX; idx++) {
    //     std::cout << model.cswm[2].u[idx][0] << " " << model.cswm[5].u[NX-1-idx][1] << std::endl;
    //     std::cout << model.cswm[5].lon[idx][0] << " " << model.cswm[5].lon[idx][0] << std::endl;
    // }
    // std::cout << model.cswm[2].um[NX/2][1] << " " << model.cswm[2].u[NX/2][1] << " " << model.cswm[2].up[NX/2][1] << std::endl;
    // std::cout << model.cswm[2].vm[NX/2][1] << " " << model.cswm[2].v[NX/2][1] << " " << model.cswm[2].vp[NX/2][1] << std::endl;


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
    double lonC = 0., latC = 0.;
    double rd = radius * acos(sin(latC) * sin(lat) + cos(latC) * cos(lat) * cos(lon-lonC));
    double r0 = radius / 3.;
    if (rd < r0) return h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return 0.;
}

double Init::GetU(double lon, double lat) {
    double u0 = 2 * M_PI * radius / (12. * 86400);
    // double alpha0 = 0.;
    // double alpha0 = M_PI / 2.;
    double alpha0 = M_PI / 4.;
    double u = u0 * (cos(alpha0) * cos(lat) + sin(alpha0) * sin(lon) * sin(lat));
    return u;
}

double Init::GetV(double lon) {
    double u0 = 2 * M_PI * radius / (12. * 86400);
    // double alpha0 = 0.;
    // double alpha0 = M_PI / 2.;
    double alpha0 = M_PI / 4.;
    double v = u0 * sin(alpha0) * cos(lon);
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
    // double alpha0 = 0.;
    // double alpha0 = M_PI / 2.;
    double alpha0 = M_PI / 4.;
    double u0 = 2 * M_PI * radius / (12. * 86400);
    double u = u0 * (cos(alpha0) * cos(lat) + sin(alpha0) * cos(lon) * sin(lat));
    return u;
}

double Init::GetVV(double lon) {
    // double alpha0 = 0.;
    // double alpha0 = M_PI / 2.;
    double alpha0 = M_PI / 4.;
    double u0 = 2 * M_PI * radius / (12. * 86400);
    double v = -u0 * sin(alpha0) * sin(lon);
    return v;
}
