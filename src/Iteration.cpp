#include "Iteration.hpp"

void Iteration::ph_pt(CSWM &model) {
    double psqrtGHU_px = 0, psqrtGHU_py = 0, dx_for_h = 0, dy_for_h = 0;
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                
                /*
                // dx_for_h = 0.5 * (model.cswm[p].x[i+1][j] - model.cswm[p].x[i-1][j]);
                // dy_for_h = 0.5 * (model.cswm[p].y[i][j+1] - model.cswm[p].y[i][j-1]);
                dx_for_h = model.cswm[p].x[i+1][j] - model.cswm[p].x[i][j];
                dy_for_h = model.cswm[p].y[i][j+1] - model.cswm[p].y[i][j];

                psqrtGHU_px = (0.125 / (model.sqrtG[i][j] * dx_for_h)) * 
                              (((model.sqrtG[i+1][j]+model.sqrtG[i][j]) * (model.cswm[p].h[i+1][j]+model.cswm[p].h[i][j]) * 
                                ((model.gUpper[i+1][j][0]+model.gUpper[i][j][0]) * model.cswm[p].u[i+1][j] + 
                                  0.25 * (model.gUpper[i+1][j][1]+model.gUpper[i][j][1]) * (model.cswm[p].v[i+1][j+1]+model.cswm[p].v[i][j+1]+model.cswm[p].v[i+1][j]+model.cswm[p].v[i][j]))) - 
                               ((model.sqrtG[i][j]+model.sqrtG[i-1][j]) * (model.cswm[p].h[i][j]+model.cswm[p].h[i-1][j]) * 
                                ((model.gUpper[i][j][0]+model.gUpper[i-1][j][0]) * model.cswm[p].u[i][j] + 
                                  0.25 * (model.gUpper[i][j][1]+model.gUpper[i-1][j][1]) * (model.cswm[p].v[i][j+1]+model.cswm[p].v[i][j]+model.cswm[p].v[i-1][j+1]+model.cswm[p].v[i-1][j]))));

                psqrtGHU_py = (0.125 / (model.sqrtG[i][j] * dy_for_h)) * 
                              (((model.sqrtG[i][j+1]+model.sqrtG[i][j]) * (model.cswm[p].h[i][j+1]+model.cswm[p].h[i][j]) * 
                                ((model.gUpper[i][j+1][3]+model.gUpper[i][j][3]) * model.cswm[p].v[i][j+1] + 
                                  0.25 * (model.gUpper[i][j+1][2]+model.gUpper[i][j][2]) * (model.cswm[p].u[i+1][j+1]+model.cswm[p].u[i][j+1]+model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]))) - 
                               ((model.sqrtG[i][j]+model.sqrtG[i][j-1]) * (model.cswm[p].h[i][j]+model.cswm[p].h[i][j-1]) * 
                                ((model.gUpper[i][j][3]+model.gUpper[i][j-1][3]) * model.cswm[p].v[i][j] + 
                                  0.25 * (model.gUpper[i][j][2]+model.gUpper[i][j-1][2]) * (model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]+model.cswm[p].u[i+1][j-1]+model.cswm[p].u[i][j-1]))));
                */

                dx_for_h = model.cswm[p].x_u[i+1][j] - model.cswm[p].x_u[i][j];
                dy_for_h = model.cswm[p].y_v[i][j+1] - model.cswm[p].y_v[i][j];

                psqrtGHU_px = (1. / (model.sqrtG[i][j] * dx_for_h)) * 
                              ((model.sqrtG_u[i+1][j] * 0.5*(model.cswm[p].h[i+1][j]+model.cswm[p].h[i][j]) * 
                                (model.gUpper_u[i+1][j][0] * model.cswm[p].u[i+1][j] + 
                                 0.25 * model.gUpper_u[i+1][j][1] * (model.cswm[p].v[i+1][j+1]+model.cswm[p].v[i][j+1]+model.cswm[p].v[i+1][j]+model.cswm[p].v[i][j]))) - 
                               (model.sqrtG_u[i][j] * 0.5*(model.cswm[p].h[i][j]+model.cswm[p].h[i-1][j]) * 
                                (model.gUpper_u[i][j][0] * model.cswm[p].u[i][j] + 
                                 0.25 * model.gUpper_u[i][j][1] * (model.cswm[p].v[i][j+1]+model.cswm[p].v[i][j]+model.cswm[p].v[i-1][j+1]+model.cswm[p].v[i-1][j]))));

                psqrtGHU_py = (1. / (model.sqrtG[i][j] * dy_for_h)) * 
                              ((model.sqrtG_v[i][j+1] * 0.5*(model.cswm[p].h[i][j+1]+model.cswm[p].h[i][j]) * 
                                (model.gUpper_v[i][j+1][3] * model.cswm[p].v[i][j+1] + 
                                 0.25 * model.gUpper_v[i][j+1][2] * (model.cswm[p].u[i+1][j+1]+model.cswm[p].u[i][j+1]+model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]))) - 
                               (model.sqrtG_v[i][j] * 0.5*(model.cswm[p].h[i][j]+model.cswm[p].h[i][j-1]) * 
                                (model.gUpper_v[i][j][3] * model.cswm[p].v[i][j] + 
                                 0.25 * model.gUpper_v[i][j][2] * (model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]+model.cswm[p].u[i+1][j-1]+model.cswm[p].u[i][j-1]))));
                
                model.cswm[p].hp[i][j] = model.cswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);

                // diffusion
                #ifdef DIFFUSION
                    model.cswm[p].hp[i][j] += D2T * Kx * (model.cswm[p].hm[i+1][j] - 2. * model.cswm[p].hm[i][j] + model.cswm[p].hm[i-1][j]) / pow((model.cswm[p].x[i+1][j] - model.cswm[p].x[i][j]), 2) + 
                                        D2T * Ky  * (model.cswm[p].hm[i][j+1] - 2. * model.cswm[p].hm[i][j] + model.cswm[p].hm[i][j-1]) / pow((model.cswm[p].y[i][j+1] - model.cswm[p].y[i][j]), 2);
			    #endif
            }
        }
    }

    // Time filter
    #ifdef TIMEFILTER
        for (int p = 0; p < 6; p++) {
            for (int i = 1; i < NX-1; i++) {
                for (int j = 1; j < NY-1; j++) {
                    model.cswm[p].h[i][j] += TIMETS * (model.cswm[p].hp[i][j] - 2 * model.cswm[p].h[i][j] + model.cswm[p].hm[i][j]);
                }
            }
        }
    #endif
    return;
}

void Iteration::pu_pt(CSWM &model) {
    double dx_for_u = 0, dy_for_u = 0;
    double pgH_px = 0, pU2_px = 0, pUV_px = 0, pV2_px = 0, rotationU = 0;
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                dx_for_u = model.cswm[p].x[i][j] - model.cswm[p].x[i-1][j];
                dy_for_u = 0.5 * (model.cswm[p].y_v[i][j+1] + model.cswm[p].y_v[i-1][j+1]) - 0.5 * (model.cswm[p].y_v[i][j] + model.cswm[p].y_v[i-1][j]); 

                pgH_px = gravity / dx_for_u * (model.cswm[p].h[i][j] - model.cswm[p].h[i-1][j]);

                pU2_px = 0.125 / dx_for_u * 
                         (model.gUpper[i][j][0] * pow((model.cswm[p].u[i+1][j] + model.cswm[p].u[i][j]), 2) -
                          model.gUpper[i-1][j][0] * pow((model.cswm[p].u[i][j] + model.cswm[p].u[i-1][j]), 2));

                pUV_px = 0.25 / dx_for_u * 
                         ((model.gUpper[i][j][1] * (model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]) * (model.cswm[p].v[i][j+1]+model.cswm[p].v[i][j])) - 
                          (model.gUpper[i-1][j][1] * (model.cswm[p].u[i][j]+model.cswm[p].u[i-1][j]) * (model.cswm[p].v[i-1][j+1]+model.cswm[p].v[i-1][j])));

                pV2_px = 0.125 / dx_for_u * 
                         (model.gUpper[i][j][3] * pow((model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j]), 2) -
                          model.gUpper[i-1][j][3] * pow((model.cswm[p].v[i-1][j+1] + model.cswm[p].v[i-1][j]), 2));

                // rotationU = (0.5*(((model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j]) - (model.cswm[p].v[i-1][j+1] + model.cswm[p].v[i-1][j])) / dx_for_u) - 
                //              0.5*((model.cswm[p].u[i][j+1] - model.cswm[p].u[i][j-1]) / dy_for_u) + 2 * model.sqrtG_u[i][j] * omega * sin(model.cswm[p].lat_u[i][j])) * 
                //             (model.gUpper_u[i][j][2] * model.cswm[p].u[i][j] + 
                //                0.25 * model.gUpper_u[i][j][3] * (model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j] + model.cswm[p].v[i-1][j+1] + model.cswm[p].v[i-1][j]));

                rotationU = (0.5*(((model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j]) - (model.cswm[p].v[i-1][j+1] + model.cswm[p].v[i-1][j])) / dx_for_u) - 
                             0.5*((model.cswm[p].u[i][j+1] - model.cswm[p].u[i][j-1]) / dy_for_u) + model.sqrtG_u[i][j] * f) * 
                            (model.gUpper_u[i][j][2] * model.cswm[p].u[i][j] + 
                               0.25 * model.gUpper_u[i][j][3] * (model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j] + model.cswm[p].v[i-1][j+1] + model.cswm[p].v[i-1][j]));

                model.cswm[p].up[i][j] = model.cswm[p].um[i][j] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);
                // model.cswm[p].up[i][j] = model.cswm[p].um[i][j] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px);
            }
        }
    }
}

void Iteration::pv_pt(CSWM &model) {
    double dx_for_v = 0, dy_for_v = 0;
    double pgH_py = 0, pU2_py = 0, pUV_py = 0, pV2_py = 0, rotationV = 0;
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                dx_for_v = 0.5 * (model.cswm[p].x_u[i+1][j] + model.cswm[p].x_u[i+1][j-1]) - 0.5 * (model.cswm[p].x_u[i][j] + model.cswm[p].x_u[i][j-1]);
                dy_for_v = model.cswm[p].y[i][j] - model.cswm[p].y[i][j-1];

                pgH_py = gravity / dy_for_v * (model.cswm[p].h[i][j] - model.cswm[p].h[i][j-1]);

                pU2_py = 0.125 / dy_for_v * (model.gUpper[i][j][0] * pow((model.cswm[p].u[i+1][j] + model.cswm[p].u[i][j]), 2) -
                                             model.gUpper[i][j-1][0] * pow((model.cswm[p].u[i+1][j-1] + model.cswm[p].u[i][j-1]), 2));

                pUV_py = 0.25 / dy_for_v * 
                         ((model.gUpper[i][j][1] * (model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]) * (model.cswm[p].v[i][j+1]+model.cswm[p].v[i][j])) - 
                          (model.gUpper[i][j-1][1] * (model.cswm[p].u[i+1][j-1]+model.cswm[p].u[i][j-1]) * (model.cswm[p].v[i][j]+model.cswm[p].v[i][j-1])));

                pV2_py = 0.125 / dy_for_v * (model.gUpper[i][j][3] * pow((model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j]), 2) -
                                             model.gUpper[i][j-1][3] * pow((model.cswm[p].v[i][j] + model.cswm[p].v[i][j-1]), 2));

                // rotationV = (0.5*((model.cswm[p].v[i+1][j] - model.cswm[p].v[i-1][j]) / dx_for_v) - 
                //              0.5*(((model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]) - (model.cswm[p].u[i+1][j-1]+model.cswm[p].u[i][j-1])) / dy_for_v) + model.sqrtG_v[i][j] * 2*omega*sin(model.cswm[p].lat_v[i][j])) * 
                //             (model.gUpper_v[i][j][1] * model.cswm[p].v[i][j] + 
                //                 0.25 * model.gUpper_v[i][j][0] * (model.cswm[p].u[i+1][j] + model.cswm[p].u[i][j] + model.cswm[p].u[i+1][j-1] + model.cswm[p].u[i][j-1]));

                rotationV = (0.5*((model.cswm[p].v[i+1][j] - model.cswm[p].v[i-1][j]) / dx_for_v) - 
                             0.5*(((model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]) - (model.cswm[p].u[i+1][j-1]+model.cswm[p].u[i][j-1])) / dy_for_v) + model.sqrtG_v[i][j] * f) * 
                            (model.gUpper_v[i][j][1] * model.cswm[p].v[i][j] + 
                                0.25 * model.gUpper_v[i][j][0] * (model.cswm[p].u[i+1][j] + model.cswm[p].u[i][j] + model.cswm[p].u[i+1][j-1] + model.cswm[p].u[i][j-1]));

                model.cswm[p].vp[i][j] = model.cswm[p].vm[i][j] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
                // model.cswm[p].vp[i][j] = model.cswm[p].vm[i][j] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py);
            }
        }
    }
}


void Iteration::Leapfrog(CSWM &model) {
    // for (int p = 0; p < 6; p++) {
    //     Output::outputGrid(p, model);
    // }
    Output::output_parameter(model);
    int n = 0;
    double timenow = 0.;
    double temp = TIMEEND / DT;
    int nmax = (int) temp;
    while (n < nmax) {
        std::cout << n << std::endl;

        // TODO: OUTPUT file
        if (n % OUTPUTINTERVAL == 0) {
            // for (int p = 0; p < 6; p++) {
            //     Output::outputPatch(n, p, model);
            // }
            Output::output_h(n, model);
            Output::output_u(n, model);
            Output::output_v(n, model);
        }

        n++;
        timenow = n * DT;

        // calculate
        ph_pt(model);
        pu_pt(model);
        pv_pt(model);
        model.BoundaryProcess(model);
        // model.ExtrapolationBoundary(model);
        // model.BoundaryTransform(model);

        // next step
        for (int p = 0; p < 6; p++) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    model.cswm[p].hm[i][j] = model.cswm[p].h[i][j];
                    model.cswm[p].h[i][j] = model.cswm[p].hp[i][j];

                    model.cswm[p].um[i][j] = model.cswm[p].u[i][j];
                    model.cswm[p].u[i][j] = model.cswm[p].up[i][j];

                    model.cswm[p].vm[i][j] = model.cswm[p].v[i][j];
                    model.cswm[p].v[i][j] = model.cswm[p].vp[i][j];
                }
            }
        }
    }
    return;
}