#include "Iteration.hpp"

void Iteration::ph_pt(CSWM &model) {
    double psqrtGHU_px, psqrtGHU_py;
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                psqrtGHU_px = (0.125 / (model.sqrtG[i][j] * (model.cswm[p].x[i+1][j] - model.cswm[p].x[i][j]))) * 
                              (((model.sqrtG[i+1][j]+model.sqrtG[i][j]) * (model.cswm[p].h[i+1][j]+model.cswm[p].h[i][j]) * 
                                ((model.gUpper[i+1][j][0]+model.gUpper[i][j][0]) * model.cswm[p].u[i+1][j] + 
                                  0.25 * (model.gUpper[i+1][j][1]+model.gUpper[i][j][1]) * (model.cswm[p].v[i+1][j+1]+model.cswm[p].v[i][j+1]+model.cswm[p].v[i+1][j]+model.cswm[p].v[i][j]))) - 
                               ((model.sqrtG[i][j]+model.sqrtG[i-1][j]) * (model.cswm[p].h[i][j]+model.cswm[p].h[i-1][j]) * 
                                ((model.gUpper[i][j][0]+model.gUpper[i-1][j][0]) * model.cswm[p].u[i][j] + 
                                  0.25 * (model.gUpper[i][j][1]+model.gUpper[i-1][j][1]) * (model.cswm[p].v[i][j+1]+model.cswm[p].v[i][j]+model.cswm[p].v[i-1][j+1]+model.cswm[p].v[i-1][j]))));

                psqrtGHU_py = (0.125 / (model.sqrtG[i][j] * (model.cswm[p].y[i][j+1] - model.cswm[p].y[i][j]))) * 
                              (((model.sqrtG[i][j+1]+model.sqrtG[i][j]) * (model.cswm[p].h[i][j+1]+model.cswm[p].h[i][j]) * 
                                ((model.gUpper[i][j+1][3]+model.gUpper[i][j][3]) * model.cswm[p].v[i][j+1] + 
                                  0.25 * (model.gUpper[i][j+1][2]+model.gUpper[i][j][2]) * (model.cswm[p].u[i+1][j+1]+model.cswm[p].u[i][j+1]+model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]))) - 
                               ((model.sqrtG[i][j]+model.sqrtG[i][j-1]) * (model.cswm[p].h[i][j]+model.cswm[p].h[i][j-1]) * 
                                ((model.gUpper[i][j][3]+model.gUpper[i][j-1][3]) * model.cswm[p].v[i][j] + 
                                  0.25 * (model.gUpper[i][j][2]+model.gUpper[i][j-1][2]) * (model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]+model.cswm[p].u[i+1][j-1]+model.cswm[p].u[i][j-1]))));

                
                model.cswm[p].hp[i][j] = model.cswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);

                // TODO: diffusion
            }
        }
    }

    // TODO: Time filter
}

void Iteration::pu_pt(CSWM &model) {
    double pgH_px = 0, pU2_px = 0, pUV_px = 0, pV2_px = 0, rotationU = 0;
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                pgH_px = gravity / (model.cswm[p].x[i+1][j] - model.cswm[p].x[i][j]) * (model.cswm[p].h[i][j] - model.cswm[p].h[i-1][j]);

                pU2_px = 0.125 / (model.cswm[p].x[i+1][j] - model.cswm[p].x[i][j]) * 
                         (model.gUpper[i][j][0] * pow((model.cswm[p].u[i+1][j] + model.cswm[p].u[i][j]), 2) -
                          model.gUpper[i-1][j][0] * pow((model.cswm[p].u[i][j] + model.cswm[p].u[i-1][j]), 2));

                pUV_px = 0.25 / (model.cswm[p].x[i+1][j] - model.cswm[p].x[i][j]) * 
                         ((model.gUpper[i][j][1] * (model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]) * (model.cswm[p].v[i][j+1]+model.cswm[p].v[i][j])) - 
                          (model.gUpper[i-1][j][1] * (model.cswm[p].u[i][j]+model.cswm[p].u[i-1][j]) * (model.cswm[p].v[i-1][j+1]+model.cswm[p].v[i-1][j])));

                pV2_px = 0.125 / (model.cswm[p].x[i+1][j] - model.cswm[p].x[i][j]) * 
                         (model.gUpper[i][j][3] * pow((model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j]), 2) -
                          model.gUpper[i-1][j][3] * pow((model.cswm[p].v[i-1][j+1] + model.cswm[p].v[i-1][j]), 2));

                rotationU = 0.25 * ((((model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j]) - (model.cswm[p].v[i-1][j+1] + model.cswm[p].v[i-1][j])) / (model.cswm[p].x[i+1][j] - model.cswm[p].x[i][j])) - 
                                     ((model.cswm[p].u[i][j+1]-model.cswm[p].u[i][j-1]) / (model.cswm[p].y[i][j+1]-model.cswm[p].y[i][j])) + (model.sqrtG[i][j]+model.sqrtG[i-1][j]) * 2*omega*sin(model.cswm[p].lat[i][j])) * 
                             ((model.gUpper[i][j][2] + model.gUpper[i-1][j][2]) * model.cswm[p].u[i][j] + 
                               0.25 * (model.gUpper[i][j][3] + model.gUpper[i-1][j][3]) * 
                                      (model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j] + model.cswm[p].v[i-1][j+1] + model.cswm[p].v[i-1][j]));

                model.cswm[p].up[i][j] = model.cswm[p].um[i][j] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);
            }
        }
    }
}

void Iteration::pv_pt(CSWM &model) {
    double pgH_py = 0, pU2_py = 0, pUV_py = 0, pV2_py = 0, rotationV = 0;
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                pgH_py = gravity / (model.cswm[p].y[i][j+1] - model.cswm[p].y[i][j]) * (model.cswm[p].h[i][j] - model.cswm[p].h[i][j-1]);

                pU2_py = 0.125 / (model.cswm[p].y[i][j+1] - model.cswm[p].y[i][j]) * 
                         (model.gUpper[i][j][0] * pow((model.cswm[p].u[i+1][j] + model.cswm[p].u[i][j]), 2) -
                          model.gUpper[i][j-1][0] * pow((model.cswm[p].u[i+1][j-1] + model.cswm[p].u[i][j-1]), 2));

                pUV_py = 0.25 / (model.cswm[p].y[i][j+1] - model.cswm[p].y[i][j]) * 
                         ((model.gUpper[i][j][1] * (model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]) * (model.cswm[p].v[i][j+1]+model.cswm[p].v[i][j])) - 
                          (model.gUpper[i][j-1][1] * (model.cswm[p].u[i+1][j-1]+model.cswm[p].u[i][j-1]) * (model.cswm[p].v[i][j]+model.cswm[p].v[i][j-1])));

                pV2_py = 0.125 / (model.cswm[p].y[i][j+1] - model.cswm[p].y[i][j]) * 
                         (model.gUpper[i][j][3] * pow((model.cswm[p].v[i][j+1] + model.cswm[p].v[i][j]), 2) -
                          model.gUpper[i][j-1][3] * pow((model.cswm[p].v[i][j] + model.cswm[p].v[i][j-1]), 2));

                rotationV = 0.25 * ((((model.cswm[p].v[i+1][j] - model.cswm[p].v[i-1][j])) / (model.cswm[p].x[i+1][j] - model.cswm[p].x[i][j])) - 
                                    (((model.cswm[p].u[i+1][j]+model.cswm[p].u[i][j]) - (model.cswm[p].u[i+1][j-1]+model.cswm[p].u[i][j-1])) / (model.cswm[p].y[i][j+1]-model.cswm[p].y[i][j])) + 
                                    (model.sqrtG[i][j]+model.sqrtG[i][j-1]) * 2*omega*sin(0.5*(model.cswm[p].lat[i][j]+model.cswm[p].lat[i][j-1]))) * 
                                   ((model.gUpper[i][j][1] + model.gUpper[i][j-1][1]) * model.cswm[p].v[i][j] + 
                                     0.25 * (model.gUpper[i][j][0] + model.gUpper[i][j-1][0]) * 
                                    (model.cswm[p].u[i+1][j] + model.cswm[p].u[i][j] + model.cswm[p].u[i+1][j-1] + model.cswm[p].u[i][j-1]));

                model.cswm[p].vp[i][j] = model.cswm[p].vm[i][j] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
            }
        }
    }
}


void Iteration::Leapfrog(CSWM &model) {
    int n = 0;
    double timenow = 0.;
    double temp = TIMEEND / DT;
    int nmax = (int) temp;
    while (n < nmax) {
        std::cout << n << std::endl;

        // TODO: OUTPUT file
        n++;
        timenow = n * DT;

        // calculate
        ph_pt(model);
        pu_pt(model);
        pv_pt(model);
        model.BoundaryProcess(model);

        // next step
        for (int p = 0; p < 6; p++) {
            for (int i = 0; i < NX-1; i++) {
                for (int j = 0; j < NY-1; j++) {
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