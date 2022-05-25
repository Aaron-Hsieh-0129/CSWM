#include "Output.hpp"

void Output::output_h(int n, CSWM &model) {
    std::fstream fouth;
    std::string hname = "../outputs/h/h_" + std::to_string(n) + ".txt";
    fouth.open(hname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                fouth << model.cswm[p].h[i][j] << " ";
            }
        }
    }
    return;
}

void Output::output_u(int n, CSWM &model) {
    std::fstream foutu;
    std::string uname = "../outputs/u/u_" + std::to_string(n) + ".txt";
    foutu.open(uname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                foutu << model.cswm[p].u[i][j] << " ";
            }
        }
    }
    return;
}

void Output::output_v(int n, CSWM &model) {
    std::fstream foutv;
    std::string vname = "../outputs/v/v_" + std::to_string(n) + ".txt";
    foutv.open(vname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                foutv << model.cswm[p].v[i][j] << " ";
            }
        }
    }
    return;
}

void Output::output_parameter(CSWM &model) {
    std::fstream foutlon;
    std::fstream foutlat;
    std::fstream foutx;
    std::fstream fouty;
    foutlon.open("../outputs/lon.txt", std::ios::out);
    foutlat.open("../outputs/lat.txt", std::ios::out);
    foutx.open("../outputs/x.txt", std::ios::out);
    fouty.open("../outputs/y.txt", std::ios::out);

    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                foutlon << model.cswm[p].lon[i][j] << " ";
                foutlat << model.cswm[p].lat[i][j] << " ";
                foutx << model.cswm[p].x[i][j] << " ";
                fouty << model.cswm[p].y[i][j] << " ";
            }
        }
    }
    return;
}