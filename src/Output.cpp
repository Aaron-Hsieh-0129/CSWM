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

void Output::outputGrid(int p, CSWM &model) {
    std::fstream foutlon;
    std::fstream foutlat;
    std::fstream foutx;
    std::fstream fouty;

    std::string lonname = "../outputs/patch" + std::to_string(p+1) + "/lon.txt";
    std::string latname = "../outputs/patch" + std::to_string(p+1) + "/lat.txt";
    std::string xname = "../outputs/patch" + std::to_string(p+1) + "/x.txt";
    std::string yname = "../outputs/patch" + std::to_string(p+1) + "/y.txt";

    foutlon.open(lonname, std::ios::out);
    foutlat.open(latname, std::ios::out);
    foutx.open(xname, std::ios::out);
    fouty.open(yname, std::ios::out);
    for (int i = 1; i < NX-1; i++) {
        for (int j = 1; j < NY-1; j++) {
            foutlon << model.cswm[p].lon[i][j] << " ";
            foutlat << model.cswm[p].lat[i][j] << " ";
            foutx << model.cswm[p].x[i][j] << " ";
            fouty << model.cswm[p].y[i][j] << " ";
        }
    }
}


void Output::outputPatch(int n, int p, CSWM &model) {
    std::fstream fouth;
    std::fstream foutu;
    std::fstream foutv;

    std::string hname = "../outputs/patch" + std::to_string(p+1) + "/h/h_" + std::to_string(n) + ".txt";
    std::string uname = "../outputs/patch" + std::to_string(p+1) + "/u/u_" + std::to_string(n) + ".txt";
    std::string vname = "../outputs/patch" + std::to_string(p+1) + "/v/v_" + std::to_string(n) + ".txt";
    fouth.open(hname, std::ios::out);
    foutu.open(uname, std::ios::out);
    foutv.open(vname, std::ios::out);

    for (int i = 1; i < NX-1; i++) {
        for (int j = 1; j < NY-1; j++) {
            fouth << model.cswm[p].h[i][j] << " ";
            foutu << model.cswm[p].u[i][j] << " ";
            foutv << model.cswm[p].v[i][j] << " ";
        }
    }
}