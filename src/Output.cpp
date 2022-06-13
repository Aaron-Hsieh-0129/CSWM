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

    std::fstream foutu_lon_lat;
    std::string u_lon_latname = "../outputs/u_lon_lat/u_lon_lat_" + std::to_string(n) + ".txt";
    foutu_lon_lat.open(u_lon_latname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                foutu << model.cswm[p].u[i][j] << " ";
                foutu_lon_lat << model.ConvertUPatch2Sphere(model, p, i, j) << " ";
            }
        }
    }

    return;
}

void Output::output_v(int n, CSWM &model) {
    std::fstream foutv;
    std::string vname = "../outputs/v/v_" + std::to_string(n) + ".txt";
    foutv.open(vname, std::ios::out);

    std::fstream foutv_lon_lat;
    std::string v_lon_latname = "../outputs/v_lon_lat/v_lon_lat_" + std::to_string(n) + ".txt";
    foutv_lon_lat.open(v_lon_latname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                foutv << model.cswm[p].v[i][j] << " ";
                foutv_lon_lat << model.ConvertVPatch2Sphere(model, p, i, j) << " ";
            }
        }
    }
    return;
}

void Output::output_parameter(CSWM &model) {
    std::fstream foutlon;
    std::fstream foutlat;
    std::fstream foutlon_withghost;
    std::fstream foutlat_withghost;
    std::fstream foutx;
    std::fstream fouty;
    // std::fstream foutA;
    // std::fstream foutAInverse;
    // std::fstream foutA_u;
    // std::fstream foutAInverse_u;
    // std::fstream foutA_v;
    // std::fstream foutAInverse_v;
    foutlon.open("../outputs/lon.txt", std::ios::out);
    foutlat.open("../outputs/lat.txt", std::ios::out);
    foutlon_withghost.open("../outputs/lon_withghost.txt", std::ios::out);
    foutlat_withghost.open("../outputs/lat_withghost.txt", std::ios::out);
    foutx.open("../outputs/x.txt", std::ios::out);
    fouty.open("../outputs/y.txt", std::ios::out);
    // foutA.open("../outputs/ConvertMatrix/A.txt", std::ios::out);
    // foutAInverse.open("../outputs/ConvertMatrix/AInverse.txt", std::ios::out);
    // foutA_u.open("../outputs/ConvertMatrix/A_u.txt", std::ios::out);
    // foutAInverse_u.open("../outputs/ConvertMatrix/AInverse_u.txt", std::ios::out);
    // foutA_v.open("../outputs/ConvertMatrix/A_v.txt", std::ios::out);
    // foutAInverse_v.open("../outputs/ConvertMatrix/AInverse_v.txt", std::ios::out);

    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                foutlon << model.cswm[p].lon[i][j] << " ";
                foutlat << model.cswm[p].lat[i][j] << " ";
                foutx << model.cswm[p].x[i][j] << " ";
                fouty << model.cswm[p].y[i][j] << " ";

                // foutA << model.cswm[p].A[i][j] << " ";
                // foutAInverse << model.cswm[p].AInverse[i][j] << " ";
                // foutA_u << model.cswm[p].A_u[i][j] << " ";
                // foutAInverse_u << model.cswm[p].AInverse_u[i][j] << " ";
                // foutA_v << model.cswm[p].A_v[i][j] << " ";
                // foutAInverse_v << model.cswm[p].AInverse_v[i][j] << " ";
            }
        }
    }

    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) { 
                foutlon_withghost << model.cswm[p].lon[i][j] << " ";
                foutlat_withghost << model.cswm[p].lat[i][j] << " ";
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