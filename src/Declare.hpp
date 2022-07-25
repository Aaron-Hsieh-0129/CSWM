#include <iostream>
#include <cmath>
#include <ctime>
#include <thread>
#include "Const.hpp"

class CSWM {
    public:
        // constructer
        class patch {
            public:
                // constructor
                patch();

                double hp[NX][NY], h[NX][NY], hm[NX][NY];
                double up[NX][NY], u[NX][NY], um[NX][NY];
                double vp[NX][NY], v[NX][NY], vm[NX][NY];
                double lon[NX][NY], lat[NX][NY];
                double lon_u[NX][NY], lat_u[NX][NY];
                double lon_v[NX][NY], lat_v[NX][NY];
                double x[NX][NY], y[NX][NY];
                double x_u[NX][NY], y_u[NX][NY];
                double x_v[NX][NY], y_v[NX][NY];

                double A[NX][NY][4], AInverse[NX][NY][4];
                double A_u[NX][NY][4], AInverse_u[NX][NY][4];
                double A_v[NX][NY][4], AInverse_v[NX][NY][4];
                double CF[NX][NY], CF_u[NX][NY], CF_v[NX][NY];
        };
        CSWM();

        patch cswm[6];
        double alpha2D[NX][NY], beta2D[NX][NY];
        double alpha2D_u[NX][NY];
        double beta2D_v[NX][NY];
        double sqrtG[NX][NY], gamma[NX][NY], gLower[NX][NY][4], gUpper[NX][NY][4];
        double sqrtG_u[NX][NY], gamma_u[NX][NY], gLower_u[NX][NY][4], gUpper_u[NX][NY][4];
        double sqrtG_v[NX][NY], gamma_v[NX][NY], gLower_v[NX][NY][4], gUpper_v[NX][NY][4];
        static void BoundaryTransform(CSWM &);
        static void BoundaryProcess(CSWM &);
        static double Interpolate(double, double, double, double, double);
        static double Interpolate4lat(double, double, double, double, double);
        static double ExtrapolateLeft(double, double, double, double, double);
        static double ExtrapolateRight(double, double, double, double, double);
        static double ConvertUPatch2Sphere(CSWM &, int, int, int);
        static double ConvertVPatch2Sphere(CSWM &, int, int, int);
        static double ConvertUSphere2Patch(CSWM &, int, int, int);
        static double ConvertVSphere2Patch(CSWM &, int, int, int);
        static double ConvertUPatch2Patch(CSWM &, int, int, int, int, int, int);
        static double ConvertVPatch2Patch(CSWM &, int, int, int, int, int, int);
        static double ConvertBU2AVPatch2Patch(CSWM &, int, int, int, int, int, int);
        static double ConvertBV2AUPatch2Patch(CSWM &, int, int, int, int, int, int);
        static void ExtrapolationBoundary(CSWM &);
        static void BoundaryWindInterpolation(CSWM &);
};
