#include <iostream>
#include <cmath>
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
                double x[NX][NY], y[NX][NY];
                double A[NX][NY][4], AInverse[NX][NY][4];
        };
        CSWM();

        patch cswm[6];
        double sqrtG[NX][NY], gamma[NX][NY], gLower[NX][NY][4], gUpper[NX][NY][4];
        static void BoundaryProcess(CSWM &);
        static double Interpolate(double A1, double A2, double y1, double y2, double B);
        static double Interpolate4lat(double A1, double A2, double y1, double y2, double B);
        static void TransformCoordinate(CSWM &);
        static void ConvertPatch2Sphere(CSWM &, int, int, int);
};
