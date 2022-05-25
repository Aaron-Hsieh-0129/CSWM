#define WORKDIRECTORY "/Users/Aaron/CSWM/src"

#define gravity (9.80665)
#define radius (6371220.)
#define omega (7.292E-5)

#define DX (2)
#define DY (2)
#define xMax (45.)
#define xMin (-45.)
#define yMax (45.)
#define yMin (-45.)
#define NX ((int) (xMax - xMin)/DX+2)
#define NY ((int) (yMax - yMin)/DY+2)
#define Kx (1000.)
#define Ky (1000.)
#define TIMETS (0.01)

#define DT (360.)
#define D2T (2. * DT)
#define TIMEEND (12. * 86400)
#define OUTPUTINTERVAL (10)

// #define DIFFUSION
// #define TIMEFILTER