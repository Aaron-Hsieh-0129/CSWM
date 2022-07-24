#define WORKDIRECTORY "/Users/Aaron/CSWM/src"

#define FILLVALUE (-999999999999999999999999990999.)
#define gravity (9.80665)
#define radius (6371220.)
#define omega (7.292E-5)
#define f (0.)

#define DX (2)
#define DY (2)
#define xMax (45.)
#define xMin (-45.)
#define yMax (45.)
#define yMin (-45.)
#define NX ((int) (xMax - xMin)/DX+2)
#define NY ((int) (yMax - yMin)/DY+2)
#define Kx (60000.)
#define Ky (60000.)
#define TIMETS (0.06)

#define DT (360.)
#define D2T (2. * DT)
#define TIMEEND (24. * 86400*6)
#define OUTPUTINTERVAL (10)

// #define HorizontalAdvection
// #define VerticalAdvection
// #define Jung
// #define Williamson
// #define ConvergenceRate
#define Geostrophy
// #define SteadyGeostrophy
#define ALPHA0 (0.)

#define DIFFUSION
#define TIMEFILTER