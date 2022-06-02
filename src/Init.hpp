#include "Declare.hpp"

class Init {
    public:
        Init();
        static void Init2d(CSWM &);

    private:
        static double GetH(double, double);
        static double GetU(double, double);
        static double GetV(double);
        static double GetUU(double, double);
        static double GetVV(double);
};