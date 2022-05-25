#include "Declare.hpp"

class Init {
    public:
        Init();
        static void Init2d(CSWM &);

    private:
        double GetH(int, int);
        double GetU(int, int);
        double GetV(int);
};