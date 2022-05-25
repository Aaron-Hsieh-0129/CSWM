#include <fstream>
#include <iomanip>
#include "Init.hpp"

class Output {
    public:
        static void output_h(int, CSWM &);
        static void output_u(int, CSWM &);
        static void output_v(int, CSWM &);
        static void output_parameter(CSWM &);
};