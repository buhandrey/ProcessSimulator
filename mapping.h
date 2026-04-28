#ifndef MAPPING_H
#define MAPPING_H

#include "ds.h"

class mapping : public ds {
public:
    mapping() { dt = 1.0; modelname = "mapping"; }
    ~mapping() = default;
    void step ();
};

void mapping::step () {
    if (initialized) {
        function (0, 0.0, 0);
        for (long n=0;n<vars.s1();n++) {
            vars.step(n, all_vars_temp.get(0, 0, n));
            for (long i=0;i<vars_tilda.s1();i++) {
                vars_tilda.step(i, n, all_vars_temp.get(0, i+1, n));
            }
        }
        marker = (marker + 1) % val.size();
        discrete_time++;
        calc_amplitude ();
        substep();
    }
    else {
        throw std::invalid_argument("Object mapping isn't defined before use\n");
    }
}

#endif