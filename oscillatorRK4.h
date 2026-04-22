//Убрать сложные методы в конец файла

#ifndef OSCILLATORRK4_H
#define OSCILLATORRK4_H

#include "ds.h"

class oscillatorRK4 : public ds {
public:
    oscillatorRK4() {
        dt = 0.0001;
        modelname = "oscillatorRK4";
    }
    ~oscillatorRK4() = default;
    void step () {
        if (initialized) {
            function (0, 0.0, 0);
            function (0, 0.5*dt, 1);
            function (1, 0.5*dt, 2);
            function (2, dt, 3);
            for (long n=0;n<vars.s1();n++) {
                vars.step(n, all_vars_temp.get(3, 0, n));
                for (long i=0;i<vars_tilda.s1();i++) {
                    vars_tilda.step(i, n, all_vars_temp.get(3, i+1, n));
                }
            }
            marker = (marker + 1) % val.size();
            discrete_time++;
            calc_amplitude ();
            substep();
        }
        else {
            throw std::invalid_argument("Object oscillatorRK4 isn't defined before use\n");
        }
    }
};

#endif