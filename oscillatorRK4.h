#ifndef OSCILLATORRK4_H
#define OSCILLATORRK4_H

#include "ds.h"

class oscillatorRK4 : public ds {
protected:
    bool dt_def;
public:
    oscillatorRK4 ();
    ~oscillatorRK4 () = default;
    void step ();
    void set_dt (double set_dt) { dt = set_dt; dt_def = true; }
};

oscillatorRK4::oscillatorRK4 () {
    dt_def = false;
    modelname = "oscillatorRK4";
}

void oscillatorRK4::step () {
    if (initialized) {
        function (0, 0.0, 0);
        function (0, 0.5*dt, 1);
        function (1, 0.5*dt, 2);
        function (2, dt, 3);
        for (long n=0;n<vars.s1();n++) {
            double new_val = vars.get(n)
                           + dt * (all_vars_temp.get(0, 0, n)
                           + 2.0*all_vars_temp.get(1, 0, n)
                           + 2.0*all_vars_temp.get(2, 0, n)
                           + all_vars_temp.get(3, 0, n))/6.0;
            vars.step(n, new_val);
            for (long i=0;i<vars_tilda.s1();i++) {
                new_val = vars_tilda.get(i, n)
                        + dt * (all_vars_temp.get(0, i+1, n)
                        + 2.0*all_vars_temp.get(1, i+1, n)
                        + 2.0*all_vars_temp.get(2, i+1, n)
                        + all_vars_temp.get(3, i+1, n))/6.0;
                vars_tilda.step(i, n, new_val);
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

#endif