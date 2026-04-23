#ifndef FHNSIMPLE_H
#define FHNSIMPLE_H
//https://www.sciencedirect.com/science/article/pii/S1571064518301088
//Chimera states in neuronal networks: A review

#include "oscillatorRK4.h"

class FHNsimple : public oscillatorRK4 {
protected:
    double par_a, par_e, x0, y0;
    bool par_a_def, par_e_def, par_buffersize_def, x0_def, y0_def;
public:
    FHNsimple ();
    ~FHNsimple() = default;
    std::string get_basename ();
    void set_par_a (double set_par);
    void set_par_e (double set_par);
    void set_buffer_size (long set_size);
    void set_x0 (double set_var0);
    void set_y0 (double set_var0);
    void try_init ();
    void function (long kin, double portion, long kout);
};

FHNsimple::FHNsimple () {
    par_a_def = false;
    par_e_def = false;
    par_buffersize_def = false;
    x0_def = false;
    y0_def = false;
    modelname = "FHNsimple";
}

std::string FHNsimple::get_basename () {
    std::stringstream name_temp;
    name_temp << "_" << modelname << "_x0_" << x0 << "_y0_" << y0 << "_a_" << par_a << "_e_" << par_e;
    return name_temp.str();
}

void FHNsimple::set_par_a (double set_par) {
    par_a = set_par;
    par_a_def = true;
    try_init ();
}

void FHNsimple::set_par_e (double set_par) {
    par_e = set_par;
    par_e_def = true;
    try_init ();
}

void FHNsimple::set_buffer_size (long set_size) {
    if (calc_lyap)
        resize_and_init (set_size, 2, 2, 5);
    else
        resize_and_init (set_size, 2, 0, 5);
    par_buffersize_def = true;
    try_init ();
}

void FHNsimple::set_x0 (double set_var0) {
    if (!par_buffersize_def)
        throw std::invalid_argument("First define buffer size (from Rulkov_map.set_x0)\n");
    x0 = set_var0;
    vars.initbuffer(0, x0);
    x0_def = true;
    try_init ();
}

void FHNsimple::set_y0 (double set_var0) {
    if (!par_buffersize_def)
        throw std::invalid_argument("First define buffer size (from Rulkov_map.set_y0)\n");
    y0 = set_var0;
    vars.initbuffer(1, y0);
    y0_def = true;
    try_init ();
}

void FHNsimple::try_init () {
    if ( par_a_def && par_e_def && par_buffersize_def && x0_def && y0_def) {
        marker = (marker + 1) % val.size();
        discrete_time++;
        calc_amplitude ();
        initialized = true;
    }
}

void FHNsimple::function (long kin, double portion, long kout) {
    double old_x = vars.get(0) + portion * all_vars_temp.get(kin, 0, 0);
    double old_y = vars.get(1) + portion * all_vars_temp.get(kin, 0, 1);
    double xsq = old_x*old_x;
    all_vars_temp.set(kout, 0, 0, (old_x - old_x*xsq/3.0 - old_y) / par_e);
    all_vars_temp.set(kout, 0, 1, old_x + par_a);
    for (long i=0;i<vars_tilda.s1();i++) {
        double old_xt = vars_tilda.get(i, 0) + portion * all_vars_temp.get(kin, i+1, 0);
        double old_yt = vars_tilda.get(i, 1) + portion * all_vars_temp.get(kin, i+1, 1);
        all_vars_temp.set(kout, i+1, 0, (old_xt * (1.0 - xsq) - old_yt) / par_e);
        all_vars_temp.set(kout, i+1, 1, old_xt);
    }
}

#endif