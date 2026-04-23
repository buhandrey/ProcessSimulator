//Убрать сложные методы в конец файла

#ifndef RULKOV_MAP_H
#define RULKOV_MAP_H

#include "mapping.h"

class Rulkov_map : public mapping {
    //https://doi.org/10.1016/j.cnsns.2010.05.007 (The effects of synaptic time delay on motifs of chemically coupled Rulkov model neurons, Igor Franović, Vladimir Miljković, 2011)
    //Dynamics of Discrete Memristor-Based Rulkov Neuron (part III), Li Jun Liu, Ying Hua Qin, and Du Qu Wei
    //https://en.wikipedia.org/wiki/Rulkov_map
protected:
    double par_a, par_m, par_s, x0, y0;
    bool par_a_def, par_m_def, par_s_def, par_buffersize_def, x0_def, y0_def;
public:
    Rulkov_map () {
        par_a_def = false;
        par_m_def = false;
        par_s_def = false;
        par_buffersize_def = false;
        x0_def = false;
        y0_def = false;
        modelname = "Rulkov_map";
    }
    ~Rulkov_map() = default;
    std::string get_basename () {
        std::stringstream name_temp;
        name_temp << "_" << modelname << "_x0_" << x0 << "_y0_" << y0 << "_a_" << par_a << "_m_" << par_m << "_s_" << par_s;
        return name_temp.str();
    }
    void set_par_a (double set_par) {
        par_a = set_par;
        par_a_def = true;
        try_init ();
    }
    void set_par_m (double set_par) {
        par_m = set_par;
        par_m_def = true;
        try_init ();
    }
    void set_par_s (double set_par) {
        par_s = set_par;
        par_s_def = true;
        try_init ();
    }
    void set_buffer_size (long set_size) {
        if (calc_lyap)
            resize_and_init (set_size, 2, 2, 1);
        else
            resize_and_init (set_size, 2, 0, 1);
        par_buffersize_def = true;
        try_init ();
    }
    void set_x0 (double set_var0) {
        if (!par_buffersize_def)
            throw std::invalid_argument("First define buffer size (from Rulkov_map.set_x0)\n");
        x0 = set_var0;
        vars.initbuffer(0, x0);
        x0_def = true;
        try_init ();
    }
    void set_y0 (double set_var0) {
        if (!par_buffersize_def)
            throw std::invalid_argument("First define buffer size (from Rulkov_map.set_y0)\n");
        y0 = set_var0;
        vars.initbuffer(1, y0);
        y0_def = true;
        try_init ();
    }
    void try_init () {
        if ( par_a_def && par_m_def && par_s_def && par_buffersize_def && x0_def && y0_def) {
            marker = (marker + 1) % val.size();
            discrete_time++;
            calc_amplitude ();
            initialized = true;
        }
    }
    void function (long kin, double portion, long kout) {
        double old_x = vars.get(0);
        double old_y = vars.get(1);
        double oneplusxsq = 1.0 + old_x*old_x;
        all_vars_temp.set(kout, 0, 0, old_y + par_a / oneplusxsq);
        all_vars_temp.set(kout, 0, 1, old_y - par_m * (old_x + par_s));
        for (long i=0;i<vars_tilda.s1();i++) {
            double old_xt = vars_tilda.get(i, 0);
            double old_yt = vars_tilda.get(i, 1);
            all_vars_temp.set(kout, i+1, 0, old_yt - 2.0*par_a*old_x*old_xt/(oneplusxsq*oneplusxsq));
            all_vars_temp.set(kout, i+1, 1, old_yt - par_m*old_xt);
        }
    }
};

#endif