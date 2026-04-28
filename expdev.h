#ifndef EXPDEV_H
#define EXPDEV_H

#include "noise.h"
#include "ran.h"

class expdev : public noise {
protected:
    double par_lambda;
    ran unit_distribution;
    double generator ();
public:
    expdev () { modelname = "exp_noise"; }
    ~expdev() = default;
    std::string get_basename ();
    void set_pars (int set_size, double set_par_lambda, long set_idum);
};

double expdev::generator() {
    double dum;
    do {
        unit_distribution.step();
        dum = unit_distribution.var();
    }
    while (dum==0.0);
    return -log(dum)/par_lambda;
}

std::string expdev::get_basename () {
    std::stringstream name_temp;
    name_temp << "_" << modelname << "_lambda_" << par_lambda << "_idum_" << idum;
    return name_temp.str();
}

void expdev::set_pars (int set_size, double set_par_lambda, long set_idum) {
    if (set_size>0) {
        val.resize(set_size);
    }
    else {
        throw std::invalid_argument("Vector size must be positive (from expdev)\n");
    }
    if (set_par_lambda>.0) {
        par_lambda = set_par_lambda;
    }
    else {
        throw std::invalid_argument("Noise scale must be positive (from expdev)\n");
    }
    if (set_idum<0) {
        idum = set_idum;
        unit_distribution.set_pars (2, 0.5, 1.0, idum);
    }
    else {
        throw std::invalid_argument("Seed idum must be positive (from expdev)\n");
    }
    initialized = true;
}

#endif