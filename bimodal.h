#ifndef BIMODAL_H
#define BIMODAL_H

#include "process.h"
#include "levy.h"
#include "ran.h"

class bimodal : public process {
protected:
    levy s1, s2;
    ran s3;
    bool inited1, inited2;
public:
    bimodal ();
    ~bimodal() = default;
    std::string get_basename ();
    void step ();
    void set_pars1 (int set_size, double set_par_mu, double set_par_alpha, double set_par_beta, double set_par_sigma, long set_idum);
    void set_pars2 (int set_size, double set_par_mu, double set_par_alpha, double set_par_beta, double set_par_sigma, long set_idum);
};

bimodal::bimodal () {
    inited1 = inited2 = false;
    s3.set_pars (2, 0.0, 1.0, -221);
    modelname = "bimodal_noise";
}

std::string bimodal::get_basename () {
    std::stringstream name_temp;
    name_temp << s1.get_basename() << s2.get_basename();
    return name_temp.str();
}

void bimodal::step () {
    if (initialized) {
        marker++;
        if (marker == val.size())
            marker = 0;
        discrete_time++;
        s3.step();
        if (s3.var() < 0.0) {
            s1.step();
            val[marker] = s1.var();
        }
        else {
            s2.step();
            val[marker] = s2.var();
        }
        substep();
    }
    else {
        throw std::invalid_argument("Object bimodal isn't defined before use\n");
    }
}

void bimodal::set_pars1 (int set_size, double set_par_mu, double set_par_alpha, double set_par_beta, double set_par_sigma, long set_idum) {
    s1.set_pars (2, set_par_mu, set_par_alpha, set_par_beta, set_par_sigma, set_idum);
    inited1 = true;
    if (inited2) {
        if (set_size>0) {
            val.resize(set_size);
        }
        else {
            throw std::invalid_argument("Vector size must be positive (from bimodal)\n");
        }
        initialized = true;
    }
}

void bimodal::set_pars2 (int set_size, double set_par_mu, double set_par_alpha, double set_par_beta, double set_par_sigma, long set_idum) {
    s2.set_pars (2, set_par_mu, set_par_alpha, set_par_beta, set_par_sigma, set_idum);
    inited2 = true;
    if (inited1) {
        if (set_size>0) {
            val.resize(set_size);
        }
        else {
            throw std::invalid_argument("Vector size must be positive (from bimodal)\n");
        }
        initialized = true;
    }
}

#endif