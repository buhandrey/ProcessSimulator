#ifndef TRIMODAL_H
#define TRIMODAL_H

#include "process.h"
#include "levy.h"
#include "ran.h"

class trimodal : public process {
protected:
    levy s1, s2, s3;
    ran s4;
    bool inited1, inited2, inited3;
public:
    trimodal () {
        inited1 = inited2 = inited3 = false;
        s4.set_pars (2, 1.5, 3.0, -221);
        modelname = "trimodal_noise";
    }
    ~trimodal() = default;
    std::string get_basename () {
        std::stringstream name_temp;
        name_temp << s1.get_basename() << s2.get_basename() << s3.get_basename();
        return name_temp.str();
    }
    void step () {
        if (initialized) {
            marker++;
            if (marker == val.size())
                marker = 0;
            discrete_time++;
            s4.step();
            double test_val = s4.var();
            if (test_val < 1.0) {
                s1.step();
                val[marker] = s1.var();
            }
            else {
                if (test_val < 2.0) {
                    s2.step();
                    val[marker] = s2.var();
                }
                else {
                    s3.step();
                    val[marker] = s3.var();
                }
            }
            substep();
        }
        else {
            throw std::invalid_argument("Object trimodal isn't defined before use");
        }
    }
    void set_pars1 (int set_size, double set_par_mu, double set_par_alpha, double set_par_beta, double set_par_sigma, long set_idum) {
        s1.set_pars (2, set_par_mu, set_par_alpha, set_par_beta, set_par_sigma, set_idum);
        inited1 = true;
        if ( (inited2)&&(inited3) ) {
            if (set_size>0) {
                val.resize(set_size);
            }
            else {
                throw std::invalid_argument("Vector size must be positive (from trimodal)\n");
            }
            initialized = true;
        }
    }
    void set_pars2 (int set_size, double set_par_mu, double set_par_alpha, double set_par_beta, double set_par_sigma, long set_idum) {
        s2.set_pars (2, set_par_mu, set_par_alpha, set_par_beta, set_par_sigma, set_idum);
        inited2 = true;
        if ( (inited1)&&(inited3) ) {
            if (set_size>0) {
                val.resize(set_size);
            }
            else {
                throw std::invalid_argument("Vector size must be positive (from trimodal)\n");
            }
            initialized = true;
        }
    }
    void set_pars3 (int set_size, double set_par_mu, double set_par_alpha, double set_par_beta, double set_par_sigma, long set_idum) {
        s3.set_pars (2, set_par_mu, set_par_alpha, set_par_beta, set_par_sigma, set_idum);
        inited3 = true;
        if ( (inited1)&&(inited2) ) {
            if (set_size>0) {
                val.resize(set_size);
            }
            else {
                throw std::invalid_argument("Vector size must be positive (from trimodal)\n");
            }
            initialized = true;
        }
    }
};

#endif