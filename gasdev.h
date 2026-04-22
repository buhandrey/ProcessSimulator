#ifndef GASDEV_H
#define GASDEV_H

#include "noise.h"
#include "ran.h"

class gasdev : public noise {
protected:
    double par_mu;
    double par_sigma;
    /*Статические параметры генератора чисел*/
    int iset;
    float gset;
    /*Статические параметры генератора чисел*/
    ran unit_doubled_distribution;
    double generator () {
        float fac, rsq, v1, v2;
        if (idum < 0) iset=0;
        if (iset == 0) {
            do {
                unit_doubled_distribution.step();
                v1=unit_doubled_distribution.var();
                unit_doubled_distribution.step();
                v2=unit_doubled_distribution.var();
                rsq=v1*v1+v2*v2;
            }
            while (rsq >= 1.0 || rsq == 0.0);
            fac=sqrt(-2.0*log(rsq)/rsq);
            gset=v1*fac;
            iset=1;
            return par_mu + v2*fac*par_sigma;
        }
        else {
            iset=0;
            return par_mu + gset*par_sigma;
        }
        return 0.0;
    }
public:
    gasdev () {
        modelname = "gasdev_noise";
    }
    ~gasdev() = default;
    std::string get_basename () {
        std::stringstream name_temp;
        name_temp << "_" << modelname << "_s_" << par_sigma << "_m_" << par_mu << "_idum_" << idum;
        return name_temp.str();
    }
    void set_pars (int set_size, double set_par_mu, double set_par_sigma, long set_idum) {
        if (set_size>0) {
            val.resize(set_size);
        }
        else {
            throw std::invalid_argument("Vector size must be positive (from gasdev)\n");
        }
        par_mu = set_par_mu;
        if (set_par_sigma>.0) {
            par_sigma = set_par_sigma;
        }
        else {
            throw std::invalid_argument("Noise scale must be positive (from gasdev)\n");
        }
        if (set_idum<0) {
            idum = set_idum;
            unit_doubled_distribution.set_pars (2, 0.0, 2.0, idum);
        }
        else {
            throw std::invalid_argument("Seed idum must be positive (from gasdev)\n");
        }
        /*Инициализация генератора чисел*/
        iset = 0;
        /*Инициализация генератора чисел*/
        initialized = true;
    }
};

#endif