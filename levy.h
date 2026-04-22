//Убрать сложные методы в конец файла

#ifndef LEVY_H
#define LEVY_H

#include "noise.h"
#include "ran.h"
#include "expdev.h"

class levy : public noise {
protected:
    double par_mu;
    double par_alpha;
    double par_beta;
    double par_sigma;
    /*Статические параметры генератора чисел*/
    float arg3, arg1_v2, b_ab, s_ab;
    /*Статические параметры генератора чисел*/
    ran pi_distribution;
    expdev exp_distribution;
    double generator () {
        pi_distribution.step();
        float uninum = pi_distribution.var();
        exp_distribution.step();
        float expnum = exp_distribution.var();
        if (fabs(par_alpha-1.0)<1.0e-10) {
            float arg1 = (0.5 * 3.141592) + (par_beta * uninum);
            float arg2 = expnum * cos(uninum);
            float x = (2.0 / 3.141592) * (arg1 * tan(uninum) - (par_beta * log(arg2 / arg1)));
            return par_mu + par_sigma * (x + arg3) / sqrt(2.0);
        }
        else {
            float arg2 = par_alpha * (uninum + b_ab);
            float n1 = sin(arg2);
            float d1 = exp(log(cos(uninum))/par_alpha);
            float n2 = cos(uninum - arg2);
            float x = s_ab * (n1/d1) * pow((n2/expnum),((1.0-par_alpha)/par_alpha));
            return par_mu + (par_sigma * x) / sqrt(2.0);
        }
        return 0.0;
    }
public:
    levy () {
        modelname = "Levy_noise";
    }
    ~levy() = default;
    std::string get_basename () {
        std::stringstream name_temp;
        name_temp << "_" << modelname << "_a_" << par_alpha << "_b_" << par_beta << "_s_" << par_sigma << "_m_" << par_mu << "_idum_" << idum;
        return name_temp.str();
    }
    void set_pars (int set_size, double set_par_mu, double set_par_alpha, double set_par_beta, double set_par_sigma, long set_idum) {
        if (set_size>0) {
            val.resize(set_size);
        }
        else {
            throw std::invalid_argument("Vector size must be positive (from levy)\n");
        }
        par_mu = set_par_mu;
        if ( (set_par_alpha>.0)&&(set_par_alpha<(2.0+1.0e-10)) ) {
            par_alpha = set_par_alpha;
        }
        else {
            throw std::invalid_argument("Noise stability parameter must be in (0;2] (from levy)\n");
        }
        if ( (set_par_beta>(-1.0+1.0e-10)) && (set_par_beta<(1.0-1.0e-10)) ) {
            par_beta = set_par_beta;
        }
        else {
            throw std::invalid_argument("Noise asymmetry paramter must be in (-1:1) (from levy)\n");
        }
        if (set_par_sigma>.0) {
            par_sigma = set_par_sigma;
        }
        else {
            throw std::invalid_argument("Noise scale must be positive (from levy)\n");
        }
        if (set_idum<0) {
            idum = set_idum;
            pi_distribution.set_pars (2, 0.0, 3.141592, idum);
            exp_distribution.set_pars (2, 1.0, idum*2-1);
        }
        else {
            throw std::invalid_argument("Seed idum must be positive (from levy)\n");
        }
        /*Инициализация генератора чисел*/
        arg3 = (2.0 / 3.141592) * par_beta * log(2.0 / 3.141592);
        arg1_v2 = 0.5 * 3.141592 * par_alpha;
        b_ab = atan(par_beta * tan(arg1_v2)) / par_alpha;
        s_ab = pow((1.0 + (par_beta*par_beta) * tan(arg1_v2)*tan(arg1_v2)), (1.0/(2.0*par_alpha)));
        /*Инициализация генератора чисел*/
        initialized = true;
    }
    void change_beta (double set_par_beta) {
        if ( (set_par_beta>(-1.0-1.0e-10)) && (set_par_beta<(1.0+1.0e-10)) ) {
            par_beta = set_par_beta;
        }
        else {
            throw std::invalid_argument("New noise asymmetry paramter must be in [-1:1] (from levy)\n");
        }
        arg3 = (2.0 / 3.141592) * par_beta * log(2.0 / 3.141592);
        b_ab = atan(par_beta * tan(arg1_v2)) / par_alpha;
        s_ab = pow((1.0 + (par_beta*par_beta) * tan(arg1_v2)*tan(arg1_v2)), (1.0/(2.0*par_alpha)));
    }
    void change_alpha (double set_par_alpha) {
        if ( (set_par_alpha>.0)&&(set_par_alpha<(2.0+1.0e-10)) ) {
            par_alpha = set_par_alpha;
        }
        else {
            throw std::invalid_argument("New noise stability parameter must be in (0;2] (from levy)\n");
        }
        arg1_v2 = 0.5 * 3.141592 * par_alpha;
        b_ab = atan(par_beta * tan(arg1_v2)) / par_alpha;
        s_ab = pow((1.0 + (par_beta*par_beta) * tan(arg1_v2)*tan(arg1_v2)), (1.0/(2.0*par_alpha)));
    }
};

#endif