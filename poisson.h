#ifndef POISSON_H
#define POISSON_H

#include "noise.h"
#include "ran.h"
#include "events.h"

class poisson : public noise {
protected:
    double par_freq;
    long prev_imp_time;
    ran unit_distribution;
    bool calc_events;
    events local_events;
    double generator ();
public:
    poisson ();
    ~poisson() = default;
    std::string get_basename ();
    void set_pars (int set_size, double set_par_freq, long set_idum);
    long get_bins_number_events () { return local_events.get_bins_number(); }
    void init_events () { calc_events = true; local_events.resize(100); local_events.init(); }
    void get_events (bool &result, double &time, double &value) { local_events.get(result, time, value); }
    void is_there_event (bool &result) { local_events.is_there_event(result); }
    double stat_get_mean_events () const {if (!calc_events) return 0.0; return local_events.stat_get_mean ();}
    double stat_get_disp_events () const {if (!calc_events) return 0.0; return local_events.stat_get_disp ();}
    double stat_get_disp_root_events () const {if (!calc_events) return 0.0; return local_events.stat_get_disp_root ();}
    double stat_get_skew_events () const {if (!calc_events) return 0.0; return local_events.stat_get_skew ();}
    double stat_get_kurt_events () const {if (!calc_events) return 0.0; return local_events.stat_get_kurt ();}
    double stat_get_correlation_events () const {if (!calc_events) return 0.0; return local_events.stat_get_correlation ();}
    double get_prob_dens_events (long i) const {if (!calc_events) return 0.0; return local_events.get_prob_dens (i);}
    double get_prob_left_events (long i) const {if (!calc_events) return 0.0; return local_events.get_prob_left (i);}
    double get_prob_right_events (long i) const {if (!calc_events) return 0.0; return local_events.get_prob_right (i);}
    double get_prob_mean_events (long i) const {if (!calc_events) return 0.0; return local_events.get_prob_mean (i);}
    double get_prob_disp_events (long i) const {if (!calc_events) return 0.0; return local_events.get_prob_disp (i);}
    double get_entropy_events () const {if (!calc_events) return 0.0; return local_events.get_entropy ();}
    double get_entropy2_events () const {if (!calc_events) return 0.0; return local_events.get_entropy2 ();}
    double extreme_part_inside_events (double sigma_in) const {if (!calc_events) return 0.0; return local_events.extreme_part_inside (sigma_in);}
    double extreme_part_outside_events (double sigma_in) const {if (!calc_events) return 0.0; return local_events.extreme_part_outside (sigma_in);}
};

double poisson::generator () {
    double dum;
    unit_distribution.step();
    dum = unit_distribution.var();
    if ((dum>(0.5-0.5*par_freq))&&(dum<(0.5+0.5*par_freq))) {
        prev_imp_time = discrete_time;
        if (calc_events) local_events.step(1.0);
        return 1.0;
    }
    else {
        if (calc_events) local_events.step(0.0);
        return 0.0;
    }
}

poisson::poisson () {
    modelname = "Poisson_signal";
    prev_imp_time = -1;
    calc_events = false;
}

std::string poisson::get_basename () {
    std::stringstream name_temp;
    name_temp << "_" << modelname << "_freq_" << par_freq << "_idum_" << idum;
    return name_temp.str();
}

void poisson::set_pars (int set_size, double set_par_freq, long set_idum) {
    if (set_size>0) {
        val.resize(set_size);
    }
    else {
        throw std::invalid_argument("Vector size must be positive (from poisson)\n");
    }
    if (set_par_freq>.0) {
        par_freq = set_par_freq;
    }
    else {
        throw std::invalid_argument("Poisson signal frequency must be positive (from poisson)\n");
    }
    if (set_idum<0) {
        idum = set_idum;
        unit_distribution.set_pars (2, 0.5, 1.0, idum);
    }
    else {
        throw std::invalid_argument("Seed idum must be negative (from poisson)\n");
    }
    initialized = true;
}

#endif