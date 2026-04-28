#ifndef OBSERVERSE_H
#define OBSERVERSE_H

#include "statistics.h"
#include "correlation.h"
#include "distribution.h"

class observers {
protected:
    bool use_statistics;
    bool use_correlation;
    bool use_prob_distribution;
    std::vector<statistics> local_stat;
    std::vector<correlation> local_corr;
    std::vector<distribution> local_distr;
public:
    observers ();
    ~observers() = default;
    void corr_init ();
    void stat_init ();
    void prob_init ();
    void prob_init (long set_size);
    void update_values (double x_val, double y_val, bool buffer_is_full);
    double stat_get_mean () const {if (!use_statistics) return 0.0; return local_stat[0].stat_get_mean ();}
    double stat_get_disp_root () const {if (!use_statistics) return 0.0; return local_stat[0].stat_get_disp_root ();}
    double stat_get_disp () const {if (!use_statistics) return 0.0; return local_stat[0].stat_get_disp ();}
    double stat_get_skew () const {if (!use_statistics) return 0.0; return local_stat[0].stat_get_skew ();}
    double stat_get_kurt () const {if (!use_statistics) return 0.0; return local_stat[0].stat_get_kurt ();}
    double stat_get_correlation () const {if (!use_correlation) return 0.0; return local_corr[0].stat_get_correlation ();}
    double get_prob_dens (long i) const {if (!use_prob_distribution) return 0.0; return local_distr[0].get_prob_dens(i);}
    double get_prob_left (long i) const {if (!use_prob_distribution) return 0.0; return local_distr[0].get_prob_left(i);}
    double get_prob_right (long i) const {if (!use_prob_distribution) return 0.0; return local_distr[0].get_prob_right(i);}
    double get_prob_mean (long i) const {if (!use_prob_distribution) return 0.0; return local_distr[0].get_prob_mean(i);}
    double get_prob_disp (long i) const {if (!use_prob_distribution) return 0.0; return local_distr[0].get_prob_disp(i);}
    double get_entropy () const {if (!use_prob_distribution) return 0.0; return local_distr[0].get_entropy();}
    double get_entropy2 () const {if (!use_prob_distribution) return 0.0; return local_distr[0].get_entropy()/log(2.0);}
    double extreme_part_inside (double sigma_in) const {if (!use_prob_distribution) return 0.0; return local_distr[0].extreme_part_inside(sigma_in);}
    double extreme_part_outside (double sigma_in) const {if (!use_prob_distribution) return 0.0; return local_distr[0].extreme_part_outside(sigma_in);}
    long get_bins_number () {if (!use_prob_distribution) return 0.0; return local_distr[0].get_bins_number();}
};

observers::observers () {
    use_statistics = false;
    use_correlation = false;
    use_prob_distribution = false;
}

void observers::corr_init () {
    local_corr.resize(1);
    use_correlation = true;
}

void observers::stat_init () {
    local_stat.resize(1);
    use_statistics = true;
}

void observers::prob_init () {
    stat_init ();
    local_distr.resize(1);
    use_prob_distribution = true;
}

void observers::prob_init (long set_size) {
    stat_init ();
    local_distr.resize(1);
    local_distr[0].resize(set_size);
    use_prob_distribution = true;
}

void observers::update_values (double x_val, double y_val, bool buffer_is_full) {
    if (use_statistics)
        local_stat[0].update_stat (x_val);
    if (use_prob_distribution) {
        local_distr[0].calc_probability (x_val);
    }
    if (use_correlation && buffer_is_full) {
        local_corr[0].update_corr (x_val, y_val);
    }
}

#endif