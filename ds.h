#ifndef DS_H
#define DS_H

#include <iostream>
#include <vector>

#include "process.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "buffer.h"

class ds : public process {
protected:
    bool calc_lyap;
    std::vector<double> lyaps_sum;
    std::vector<double> lyaps_time;
    std::vector<double> lyaps_time_prev;
    array3d<double> all_vars_temp;
    array1d<buffer> vars;
    array2d<buffer> vars_tilda;
    double dt;
public:
    ds ();
    ~ds() = default;
    virtual void function (long, double, long) = 0;
    void resize_and_init (long set_size, long set_nvars, long set_nLCE, long set_nk);
    void calc_amplitude ();
    void normalization_init ();
    void normalization_reset ();
    void GramSmidt ();
    void lyap_init () { calc_lyap = true; }
    bool is_LCE_calced () { return calc_lyap; }
    double l(long i) { if (i<vars_tilda.s1()) return lyaps_sum[i]/lyaps_time[i]; return 0.0; }
    double getdt () { return dt; }
    double get_vars (long i1) { return vars.get(i1); }
    double get_time (long i1) { return vars.get_time(i1); }
    double get_entropy_vars (long i1) { return vars.get_entropy(i1); }
    double get_entropy2_vars (long i1) { return vars.get_entropy2(i1); }
    double get_prob_dens_vars (long i1, long win_number) { return vars.get_prob_dens(i1, win_number); }
    double get_prob_right_vars (long i1, long win_number) { return vars.get_prob_right(i1, win_number); }
    double get_prob_left_vars (long i1, long win_number) { return vars.get_prob_left(i1, win_number); }
    void prob_init_vars (long i1, long set_win_size) { vars.prob_init(i1, set_win_size); }
    void update_stat_vars (long i1) { vars.update_stat(i1); }
    long lyapnumber() { return vars_tilda.s1(); }
    long varnumber() { return vars.s1(); }
    void init_minima (long i1) { vars.init_minima(i1); }
    void init_minimaRR (long i1) { vars.init_minimaRR(i1); }
    void init_maxima (long i1) { vars.init_maxima(i1); }
    void init_maximaRR (long i1) { vars.init_maximaRR(i1); }
    void init_minima (long i1, double value_to_send) { vars.init_minima(i1, value_to_send); }
    void init_minimaRR (long i1, double value_to_send) { vars.init_minimaRR(i1, value_to_send); }
    void init_maxima (long i1, double value_to_send) { vars.init_maxima(i1, value_to_send); }
    void init_maximaRR (long i1, double value_to_send) { vars.init_maximaRR(i1, value_to_send); }
    void get_minima (long i1, bool &result, double &time, double &value) { vars.get_minima(i1, result, time, value); }
    void get_minimaRR (long i1, bool &result, double &time, double &value) { vars.get_minimaRR(i1, result, time, value); }
    void get_maxima (long i1, bool &result, double &time, double &value) { vars.get_maxima(i1, result, time, value); }
    void get_maximaRR (long i1, bool &result, double &time, double &value) { vars.get_maximaRR(i1, result, time, value); }
    long get_bins_number (long i1) { return vars.get_bins_number(i1); }
    long get_bins_number_minima (long i1) { return vars.get_bins_number_minima(i1); }
    long get_bins_number_minimaRR (long i1) { return vars.get_bins_number_minimaRR(i1); }
    long get_bins_number_maxima (long i1) { return vars.get_bins_number_maxima(i1); }
    long get_bins_number_maximaRR (long i1) { return vars.get_bins_number_maximaRR(i1); }
    double stat_get_mean_minima (long i1) const {return vars.stat_get_mean_minima (i1);}
    double stat_get_mean_maxima (long i1) const {return vars.stat_get_mean_maxima (i1);}
    double stat_get_mean_minimaRR (long i1) const {return vars.stat_get_mean_minimaRR (i1);}
    double stat_get_mean_maximaRR (long i1) const {return vars.stat_get_mean_maximaRR (i1);}
    double stat_get_disp_minima (long i1) const {return vars.stat_get_disp_minima (i1);}
    double stat_get_disp_maxima (long i1) const {return vars.stat_get_disp_maxima (i1);}
    double stat_get_disp_minimaRR (long i1) const {return vars.stat_get_disp_minimaRR (i1);}
    double stat_get_disp_maximaRR (long i1) const {return vars.stat_get_disp_maximaRR (i1);}
    double stat_get_disp_root_minima (long i1) const {return vars.stat_get_disp_root_minima (i1);}
    double stat_get_disp_root_maxima (long i1) const {return vars.stat_get_disp_root_maxima (i1);}
    double stat_get_disp_root_minimaRR (long i1) const {return vars.stat_get_disp_root_minimaRR (i1);}
    double stat_get_disp_root_maximaRR (long i1) const {return vars.stat_get_disp_root_maximaRR (i1);}
    double stat_get_skew_minima (long i1) const {return vars.stat_get_skew_minima (i1);}
    double stat_get_skew_maxima (long i1) const {return vars.stat_get_skew_maxima (i1);}
    double stat_get_skew_minimaRR (long i1) const {return vars.stat_get_skew_minimaRR (i1);}
    double stat_get_skew_maximaRR (long i1) const {return vars.stat_get_skew_maximaRR (i1);}
    double stat_get_kurt_minima (long i1) const {return vars.stat_get_kurt_minima (i1);}
    double stat_get_kurt_maxima (long i1) const {return vars.stat_get_kurt_maxima (i1);}
    double stat_get_kurt_minimaRR (long i1) const {return vars.stat_get_kurt_minimaRR (i1);}
    double stat_get_kurt_maximaRR (long i1) const {return vars.stat_get_kurt_maximaRR (i1);}
    double stat_get_correlation_minima (long i1) const {return vars.stat_get_correlation_minima (i1);}
    double stat_get_correlation_maxima (long i1) const {return vars.stat_get_correlation_maxima (i1);}
    double stat_get_correlation_minimaRR (long i1) const {return vars.stat_get_correlation_minimaRR (i1);}
    double stat_get_correlation_maximaRR (long i1) const {return vars.stat_get_correlation_maximaRR (i1);}
    double get_prob_dens_minima (long i1, long i) const {return vars.get_prob_dens_minima (i1, i);}
    double get_prob_dens_maxima (long i1, long i) const {return vars.get_prob_dens_maxima (i1, i);}
    double get_prob_dens_minimaRR (long i1, long i) const {return vars.get_prob_dens_minimaRR (i1, i);}
    double get_prob_dens_maximaRR (long i1, long i) const {return vars.get_prob_dens_maximaRR (i1, i);}
    double get_prob_left_minima (long i1, long i) const {return vars.get_prob_left_minima (i1, i);}
    double get_prob_left_maxima (long i1, long i) const {return vars.get_prob_left_maxima (i1, i);}
    double get_prob_left_minimaRR (long i1, long i) const {return vars.get_prob_left_minimaRR (i1, i);}
    double get_prob_left_maximaRR (long i1, long i) const {return vars.get_prob_left_maximaRR (i1, i);}
    double get_prob_right_minima (long i1, long i) const {return vars.get_prob_right_minima (i1, i);}
    double get_prob_right_maxima (long i1, long i) const {return vars.get_prob_right_maxima (i1, i);}
    double get_prob_right_minimaRR (long i1, long i) const {return vars.get_prob_right_minimaRR (i1, i);}
    double get_prob_right_maximaRR (long i1, long i) const {return vars.get_prob_right_maximaRR (i1, i);}
    double get_prob_mean_minima (long i1, long i) const {return vars.get_prob_mean_minima (i1, i);}
    double get_prob_mean_maxima (long i1, long i) const {return vars.get_prob_mean_maxima (i1, i);}
    double get_prob_mean_minimaRR (long i1, long i) const {return vars.get_prob_mean_minimaRR (i1, i);}
    double get_prob_mean_maximaRR (long i1, long i) const {return vars.get_prob_mean_maximaRR (i1, i);}
    double get_prob_disp_minima (long i1, long i) const {return vars.get_prob_disp_minima (i1, i);}
    double get_prob_disp_maxima (long i1, long i) const {return vars.get_prob_disp_maxima (i1, i);}
    double get_prob_disp_minimaRR (long i1, long i) const {return vars.get_prob_disp_minimaRR (i1, i);}
    double get_prob_disp_maximaRR (long i1, long i) const {return vars.get_prob_disp_maximaRR (i1, i);}
    double get_entropy_minima (long i1) const {return vars.get_entropy_minima (i1);}
    double get_entropy_maxima (long i1) const {return vars.get_entropy_maxima (i1);}
    double get_entropy_minimaRR (long i1) const {return vars.get_entropy_minimaRR (i1);}
    double get_entropy_maximaRR (long i1) const {return vars.get_entropy_maximaRR (i1);}
    double get_entropy2_minima (long i1) const {return vars.get_entropy2_minima (i1);}
    double get_entropy2_maxima (long i1) const {return vars.get_entropy2_maxima (i1);}
    double get_entropy2_minimaRR (long i1) const {return vars.get_entropy2_minimaRR (i1);}
    double get_entropy2_maximaRR (long i1) const {return vars.get_entropy2_maximaRR (i1);}
    double extreme_part_inside_minima (long i1, double sigma_in) const {return vars.extreme_part_inside_minima (i1, sigma_in);}
    double extreme_part_inside_maxima (long i1, double sigma_in) const {return vars.extreme_part_inside_maxima (i1, sigma_in);}
    double extreme_part_inside_minimaRR (long i1, double sigma_in) const {return vars.extreme_part_inside_minimaRR (i1, sigma_in);}
    double extreme_part_inside_maximaRR (long i1, double sigma_in) const {return vars.extreme_part_inside_maximaRR (i1, sigma_in);}
    double extreme_part_outside_minima (long i1, double sigma_in) const {return vars.extreme_part_outside_minima (i1, sigma_in);}
    double extreme_part_outside_maxima (long i1, double sigma_in) const {return vars.extreme_part_outside_maxima (i1, sigma_in);}
    double extreme_part_outside_minimaRR (long i1, double sigma_in) const {return vars.extreme_part_outside_minimaRR (i1, sigma_in);}
    double extreme_part_outside_maximaRR (long i1, double sigma_in) const {return vars.extreme_part_outside_maximaRR (i1, sigma_in);}
};

ds::ds () {
    calc_lyap=false;
    dt = 0.0;
    modelname = "ds";
}

void ds::resize_and_init (long set_size, long set_nvars, long set_nLCE, long set_nk) {
    if (set_size<1)
        throw std::invalid_argument("Vector size must be positive (from ds)\n");
    if (set_nvars<1)
        throw std::invalid_argument("Model size must be positive (from ds)\n");
    if (set_nk<1)
        throw std::invalid_argument("Temporal buffer must be positive (from ds)\n");
    if (set_nLCE<0)
        throw std::invalid_argument("LCE buffer must be non-negative (from ds)\n");
    val.resize(set_size, 0.0);
    vars.resize(set_size, set_nvars);
    all_vars_temp.resize(set_nk, 1+set_nLCE, set_nvars);
    if ((calc_lyap)&&(set_nLCE>0)) {
        vars_tilda.resize(set_size, set_nLCE, set_nvars);
        lyaps_sum.resize(set_nLCE, 0.0);
        lyaps_time.resize(set_nLCE, 0.0);
        lyaps_time_prev.resize(set_nLCE, 0.0);
    }
}

void ds::calc_amplitude () {
    double sumsq = 0.0;
    for (long n=0;n<vars.s1();n++) {
        sumsq += vars.get(n)*vars.get(n);
    }
    val[marker] = sqrt(sumsq);
}

void ds::normalization_init () {
    for (long l=0;l<vars_tilda.s1();l++) {
        for (long n=0;n<vars_tilda.s2();n++) {
            for (long h=0;h<vars_tilda.s0();h++) {
                if (l==n) vars_tilda.step(l, n, 1.0);
                else vars_tilda.step(l, n, 0.0);
            }
        }
        lyaps_sum[l] = 0.0;
        lyaps_time[l] = 0.0;
        lyaps_time_prev[l] = discrete_time;
    }
}

void ds::normalization_reset () {
    GramSmidt ();
    for (long i=0;i<vars_tilda.s1();i++) {
        lyaps_sum[i] = 0.0;
        lyaps_time[i] = 0.0;
        lyaps_time_prev[i] = discrete_time;
    }
}

void ds::GramSmidt () {
    for (long l=0;l<vars_tilda.s1();l++) {
        for (long j=0;j<l;j++) {
            double prod = 0.0;
            for (long h=0;h<vars_tilda.s0();h++) {
                for (long n=0;n<vars_tilda.s2();n++) {
                    prod += vars_tilda.get_delayed(h, l, n) * vars_tilda.get_delayed(h, j, n);
                }
            }
            for (long h=0;h<vars_tilda.s0();h++) {
                for (long n=0;n<vars_tilda.s2();n++) {
                    vars_tilda.correct(h, l, n, vars_tilda.get_delayed(h, l, n) - prod * vars_tilda.get_delayed(h, j, n));
                }
            }
        }
        double normsq = 0.0;
        for (long h=0;h<vars_tilda.s0();h++) {
            for (long n=0;n<vars_tilda.s2();n++) {
                double temp_val = vars_tilda.get_delayed(h, l, n);
                normsq += temp_val * temp_val;
            }
        }
        double norm = sqrt(normsq);
        for (long h=0;h<vars_tilda.s0();h++) {
            for (long n=0;n<vars_tilda.s2();n++) {
                vars_tilda.correct(h, l, n, vars_tilda.get_delayed(h, l, n)/norm);
            }
        }
        lyaps_sum[l] += log(norm);
        lyaps_time[l] += (discrete_time-lyaps_time_prev[l])*dt;
        lyaps_time_prev[l] = discrete_time;
    }
}

#endif