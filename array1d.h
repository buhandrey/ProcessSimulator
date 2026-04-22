#ifndef ARRAY1D_H
#define ARRAY1D_H

#include <iostream>
#include <vector>

template <typename array_element>
class array1d {
protected:
    long size0;
    long size1;
    std::vector<array_element> val;
public:
    array1d () {
        size0 = 0;
        size1 = 0;
    }
    ~array1d() = default;
    void resize(long set_size0, long set_size1) {
        size0 = set_size0;
        size1 = set_size1;
        val.resize(size1);
        for (long s1=0;s1<size1;s1++) {
            val[s1].resize(size0);
            val[s1].init();
        }
    }
    long s0 () { return size0; }
    long s1 () { return size1; }
    double get (long i1) { return val[i1].var(); }
    double get_delayed (long delay, long i1) { return val[i1].var(delay); }
    double get_time (long i1) { return val[i1].time_int(); }
    double get_entropy (long i1) { return val[i1].get_entropy(); }
    double get_entropy2 (long i1) { return val[i1].get_entropy2(); }
    double get_prob_dens (long i1, long win_number) { return val[i1].get_prob_dens(win_number); }
    double get_prob_right (long i1, long win_number) { return val[i1].get_prob_right(win_number); }
    double get_prob_left (long i1, long win_number) { return val[i1].get_prob_left(win_number); }
    void initbuffer (long i1, double value) { val[i1].initwithval(value); }
    void step (long i1, double value) { val[i1].step(value); }
    void correct (long i1, double value) { val[i1].correct(value); }
    void correct (long delay, long i1, double value) { val[i1].correct(delay, value); }
    void prob_init (long i1, long set_win_size) { val[i1].prob_init(set_win_size); }
    void update_stat (long i1) { val[i1].update_stat(); }
    long get_bins_number (long i1) { return val[i1].get_bins_number(); }
    long get_bins_number_minima (long i1) { return val[i1].get_bins_number_minima(); }
    long get_bins_number_minimaRR (long i1) { return val[i1].get_bins_number_minimaRR(); }
    long get_bins_number_maxima (long i1) { return val[i1].get_bins_number_maxima(); }
    long get_bins_number_maximaRR (long i1) { return val[i1].get_bins_number_maximaRR(); }
    void init_minima (long var) { val[var].init_minima(); }
    void init_minimaRR (long var) { val[var].init_minimaRR(); }
    void init_maxima (long var) { val[var].init_maxima(); }
    void init_maximaRR (long var) { val[var].init_maximaRR(); }
    void init_minima (long var, double set_limit_value_new) { val[var].init_minima(set_limit_value_new); }
    void init_minimaRR (long var, double set_limit_value_new) { val[var].init_minimaRR(set_limit_value_new); }
    void init_maxima (long var, double set_limit_value_new) { val[var].init_maxima(set_limit_value_new); }
    void init_maximaRR (long var, double set_limit_value_new) { val[var].init_maximaRR(set_limit_value_new); }
    void get_minima (long var, bool &result, double &time, double &value) { val[var].get_minima(result, time, value); }
    void get_minimaRR (long var, bool &result, double &time, double &value) { val[var].get_minimaRR(result, time, value); }
    void get_maxima (long var, bool &result, double &time, double &value) { val[var].get_maxima(result, time, value); }
    void get_maximaRR (long var, bool &result, double &time, double &value) { val[var].get_maximaRR(result, time, value); }
    double stat_get_mean_minima (long var) const {return val[var].stat_get_mean_minima ();}
    double stat_get_mean_maxima (long var) const {return val[var].stat_get_mean_maxima ();}
    double stat_get_mean_minimaRR (long var) const {return val[var].stat_get_mean_minimaRR ();}
    double stat_get_mean_maximaRR (long var) const {return val[var].stat_get_mean_maximaRR ();}
    double stat_get_disp_minima (long var) const {return val[var].stat_get_disp_minima ();}
    double stat_get_disp_maxima (long var) const {return val[var].stat_get_disp_maxima ();}
    double stat_get_disp_minimaRR (long var) const {return val[var].stat_get_disp_minimaRR ();}
    double stat_get_disp_maximaRR (long var) const {return val[var].stat_get_disp_maximaRR ();}
    double stat_get_disp_root_minima (long var) const {return val[var].stat_get_disp_root_minima ();}
    double stat_get_disp_root_maxima (long var) const {return val[var].stat_get_disp_root_maxima ();}
    double stat_get_disp_root_minimaRR (long var) const {return val[var].stat_get_disp_root_minimaRR ();}
    double stat_get_disp_root_maximaRR (long var) const {return val[var].stat_get_disp_root_maximaRR ();}
    double stat_get_skew_minima (long var) const {return val[var].stat_get_skew_minima ();}
    double stat_get_skew_maxima (long var) const {return val[var].stat_get_skew_maxima ();}
    double stat_get_skew_minimaRR (long var) const {return val[var].stat_get_skew_minimaRR ();}
    double stat_get_skew_maximaRR (long var) const {return val[var].stat_get_skew_maximaRR ();}
    double stat_get_kurt_minima (long var) const {return val[var].stat_get_kurt_minima ();}
    double stat_get_kurt_maxima (long var) const {return val[var].stat_get_kurt_maxima ();}
    double stat_get_kurt_minimaRR (long var) const {return val[var].stat_get_kurt_minimaRR ();}
    double stat_get_kurt_maximaRR (long var) const {return val[var].stat_get_kurt_maximaRR ();}
    double stat_get_correlation_minima (long var) const {return val[var].stat_get_correlation_minima ();}
    double stat_get_correlation_maxima (long var) const {return val[var].stat_get_correlation_maxima ();}
    double stat_get_correlation_minimaRR (long var) const {return val[var].stat_get_correlation_minimaRR ();}
    double stat_get_correlation_maximaRR (long var) const {return val[var].stat_get_correlation_maximaRR ();}
    double get_prob_dens_minima (long var, long i) const {return val[var].get_prob_dens_minima (i);}
    double get_prob_dens_maxima (long var, long i) const {return val[var].get_prob_dens_maxima (i);}
    double get_prob_dens_minimaRR (long var, long i) const {return val[var].get_prob_dens_minimaRR (i);}
    double get_prob_dens_maximaRR (long var, long i) const {return val[var].get_prob_dens_maximaRR (i);}
    double get_prob_left_minima (long var, long i) const {return val[var].get_prob_left_minima (i);}
    double get_prob_left_maxima (long var, long i) const {return val[var].get_prob_left_maxima (i);}
    double get_prob_left_minimaRR (long var, long i) const {return val[var].get_prob_left_minimaRR (i);}
    double get_prob_left_maximaRR (long var, long i) const {return val[var].get_prob_left_maximaRR (i);}
    double get_prob_right_minima (long var, long i) const {return val[var].get_prob_right_minima (i);}
    double get_prob_right_maxima (long var, long i) const {return val[var].get_prob_right_maxima (i);}
    double get_prob_right_minimaRR (long var, long i) const {return val[var].get_prob_right_minimaRR (i);}
    double get_prob_right_maximaRR (long var, long i) const {return val[var].get_prob_right_maximaRR (i);}
    double get_prob_mean_minima (long var, long i) const {return val[var].get_prob_mean_minima (i);}
    double get_prob_mean_maxima (long var, long i) const {return val[var].get_prob_mean_maxima (i);}
    double get_prob_mean_minimaRR (long var, long i) const {return val[var].get_prob_mean_minimaRR (i);}
    double get_prob_mean_maximaRR (long var, long i) const {return val[var].get_prob_mean_maximaRR (i);}
    double get_prob_disp_minima (long var, long i) const {return val[var].get_prob_disp_minima (i);}
    double get_prob_disp_maxima (long var, long i) const {return val[var].get_prob_disp_maxima (i);}
    double get_prob_disp_minimaRR (long var, long i) const {return val[var].get_prob_disp_minimaRR (i);}
    double get_prob_disp_maximaRR (long var, long i) const {return val[var].get_prob_disp_maximaRR (i);}
    double get_entropy_minima (long var) const {return val[var].get_entropy_minima ();}
    double get_entropy_maxima (long var) const {return val[var].get_entropy_maxima ();}
    double get_entropy_minimaRR (long var) const {return val[var].get_entropy_minimaRR ();}
    double get_entropy_maximaRR (long var) const {return val[var].get_entropy_maximaRR ();}
    double get_entropy2_minima (long var) const {return val[var].get_entropy2_minima ();}
    double get_entropy2_maxima (long var) const {return val[var].get_entropy2_maxima ();}
    double get_entropy2_minimaRR (long var) const {return val[var].get_entropy2_minimaRR ();}
    double get_entropy2_maximaRR (long var) const {return val[var].get_entropy2_maximaRR ();}
    double extreme_part_inside_minima (long var, double sigma_in) const {return val[var].extreme_part_inside_minima (sigma_in);}
    double extreme_part_inside_maxima (long var, double sigma_in) const {return val[var].extreme_part_inside_maxima (sigma_in);}
    double extreme_part_inside_minimaRR (long var, double sigma_in) const {return val[var].extreme_part_inside_minimaRR (sigma_in);}
    double extreme_part_inside_maximaRR (long var, double sigma_in) const {return val[var].extreme_part_inside_maximaRR (sigma_in);}
    double extreme_part_outside_minima (long var, double sigma_in) const {return val[var].extreme_part_outside_minima (sigma_in);}
    double extreme_part_outside_maxima (long var, double sigma_in) const {return val[var].extreme_part_outside_maxima (sigma_in);}
    double extreme_part_outside_minimaRR (long var, double sigma_in) const {return val[var].extreme_part_outside_minimaRR (sigma_in);}
    double extreme_part_outside_maximaRR (long var, double sigma_in) const {return val[var].extreme_part_outside_maximaRR (sigma_in);}
};

#endif