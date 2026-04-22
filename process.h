#ifndef PROCESS_H
#define PROCESS_H

#include <iostream>
#include <vector>

#include "observers.h"
#include "minima.h"
#include "minimaRR.h"
#include "maxima.h"
#include "maximaRR.h"

class process : public observers {
protected:
    long marker;
    long discrete_time;
    bool initialized;
    std::vector<double> val;
    bool calc_minima;
    bool calc_minimaRR;
    bool calc_maxima;
    bool calc_maximaRR;
    minima local_minima;
    minimaRR local_minimaRR;
    maxima local_maxima;
    maximaRR local_maximaRR;
    std::string modelname;
public:
    process () {
        marker = -1;
        discrete_time = -1;
        initialized = false;
        calc_minima = false;
        calc_minimaRR = false;
        calc_maxima = false;
        calc_maximaRR = false;
        modelname = "process";
    }
    ~process() = default;
    std::string get_modelname () { return modelname;}
    virtual std::string get_basename () = 0;
    void update_stat () {
        bool buffer_is_full = false;
        long marker_y = marker + 1;
        if (marker_y == val.size())
            marker_y = 0;
        if (discrete_time>val.size()) {
            buffer_is_full = true;
        }
        update_values(val[marker], val[marker_y], buffer_is_full);
    }
    virtual void step () = 0;
    double var() const { return val[marker]; }
    double var(long delay) const { if (delay<val.size()) { long d_marker = marker - delay; if (d_marker<0) d_marker+=val.size(); return val[d_marker];} else throw std::invalid_argument("Delay is larger than the buffer size (from process)\n"); return 0.0; }
    long time_int() const { return discrete_time; }
    long size() {
        return val.size();
    }
    void resize (long set_size) {
        val.resize(set_size, 0.0);
    }
    void init () {
        if (val.size()>1)
            initialized = true;
        else
            throw std::invalid_argument("The size of val should be not less than 2 (from process)\n");
    }
    void substep () {
        if (calc_minima)
            local_minima.step(val[marker]);
        if (calc_minimaRR)
            local_minimaRR.step(val[marker]);
        if (calc_maxima)
            local_maxima.step(val[marker]);
        if (calc_maximaRR)
            local_maximaRR.step(val[marker]);
    }
    long get_bins_number_minima () { return local_minima.get_bins_number(); }
    long get_bins_number_minimaRR () { return local_minimaRR.get_bins_number(); }
    long get_bins_number_maxima () { return local_maxima.get_bins_number(); }
    long get_bins_number_maximaRR () { return local_maximaRR.get_bins_number(); }
    void init_minima () { calc_minima = true; local_minima.resize(100); local_minima.init(); }
    void init_minimaRR () { calc_minimaRR = true; local_minimaRR.resize(100); local_minimaRR.init(); }
    void init_maxima () { calc_maxima = true; local_maxima.resize(100); local_maxima.init(); }
    void init_maximaRR () { calc_maximaRR = true; local_maximaRR.resize(100); local_maximaRR.init(); }
    void init_minima (double limit_value_new) { init_minima (); local_minima.set_limit_value(limit_value_new); }
    void init_minimaRR (double limit_value_new) { init_minimaRR (); local_minimaRR.set_limit_value(limit_value_new); }
    void init_maxima (double limit_value_new) { init_maxima (); local_maxima.set_limit_value(limit_value_new); }
    void init_maximaRR (double limit_value_new) { init_maximaRR (); local_maximaRR.set_limit_value(limit_value_new); }
    void get_minima (bool &result, double &time, double &value) { local_minima.get(result, time, value); }
    void get_minimaRR (bool &result, double &time, double &value) { local_minimaRR.get(result, time, value); }
    void get_maxima (bool &result, double &time, double &value) { local_maxima.get(result, time, value); }
    void get_maximaRR (bool &result, double &time, double &value) { local_maximaRR.get(result, time, value); }
    double stat_get_mean_minima () const {if (!calc_minima) return 0.0; return local_minima.stat_get_mean ();}
    double stat_get_mean_maxima () const {if (!calc_maxima) return 0.0; return local_maxima.stat_get_mean ();}
    double stat_get_mean_minimaRR () const {if (!calc_minimaRR) return 0.0; return local_minimaRR.stat_get_mean ();}
    double stat_get_mean_maximaRR () const {if (!calc_maximaRR) return 0.0; return local_maximaRR.stat_get_mean ();}
    double stat_get_disp_minima () const {if (!calc_minima) return 0.0; return local_minima.stat_get_disp ();}
    double stat_get_disp_maxima () const {if (!calc_maxima) return 0.0; return local_maxima.stat_get_disp ();}
    double stat_get_disp_minimaRR () const {if (!calc_minimaRR) return 0.0; return local_minimaRR.stat_get_disp ();}
    double stat_get_disp_maximaRR () const {if (!calc_maximaRR) return 0.0; return local_maximaRR.stat_get_disp ();}
    double stat_get_disp_root_minima () const {if (!calc_minima) return 0.0; return local_minima.stat_get_disp_root ();}
    double stat_get_disp_root_maxima () const {if (!calc_maxima) return 0.0; return local_maxima.stat_get_disp_root ();}
    double stat_get_disp_root_minimaRR () const {if (!calc_minimaRR) return 0.0; return local_minimaRR.stat_get_disp_root ();}
    double stat_get_disp_root_maximaRR () const {if (!calc_maximaRR) return 0.0; return local_maximaRR.stat_get_disp_root ();}
    double stat_get_skew_minima () const {if (!calc_minima) return 0.0; return local_minima.stat_get_skew ();}
    double stat_get_skew_maxima () const {if (!calc_maxima) return 0.0; return local_maxima.stat_get_skew ();}
    double stat_get_skew_minimaRR () const {if (!calc_minimaRR) return 0.0; return local_minimaRR.stat_get_skew ();}
    double stat_get_skew_maximaRR () const {if (!calc_maximaRR) return 0.0; return local_maximaRR.stat_get_skew ();}
    double stat_get_kurt_minima () const {if (!calc_minima) return 0.0; return local_minima.stat_get_kurt ();}
    double stat_get_kurt_maxima () const {if (!calc_maxima) return 0.0; return local_maxima.stat_get_kurt ();}
    double stat_get_kurt_minimaRR () const {if (!calc_minimaRR) return 0.0; return local_minimaRR.stat_get_kurt ();}
    double stat_get_kurt_maximaRR () const {if (!calc_maximaRR) return 0.0; return local_maximaRR.stat_get_kurt ();}
    double stat_get_correlation_minima () const {if (!calc_minima) return 0.0; return local_minima.stat_get_correlation ();}
    double stat_get_correlation_maxima () const {if (!calc_maxima) return 0.0; return local_maxima.stat_get_correlation ();}
    double stat_get_correlation_minimaRR () const {if (!calc_minimaRR) return 0.0; return local_minimaRR.stat_get_correlation ();}
    double stat_get_correlation_maximaRR () const {if (!calc_maximaRR) return 0.0; return local_maximaRR.stat_get_correlation ();}
    double get_prob_dens_minima (long i) const {if (!calc_minima) return 0.0; return local_minima.get_prob_dens (i);}
    double get_prob_dens_maxima (long i) const {if (!calc_maxima) return 0.0; return local_maxima.get_prob_dens (i);}
    double get_prob_dens_minimaRR (long i) const {if (!calc_minimaRR) return 0.0; return local_minimaRR.get_prob_dens (i);}
    double get_prob_dens_maximaRR (long i) const {if (!calc_maximaRR) return 0.0; return local_maximaRR.get_prob_dens (i);}
    double get_prob_left_minima (long i) const {if (!calc_minima) return 0.0; return local_minima.get_prob_left (i);}
    double get_prob_left_maxima (long i) const {if (!calc_maxima) return 0.0; return local_maxima.get_prob_left (i);}
    double get_prob_left_minimaRR (long i) const {if (!calc_minimaRR) return 0.0; return local_minimaRR.get_prob_left (i);}
    double get_prob_left_maximaRR (long i) const {if (!calc_maximaRR) return 0.0; return local_maximaRR.get_prob_left (i);}
    double get_prob_right_minima (long i) const {if (!calc_minima) return 0.0; return local_minima.get_prob_right (i);}
    double get_prob_right_maxima (long i) const {if (!calc_maxima) return 0.0; return local_maxima.get_prob_right (i);}
    double get_prob_right_minimaRR (long i) const {if (!calc_minimaRR) return 0.0; return local_minimaRR.get_prob_right (i);}
    double get_prob_right_maximaRR (long i) const {if (!calc_maximaRR) return 0.0; return local_maximaRR.get_prob_right (i);}
    double get_prob_mean_minima (long i) const {if (!calc_minima) return 0.0; return local_minima.get_prob_mean (i);}
    double get_prob_mean_maxima (long i) const {if (!calc_maxima) return 0.0; return local_maxima.get_prob_mean (i);}
    double get_prob_mean_minimaRR (long i) const {if (!calc_minimaRR) return 0.0; return local_minimaRR.get_prob_mean (i);}
    double get_prob_mean_maximaRR (long i) const {if (!calc_maximaRR) return 0.0; return local_maximaRR.get_prob_mean (i);}
    double get_prob_disp_minima (long i) const {if (!calc_minima) return 0.0; return local_minima.get_prob_disp (i);}
    double get_prob_disp_maxima (long i) const {if (!calc_maxima) return 0.0; return local_maxima.get_prob_disp (i);}
    double get_prob_disp_minimaRR (long i) const {if (!calc_minimaRR) return 0.0; return local_minimaRR.get_prob_disp (i);}
    double get_prob_disp_maximaRR (long i) const {if (!calc_maximaRR) return 0.0; return local_maximaRR.get_prob_disp (i);}
    double get_entropy_minima () const {if (!calc_minima) return 0.0; return local_minima.get_entropy ();}
    double get_entropy_maxima () const {if (!calc_maxima) return 0.0; return local_maxima.get_entropy ();}
    double get_entropy_minimaRR () const {if (!calc_minimaRR) return 0.0; return local_minimaRR.get_entropy ();}
    double get_entropy_maximaRR () const {if (!calc_maximaRR) return 0.0; return local_maximaRR.get_entropy ();}
    double get_entropy2_minima () const {if (!calc_minima) return 0.0; return local_minima.get_entropy2 ();}
    double get_entropy2_maxima () const {if (!calc_maxima) return 0.0; return local_maxima.get_entropy2 ();}
    double get_entropy2_minimaRR () const {if (!calc_minimaRR) return 0.0; return local_minimaRR.get_entropy2 ();}
    double get_entropy2_maximaRR () const {if (!calc_maximaRR) return 0.0; return local_maximaRR.get_entropy2 ();}
    double extreme_part_inside_minima (double sigma_in) const {if (!calc_minima) return 0.0; return local_minima.extreme_part_inside (sigma_in);}
    double extreme_part_inside_maxima (double sigma_in) const {if (!calc_maxima) return 0.0; return local_maxima.extreme_part_inside (sigma_in);}
    double extreme_part_inside_minimaRR (double sigma_in) const {if (!calc_minimaRR) return 0.0; return local_minimaRR.extreme_part_inside (sigma_in);}
    double extreme_part_inside_maximaRR (double sigma_in) const {if (!calc_maximaRR) return 0.0; return local_maximaRR.extreme_part_inside (sigma_in);}
    double extreme_part_outside_minima (double sigma_in) const {if (!calc_minima) return 0.0; return local_minima.extreme_part_outside (sigma_in);}
    double extreme_part_outside_maxima (double sigma_in) const {if (!calc_maxima) return 0.0; return local_maxima.extreme_part_outside (sigma_in);}
    double extreme_part_outside_minimaRR (double sigma_in) const {if (!calc_minimaRR) return 0.0; return local_minimaRR.extreme_part_outside (sigma_in);}
    double extreme_part_outside_maximaRR (double sigma_in) const {if (!calc_maximaRR) return 0.0; return local_maximaRR.extreme_part_outside (sigma_in);}
};

#endif