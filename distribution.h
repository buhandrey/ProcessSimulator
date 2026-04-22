#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

class distribution {
protected:
    std::vector<long> prob_histogram;
    std::vector<double> prob_mins;
    std::vector<double> prob_maxs;
    std::vector<double> prob_M1;
    std::vector<double> prob_M2;
    long distr_size;
    long prob_filled;
    long stat_counter;
    double global_width;
    double global_stat_M1_x;
    double global_stat_M2_x;
public:
    distribution () {
        distr_size = 201;
        prob_filled = 0;
        stat_counter = 0;
        global_width = 0.0;
        global_stat_M1_x = 0.0;
        global_stat_M2_x = 0.0;
        resize_arrays ();
    }
    void resize (long distr_size_manual) {
        if (distr_size_manual>1)
            distr_size = distr_size_manual;
        else
            distr_size = 201;
        resize_arrays ();
    }
    ~distribution() = default;
    int calc_probability (double x) {
        const double n = static_cast<double>((stat_counter+1));
        double delta_x = x - global_stat_M1_x;
        double delta_x_n = delta_x/n;
        double term1_x = delta_x * delta_x_n * (n - 1);
        global_stat_M2_x += term1_x;
        global_stat_M1_x += delta_x_n;
        if (prob_filled==0) {
            init_interval(x, prob_filled);
            return 0;
        }
        for (long i=0;i<prob_filled;i++) {
            if ( (x>=(prob_mins[i]-1.0e-6)) && (x<=(prob_maxs[i]+1.0e-6)) ) {
                add_val_to_interval (x, i);
                return 0;
            }
        }
        if ((prob_filled>20)) {
            for (long i=0;i<prob_filled;i++) {
                bool x_is_not_lapped = false;
                if ((i==0) && (x<(prob_M1[i]))) x_is_not_lapped = true;
                if ((i==(prob_filled-1)) && (x>(prob_M1[i]))) x_is_not_lapped = true;
                if ((i>0) && (i<(prob_filled-1)) && ((x>(prob_M1[i-1]))&&(x<(prob_M1[i+1])))) x_is_not_lapped = true;
                if (x_is_not_lapped) {
                    if (fabs(x-prob_M1[i])<(0.5*global_width/prob_filled)) {
                        add_val_to_interval (x, i);
                        return 0;
                    }
                }
            }
        }
        bool interval_added = add_interval(x);
        if (interval_added) return 0;
        long nearest_i = -1;
        double min_dx = 1.0e20;
        long nearest_h_left_i = -1;
        double min_dh = 1.0e20;
        bool x_at_left = false;
        for (long i=0;(i<prob_filled)&&(!x_at_left);i++) {
            if (i>0) {
                if ((x<prob_M1[i])&&(x>prob_M1[i-1])) x_at_left = true;
            }
            else {
                if (x<prob_M1[i]) x_at_left = true;
            }
            if (x_at_left) {
                double dx = prob_M1[i]-x;
                double min_dx_test = (prob_histogram[i]+1.0)*dx*dx;
                if (min_dx_test<min_dx) {
                    min_dx = min_dx_test;
                    nearest_i = i;
                }
            }
        }
        bool x_at_right = false;
        for (long i=0;(i<prob_filled)&&(!x_at_right);i++) {
            if (i<(prob_filled-1)) {
                if ((x>prob_M1[i])&&(x<prob_M1[i+1])) x_at_right = true;
            }
            else {
                if (x>prob_M1[i]) x_at_right = true;
            }
            if (x_at_right) {
                double dx = x-prob_M1[i];
                double min_dx_test = (prob_histogram[i]+1.0)*dx*dx;
                if (min_dx_test<min_dx) {
                    min_dx = min_dx_test;
                    nearest_i = i;
                }
            }
        }
        for (long i=0;(i<(prob_filled-1));i++) {
            double n_i = prob_histogram[i];
            double n_j = prob_histogram[i+1];
            double dmu = prob_M1[i+1] - prob_M1[i];
            double min_dh_test = (n_i * n_j / (n_i + n_j)) * dmu * dmu;
            if (min_dh_test<min_dh) {
                min_dh = min_dh_test;
                nearest_h_left_i = i;
            }
        }
        bool bad_intervals = false;
        if (nearest_h_left_i>=0) {
            if ((x>prob_mins[nearest_h_left_i]) && (x<prob_maxs[nearest_h_left_i+1]))
                bad_intervals = true;
        }
        if ((min_dx<min_dh)||(bad_intervals)) {
            add_val_to_interval (x, nearest_i);
            return 0;
        }
        if (min_dh<1.0e10) {
            merge_and_init(x, nearest_h_left_i);
            return 0;
        }
        return 0;
    }
    double get_prob_dens (long i) const {
        if ((i<distr_size)&&(i>=0)) {
            double width = prob_maxs[i] - prob_mins[i];
            if ((width <= 0.0)||(prob_filled == 0)) return 0;
            return (prob_histogram[i])/(stat_counter*width);
        }
        else
            return 0.0;
    }
    double get_prob_left (long i) const {
        if ((i<distr_size)&&(i>=0)) {
            return prob_mins[i];
        }
        else
            return 0.0;
    }
    double get_prob_right (long i) const {
        if ((i<distr_size)&&(i>=0)) {
            return prob_maxs[i];
        }
        else
            return 0.0;
    }
    double get_prob_mean (long i) const {
        if ((i<distr_size)&&(i>=0)) {
            return prob_M1[i];
        }
        else
            return 0.0;
    }
    double get_prob_disp (long i) const {
        if ((i<distr_size)&&(i>=0)) {
            return prob_M2[i]/(1.0*(prob_histogram[i]-1.0));
        }
        else
            return 0.0;
    }
    double get_entropy () const {
        double entropy_sum = 0.0;
        for (long i=0;i<prob_filled;i++) {
            double width = prob_maxs[i] - prob_mins[i];
            if (width <= 0.0) continue;
            double p_i = (1.0*prob_histogram[i])/(1.0*stat_counter);
            entropy_sum -= p_i * log(p_i) - p_i*log(width);
        }
        return entropy_sum;
    }
    double extreme_part_inside (double sigma_in) const {
        if (stat_counter<=1) return NAN;
        double sigma = sqrt(global_stat_M2_x/(stat_counter-1.0));
        double left = global_stat_M1_x - sigma_in * sigma;
        double right = global_stat_M1_x + sigma_in * sigma;
        double prob_inside = 0.0;
        for (long i=0;i<prob_filled;i++) {
            double bin_l = prob_mins[i];
            double bin_r = prob_maxs[i];
            double width = bin_r - bin_l; if (width <= 0.0) continue;
            double p_i = (1.0 * prob_histogram[i])/(1.0*stat_counter);
            if (bin_r <= left || bin_l >= right) continue;
            if (bin_l >= left && bin_r <= right) {
                prob_inside += p_i; continue;
            }
            double overlap_l = (bin_l > left) ? bin_l : left;
            double overlap_r = (bin_r < right) ? bin_r : right;
            double overlap = overlap_r - overlap_l;
            if (overlap > 0.0)
                prob_inside += p_i * (overlap / width);
        }
        return prob_inside;
    }
    double extreme_part_outside (double sigma_out) const {
        return 1.0 - extreme_part_inside(sigma_out);
    }
    long get_bins_number () {
        return prob_filled;
    }
protected:
    void resize_arrays () {
        prob_histogram.resize(distr_size, 0);
        prob_mins.resize(distr_size, .0);
        prob_maxs.resize(distr_size, .0);
        prob_M1.resize(distr_size, .0);
        prob_M2.resize(distr_size, .0);
    }
    void shift_arrays (long from, long to) {
        prob_histogram[to] = prob_histogram[from];
        prob_mins[to] = prob_mins[from];
        prob_maxs[to] = prob_maxs[from];
        prob_M1[to] = prob_M1[from];
        prob_M2[to] = prob_M2[from];
    }
    void init_interval (double x, long i) {
        reinit_interval (x, i);
        prob_filled++;
    }
    void reinit_interval (double x, long i) {
        global_width += 2.0e-6;
        prob_histogram[i] = 1;
        prob_mins[i] = x-1.0e-6;
        prob_maxs[i] = x+1.0e-6;
        prob_M1[i] = x;
        prob_M2[i] = 0.0;
        stat_counter++;
    }
    bool add_interval (double x) {
        if (prob_filled<distr_size) {
            long i_right = prob_filled;
            for (long i=(prob_filled-1); i>=0; i--) {
                if (x<(prob_mins[i]-1.0e-6)) i_right = i;
            }
            if (i_right>0) if (x<(prob_maxs[i_right-1]+1.0e-6)) return false;
            for (long i=prob_filled; i>i_right; i--) {
                shift_arrays(i-1, i);
            }
            init_interval(x, i_right);
            return true;
        }
        return false;
    }
    void add_val_to_interval (double x, long i) {
        prob_histogram[i]++;
        if (x<prob_mins[i]) {
            global_width += prob_mins[i]-x;
            prob_mins[i] = x;
        }
        if (x>prob_maxs[i]) {
            global_width += x-prob_maxs[i];
            prob_maxs[i] = x;
        }
        double old_mid = prob_M1[i];
        prob_M1[i] += (x - old_mid) / (1.0*prob_histogram[i]);
        prob_M2[i] += (x - old_mid) * (x - prob_M1[i]);
        stat_counter++;
    }
    void merge_and_init (double x, long nearest_h_left_i) {
        global_width += prob_mins[nearest_h_left_i+1] - prob_maxs[nearest_h_left_i];
        double n_new = prob_histogram[nearest_h_left_i] + prob_histogram[nearest_h_left_i+1];
        double n_prod = prob_histogram[nearest_h_left_i] * prob_histogram[nearest_h_left_i+1];
        double delta = prob_M1[nearest_h_left_i+1]-prob_M1[nearest_h_left_i];
        double disp_new = prob_M2[nearest_h_left_i+1] + prob_M2[nearest_h_left_i] + n_prod*delta*delta/n_new;
        double mean_new = prob_M1[nearest_h_left_i] + delta*prob_histogram[nearest_h_left_i+1]/n_new;
        long replacing_i = nearest_h_left_i;
        if (x<prob_M1[nearest_h_left_i]) {
            prob_M2[nearest_h_left_i+1] = disp_new;
            prob_M1[nearest_h_left_i+1] = mean_new;
            prob_histogram[nearest_h_left_i+1]+=prob_histogram[nearest_h_left_i];
            prob_mins[nearest_h_left_i+1] = prob_mins[nearest_h_left_i];
            for (;(replacing_i>0)&&(prob_M1[replacing_i-1]>x);replacing_i--) {
                shift_arrays(replacing_i-1, replacing_i);
            }
        }
        else {
            prob_M2[nearest_h_left_i] = disp_new;
            prob_M1[nearest_h_left_i] = mean_new;
            prob_histogram[nearest_h_left_i]+=prob_histogram[nearest_h_left_i+1];
            prob_maxs[nearest_h_left_i] = prob_maxs[nearest_h_left_i+1];
            for (replacing_i = (nearest_h_left_i+1); (replacing_i<(prob_filled-1))&&(prob_M1[replacing_i+1]<x); replacing_i++) {
                shift_arrays(replacing_i+1, replacing_i);
            }
        }
        reinit_interval(x, replacing_i);
        stat_counter++;
    }
};

#endif