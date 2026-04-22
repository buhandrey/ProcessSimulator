#ifndef STATISTICS_H
#define STATISTICS_H

#include <cmath>

class statistics {
protected:
    long long stat_counter;
    double stat_M1_x;
    double stat_M2_x;
    double stat_M3_x;
    double stat_M4_x;
public:
    statistics () {
        stat_counter = 0;
        stat_M1_x = 0.0;
        stat_M2_x = 0.0;
        stat_M3_x = 0.0;
        stat_M4_x = 0.0;
    }
    ~statistics() = default;
    void update_stat (double x) {
        stat_counter++;
        const double n = static_cast<double>(stat_counter);
        double delta_x = x - stat_M1_x;
        double delta_x_n = delta_x/n;
        double delta_x_n2 = delta_x_n * delta_x_n;
        double term1_x = delta_x * delta_x_n * (n - 1);
        if (stat_counter>=4)
            stat_M4_x += delta_x_n2 * term1_x * (n * n - 3.0 * n + 3.0) + 6.0 * delta_x_n2 * stat_M2_x - 4.0 * delta_x_n * stat_M3_x;
        if (stat_counter>=3)
            stat_M3_x += term1_x * delta_x_n * (n-2.0) - 3.0 * stat_M2_x * delta_x_n;
        stat_M2_x += term1_x;
        stat_M1_x += delta_x_n;
    }
    double stat_get_mean () const {
        return stat_M1_x;
    }
    double stat_get_disp_root () const {
        if (stat_counter<=1) return NAN;
        return sqrt(stat_M2_x/(stat_counter-1.0));
    }
    double stat_get_disp () const {
        if (stat_counter<=1) return NAN;
        return stat_M2_x/(stat_counter-1.0);
    }
    double stat_get_skew () const {
        if ((stat_counter<=2)||(stat_M2_x==0)) return NAN;
        double stat_M2_x_inv = 1.0/stat_M2_x;
        return sqrt(1.0*stat_counter*stat_M2_x_inv)*stat_M3_x*stat_M2_x_inv;
    }
    double stat_get_kurt () const {
        if ((stat_counter<=3)||(stat_M2_x==0)) return NAN;
        double stat_M2_x_inv = 1.0/stat_M2_x;
        return (1.0*stat_counter*stat_M2_x_inv)*(stat_M4_x*stat_M2_x_inv) - 3.0;
    }
};

#endif