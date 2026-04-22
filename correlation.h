#ifndef CORRELATION_H
#define CORRELATION_H

#include <cmath>

class correlation {
protected:
    long long stat_counter;
    double stat_M1_x;
    double stat_M2_x;
    double stat_M1_y;
    double stat_M2_y;
    double stat_M2_xy;
public:
    correlation () {
        stat_counter = 0;
        stat_M1_x = 0.0;
        stat_M2_x = 0.0;
        stat_M1_y = 0.0;
        stat_M2_y = 0.0;
        stat_M2_xy = 0.0;
    }
    ~correlation() = default;
    void update_corr (double x, double y) {
        stat_counter++;
        const double n = static_cast<double>(stat_counter);
        double delta_x = x - stat_M1_x;
        double delta_y = y - stat_M1_y;
        double delta_x_n = delta_x/n;
        double delta_y_n = delta_y/n;
        double term1_x = delta_x * delta_x_n * (n - 1);
        double term1_y = delta_y * delta_y_n * (n - 1);
        stat_M2_x += term1_x;
        stat_M2_y += term1_y;
        stat_M1_x += delta_x_n;
        stat_M1_y += delta_y_n;
        stat_M2_xy += delta_x * (y - stat_M1_y);
    }
    double stat_get_correlation () const {
        if ((stat_M2_x == 0) || (stat_M2_y == 0)) return NAN;
        return stat_M2_xy/sqrt(stat_M2_x*stat_M2_y);
    }
};

#endif 