#ifndef EXTREMA_H
#define EXTREMA_H

#include <iostream>
#include <vector>

#include "observers.h"

class extrema : public observers {
protected:
    long marker;
    long marker_raw;
    long discrete_time;
    long discrete_time_raw;
    bool initialized;
    bool lasttimeisfound;
    double limit_value;
    std::vector<double> val;
    std::vector<double> val_raw;
public:
    extrema () {
        marker = -1;
        marker_raw = -1;
        discrete_time = -1;
        discrete_time_raw = -1;
        initialized = false;
        lasttimeisfound = false;
    }
    ~extrema() = default;
    void step () {
        throw std::invalid_argument("Object extrema isn't designed to use empty step()\n");
    }
    long size() {
        return val.size();
    }
    void resize (long set_size) {
        val.resize(set_size, 0.0);
        val_raw.resize(set_size, 0.0);
    }
    void init () {
        if (val.size()>1) {
            initialized = true;
            corr_init();
            stat_init();
            prob_init();
        }
        else
            throw std::invalid_argument("The size of val should be not less than 2 (from maximaRR)\n");
    }
    void set_limit_value (double limit_value_new) {
        limit_value = limit_value_new;
    }

};

#endif