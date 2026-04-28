#ifndef MAXIMARR_H
#define MAXIMARR_H

#include "extrema.h"

class maximaRR : public extrema {
public:
    maximaRR () { limit_value = -1.0e12; }
    ~maximaRR() = default;
    void step (double new_value);
    void update_stat ();
    void get (bool &result, double &time, double &value);
};

void maximaRR::step (double new_value) {
    if (initialized) {
        long max_size = val_raw.size();
        marker_raw = (marker_raw + 1) % max_size;
        discrete_time_raw++;
        val_raw[marker_raw] = new_value;
        if (discrete_time_raw>=max_size) {
            long half_delay = max_size / 2;
            long half_marker_raw = (marker_raw + half_delay) % max_size;
            bool maximum_is_found = true;
            for (long h=0;h<val_raw.size();h++)
                if (h!=half_marker_raw)
                    if (val_raw[h]>(val_raw[half_marker_raw]-1.0e-12))
                        maximum_is_found = false;
            if ( (maximum_is_found)&&(val_raw[half_marker_raw]>limit_value) ) {
                lasttimeisfound = true;
                marker = (marker + 1) % max_size;
                discrete_time++;
                val[marker] = discrete_time_raw - half_delay;
                if (discrete_time==0)
                    lasttimeisfound = false;
                else
                    update_stat();
            }
            else
                lasttimeisfound = false;
        }
    }
    else {
        throw std::invalid_argument("Object buffer isn't defined before use\n");
    }
}

void maximaRR::update_stat () {
    long marker_y = marker + 2;
    if (marker_y >= val.size()) marker_y -= val.size();
    long d_marker_y = marker_y - 1;
    if (d_marker_y < 0) d_marker_y += val.size();
    long d_marker = marker - 1;
    if (d_marker < 0) d_marker += val.size();
    bool buffer_is_full = false;
    if (discrete_time>val.size()) {
        buffer_is_full = true;
    }
    update_values(val[marker]-val[d_marker], val[marker_y]-val[d_marker_y], buffer_is_full);
}

void maximaRR::get (bool &result, double &time, double &value) {
    result = false;
    if (lasttimeisfound) {
        long d_marker = marker - 1;
        if (d_marker<0) d_marker += val.size();
        value = val[marker] - val[d_marker];
        time = discrete_time_raw - (val.size() / 2);
        result = true;
    }
}

#endif