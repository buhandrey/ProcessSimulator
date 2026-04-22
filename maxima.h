#ifndef MAXIMA_H
#define MAXIMA_H

#include "extrema.h"

class maxima : public extrema {
public:
    maxima () {
        limit_value = -1.0e12;
    }
    ~maxima() = default;
    void step (double new_value) {
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
                    val[marker] = val_raw[half_marker_raw];
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
    void get (bool &result, double &time, double &value) {
        result = false;
        if (lasttimeisfound) {
            value = val[marker];
            time = discrete_time_raw - (val.size() / 2);
            result = true;
        }
    }
};

#endif