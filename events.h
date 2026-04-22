//Убрать сложные методы в конец файла

#ifndef EVENTS_H
#define EVENTS_H

#include "extrema.h"

class events : public extrema {
public:
    events () = default;
    ~events() = default;
    void step (double new_value) {
        if (initialized) {
            long max_size = val_raw.size();
            marker_raw = (marker_raw + 1) % max_size;
            discrete_time_raw++;
            val_raw[marker_raw] = new_value;
            if (new_value>0.5) {
                lasttimeisfound = true;
                marker = (marker + 1) % max_size;
                discrete_time++;
                val[marker] = discrete_time_raw;
                if (discrete_time==0)
                    lasttimeisfound = false;
                else
                    update_stat();
            }
            else
                lasttimeisfound = false;
        }
        else {
            throw std::invalid_argument("Object events isn't defined before use\n");
        }
    }
    void update_stat () {
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
    void get (bool &result, double &time, double &value) {
        result = false;
        if (lasttimeisfound) {
            long d_marker = marker - 1;
            if (d_marker<0) d_marker += val.size();
            value = val[marker] - val[d_marker];
            time = discrete_time_raw - (val.size() / 2);
            result = true;
        }
    }
};

#endif