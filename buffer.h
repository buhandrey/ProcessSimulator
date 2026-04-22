#ifndef BUFFER_H
#define BUFFER_H

#include <iostream>
#include <sstream>
#include <string>

#include "process.h"

class buffer : public process {
public:
    buffer() {
        modelname = "buffer";
    };
    ~buffer() = default;
    std::string get_basename () {
        std::stringstream name_temp;
        name_temp << "_" << modelname;
        return name_temp.str();
    }
    void correct (long delay, double new_val) {
        long d_marker = marker - delay;
        if (d_marker<0) d_marker+=val.size();
        val[d_marker] = new_val;
    }
    void correct (double new_val) {
        val[marker] = new_val;
    }
    void step () {
        throw std::invalid_argument("Object buffer isn't designed to use empty step()\n");
    }
    void step (double new_value) {
        if (initialized) {
            marker = (marker + 1) % val.size();
            discrete_time++;
            val[marker] = new_value;
            substep();
        }
        else {
            throw std::invalid_argument("Object buffer isn't defined before use\n");
        }
    }
    void initwithval (double new_value) {
        for (long h=0;h<val.size();h++)
            val[h] = new_value;
        step(new_value);
    }
};

#endif