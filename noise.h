#ifndef NOISE_H
#define NOISE_H

#include "process.h"

class noise : public process {
protected:
    long idum; /*Зерно шума*/
    virtual double generator () = 0;
public:
    noise() { modelname = "noise"; };
    ~noise() = default;
    void step ();
};

void noise::step () {
    if (initialized) {
        marker++;
        if (marker == val.size())
            marker = 0;
        discrete_time++;
        val[marker] = generator ();
        substep();
    }
    else {
        throw std::invalid_argument("Object noise isn't defined before use\n");
    }
}

#endif