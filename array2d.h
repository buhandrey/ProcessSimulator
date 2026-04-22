#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <iostream>
#include <vector>

template <typename array_element>
class array2d {
protected:
    long size0;
    long size1;
    long size2;
    long s1s2;
    std::vector<array_element> val;
public:
    array2d () {
        size0 = 0;
        size1 = 0;
        size2 = 0;
        s1s2 = 0;
    }
    ~array2d() = default;
    void resize(long set_size0, long set_size1, long set_size2) {
        size0 = set_size0;
        size1 = set_size1;
        size2 = set_size2;
        s1s2 = size1 * size2;
        val.resize(s1s2);
        for (long s=0;s<s1s2;s++) {
            val[s].resize(size0);
            val[s].init();
        }
    }
    long s0 () { return size0; }
    long s1 () { return size1; }
    long s2 () { return size2; }
    long s12 () { return s1s2; }
    double get (long i) {
        if (i>=s1s2)
            throw std::invalid_argument("Out of range 1 (from array2d)\n");
        return val[i].var();
    }
    double get_delayed (long delay, long i) {
        return val[i].var(delay);
    }
    double get (long i1, long i2) {
        return val[i1*size2 + i2].var();
    }
    double get_delayed (long delay, long i1, long i2) {
        return val[i1*size2 + i2].var(delay);
    }
    void initbuffer (long i, double value) {
        val[i].initwithval(value);
    }
    void initbuffer (long i1, long i2, double value) {
        val[i1*size2 + i2].initwithval(value);
    }
    void step (long i, double value) {
        val[i].step(value);
    }
    void step (long i1, long i2, double value) {
        val[i1*size2 + i2].step(value);
    }
    void correct (long i, double value) {
        val[i].correct(value);
    }
    void correct (long i1, long i2, double value) {
        val[i1*size2 + i2].correct(value);
    }
    void correct (long delay, long i1, long i2, double value) {
        val[i1*size2 + i2].correct(delay, value);
    }
};

#endif