#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <iostream>
#include <vector>

template <typename array_element>
class array3d {
protected:
    long size1;
    long size2;
    long size3;
    long s2s3;
    long s1s2s3;
    std::vector<array_element> val;
public:
    array3d () {
        size1 = 0;
        size2 = 0;
        size3 = 0;
        s2s3 = 0;
        s1s2s3 = 0;
    }
    ~array3d() = default;
    void resize(long set_size1, long set_size2, long set_size3) {
        size1 = set_size1;
        size2 = set_size2;
        size3 = set_size3;
        s2s3 = size2 * size3;
        s1s2s3 = size1 * size2 * size3;
        val.resize(s1s2s3, 0.0);
    }
    long s1 () { return size1; }
    long s2 () { return size2; }
    long s3 () { return size3; }
    long s23 () { return s2s3; }
    long s123 () { return s1s2s3; }
    double get (long i) {
        if (i>=s1s2s3)
            throw std::invalid_argument("Out of range 1 (from array3d)\n");
        return val[i];
    }
    double get (long i1, long i2, long i3) {
        if (i1>=size1 || i2>=size2 || i3>=size3)
            throw std::invalid_argument("Out of range 2 (from array3d)\n");
        return val[i1*s2s3 + i2*size3 + i3];
    }
    void set (long i, double value) {
        if (i>=s1s2s3)
            throw std::invalid_argument("Out of range 3 (from array3d)\n");
        val[i] = value;
    }
    void set (long i1, long i2, long i3, double value) {
        if (i1>=size1 || i2>=size2 || i3>=size3)
            throw std::invalid_argument("Out of range 4 (from array3d)\n");
        val[i1*s2s3 + i2*size3 + i3] = value;
    }
};

#endif