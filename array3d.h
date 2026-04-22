//Убрать сложные методы в конец файла

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
    double get (long i) { return val[i]; }
    double get (long i1, long i2, long i3) { return val[i1*s2s3 + i2*size3 + i3]; }
    void set (long i, double value) { val[i] = value; }
    void set (long i1, long i2, long i3, double value) { val[i1*s2s3 + i2*size3 + i3] = value; }
};

#endif