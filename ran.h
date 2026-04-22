#ifndef RAN_H
#define RAN_H

#include "noise.h"

class ran : public noise {
protected:
    double par_mu; /*Параметр среднего значения*/
    double par_range; /*Параметр диапазона значений*/
    /*Статические параметры генератора чисел*/
    long MBIG, MSEED, MZ, ma[56];
    double FAC;
    int inext, inextp, iff;
    /*Статические параметры генератора чисел*/
    double generator () {
        long mj,mk;
        int i,ii,k;
        if (idum < 0 || iff == 0) {
            iff=1;
            mj=labs(MSEED-labs(idum));
            mj %= MBIG;
            ma[55]=mj;
            mk=1;
            for (i=1;i<=54;i++) {
                ii=(21*i) % 55;
                ma[ii]=mk;
                mk=mj-mk;
                if (mk < MZ) mk += MBIG;
                mj=ma[ii];
            }
            for (k=1;k<=4;k++)
                for (i=1;i<=55;i++) {
                    ma[i] -= ma[1+(i+30) % 55];
                    if (ma[i] < MZ) ma[i] += MBIG;
                }
                inext=0;
            inextp=31;
            idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return par_mu + par_range*(mj*FAC-0.5);
    }
public:
    ran () {
        modelname = "ran_noise";
    }
    ~ran() = default;
    std::string get_basename () {
        std::stringstream name_temp;
        name_temp << "_" << modelname << "_range_" << par_range << "_m_" << par_mu << "_idum_" << idum;
        return name_temp.str();
    }
    void set_pars (int set_size, double set_par_mu, double set_par_range, long set_idum) {
        if (set_size>0) {
            val.resize(set_size);
        }
        else {
            throw std::invalid_argument("Vector size must be positive (from ran)\n");
        }
        par_mu = set_par_mu;
        if (set_par_range>.0) {
            par_range = set_par_range;
        }
        else {
            throw std::invalid_argument("Noise range must be positive (from ran)\n");
        }
        if (set_idum<0) {
            idum = set_idum;
        }
        else {
            throw std::invalid_argument("Seed idum must be positive (from ran)\n");
        }
        /*Инициализация генератора чисел*/
        MBIG = 1000000000;
        MSEED = 161803398;
        MZ = 0;
        FAC = (1.0/MBIG);
        inext = 0, inextp = 0;
        for (int i=0;i<56;i++) ma[i] = 0;
        iff = 0;
        /*Инициализация генератора чисел*/
        initialized = true;
    }
};

#endif