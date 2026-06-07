#ifndef FHN_ENS_H
#define FHN_ENS_H

#include <omp.h>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

struct FHNState {
    double x, y;
};

struct FHNDeriv {
    double dx, dy;
};

struct Event {
    double t;
    int neuron;
};

struct FHNParams {
    double R = 0.1;
    double I = 0.5;
    double alpha = 0.7;
    double beta = 0.8;
    double epsilon = 0.08;
    double gauss_y_s = 0.0;
    double poisson_x_f = 1.0e-3;
    double poisson_x_a = 0.0;
};

inline FHNDeriv rhs (const FHNState& s, const FHNParams& p, const FHNState& in) {
    return {
        (s.x * (1.0 - s.x * s.x / 3.0) - s.y + p.R * p.I + in.x) / p.epsilon,
        s.x + p.alpha - p.beta * s.y + in.y
    };
}

class FHNbuffer {
    int size;
    int marker;
    bool filled;
    double th_x;
    double refractory;
    double lastspiketime;
    std::vector<FHNState> h;
    std::vector<double> time;
public:
    FHNbuffer () {
        size = 20;
        marker = size-1;
        filled = false;
        th_x = 0.5;
        refractory = 1.0;
        lastspiketime = -1.0e3;
        h.resize(size);
        time.resize(size);
    }
    void push (const FHNState& val_in, double t_in) {
        if (marker == 0) filled = true;
        marker = (marker + 1) % size;
        h[marker] = val_in;
        time[marker] = t_in;
    }
    void correct (const FHNState& val_in) {
        h[marker] = val_in;
    }
    bool detect_max_x (double &x_out, double &t_out) {
        int n = filled ? size : marker+1;
        if (n < size) return false;
        int check = (marker + size/2) % size;
        if (h[check].x<th_x) return false;
        int check_left = (check - 1 + size) % size;
        int check_right = (check + 1) % size;
        double dxleft = h[check].x - h[check_left].x;
        double dxright = h[check].x - h[check_right].x;
        if ((dxleft<0.0) || (dxright<0.0)) return false;
        for (int shift = 1-size/2; shift<(size/2); shift++) {
            int compared = (check + shift + size) % size;
            if (shift<0) {
                if (fabs(h[compared].x-h[check].x)<dxleft) return false;
            }
            if (shift>0) {
                if (fabs(h[compared].x-h[check].x)<dxright) return false;
            }
        }
        if ((time[check]-lastspiketime) < refractory) return false;
        x_out = h[check].x;
        t_out = time[check];
        lastspiketime = time[check];
        return true;
    }
};

class FHNens {
public:
    long N1 = 0, N2 = 0, M = 0, TN;
    std::vector<FHNState> tmp;
    std::vector<FHNState> state;
    std::vector<FHNParams> pars;
    std::vector<FHNbuffer> hist;
    std::vector<long> n_spikes;
    std::vector<long> n_input;
    std::vector<double> tx_event;
    std::vector<double> dW;
    std::vector<int> coup_sigma_P_x;
    std::vector<double> coup_sigma_x;
    std::vector<double> coup_gamma_x;
    std::vector<Event> poisson_events;
    std::vector<Event> layer1_events;
    std::vector<Event> layer2_events;
    FHNens () { N1 = 0; N2 = 0; M = 0; TN = 0; }
    void resize (int newN1, int newN2, int newM) {
        N1 = newN1; N2 = newN2; M = newM;
        TN = N1 * N2 * M;
        tmp.resize(TN);
        state.resize(TN);
        pars.resize(TN);
        hist.resize(TN);
        tx_event.resize(TN);
        dW.resize(TN);
        n_spikes.assign(TN, 0);
        n_input.assign(TN, 0);
        coup_sigma_x.resize(M);
        coup_sigma_P_x.resize(M);
        coup_gamma_x.resize(M*M);
        poisson_events.clear();
        layer1_events.clear();
        layer2_events.clear();
    }
    long d3to1 (long n1, long n2, long m) {
        return (n1 * N2 + n2) * M + m;
    }
    long d2to1 (long m1, long m2) {
        return m1 * M + m2;
    }
    void randomIC (double rx, double ry) {
        std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<double> dist(0.0, 2.0 * M_PI);
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    long idx = d3to1 (n1, n2, m);
                    double theta = dist(gen);
                    state[idx].x = rx * std::cos(theta);
                    state[idx].y = ry * std::sin(theta);
                    hist[idx].push(state[idx], 0.0);
                }
            }
        }
    }
    void equilibriumIC () {
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    long idx = d3to1 (n1, n2, m);
                    state[idx].x = -pars[idx].alpha;
                    state[idx].y = pars[idx].alpha * (pars[idx].alpha*pars[idx].alpha/3.0 - 1.0);
                    hist[idx].push(state[idx], 0.0);
                }
            }
        }
    }
    void setFHNpars (double new_alpha, double new_beta, double new_epsilon, double new_I, double new_R) {
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    long idx = d3to1 (n1, n2, m);
                    pars[idx].alpha = new_alpha;
                    pars[idx].beta = new_beta;
                    pars[idx].epsilon = new_epsilon;
                    pars[idx].I = new_I;
                    pars[idx].R = new_R;
                }
            }
        }
    }
    void set_poisson_pars (double new_f, double new_a) {
        std::mt19937 gen(std::random_device{}());
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    long idx = d3to1 (n1, n2, m);
                    pars[idx].poisson_x_f = new_f;
                    pars[idx].poisson_x_a = new_a;
                    auto &p = pars[idx];
                    std::exponential_distribution<double> exp(p.poisson_x_f);
                    tx_event[idx] = exp(gen);
                }
            }
        }
    }
    void set_poisson_pars (double new_f, double new_a, long layernum) {
        std::mt19937 gen(std::random_device{}());
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                int m = layernum; {
                    long idx = d3to1 (n1, n2, m);
                    pars[idx].poisson_x_f = new_f;
                    pars[idx].poisson_x_a = new_a;
                    auto &p = pars[idx];
                    std::exponential_distribution<double> exp(p.poisson_x_f);
                    tx_event[idx] = exp(gen);
                }
            }
        }
    }
    void set_coup_sigma_x (double new_sigma) {
        for (int m = 0; m < M; m++) set_coup_sigma_x (new_sigma, m);
    }
    void set_coup_sigma_x (double new_sigma, int m) {
        coup_sigma_x[m] = new_sigma;
    }
    void set_coup_sigma_P_x (int new_P) {
        for (int m = 0; m < M; m++) set_coup_sigma_P_x (new_P, m);
    }
    void set_coup_sigma_P_x (int new_P, int m) {
        coup_sigma_P_x[m] = new_P;
    }
    void set_coup_gamma_x (double new_gamma) {
        for (int m1 = 0; m1 < M; m1++) for (int m2 = 0; m2 < M; m2++) set_coup_gamma_x (new_gamma, m1, m2);
    }
    void set_coup_gamma_x (double new_gamma, int m1, int m2) {
        coup_gamma_x[d2to1(m1, m2)] = new_gamma;
    }
    double coupling_x (int n1, int m, const std::vector<FHNState>& S) {
        double cs = 0.0;
        if (N2==1) {
            for (int n1s = n1 - coup_sigma_P_x[m]; n1s <= n1 + coup_sigma_P_x[m]; n1s++) {
                int j = (n1s + N1) % N1;
                if (j != n1) cs += coup_sigma_x[m] * (S[d3to1 (j, 0, m)].x - S[d3to1 (n1, 0, m)].x);
            }
        }
        if (coup_sigma_P_x[m] > 0)
            cs /= (2.0*coup_sigma_P_x[m]);
        else
            cs = 0.0;
        double cg = 0.0;
        for (int m2 = 0; m2 < M; m2++)
            cg += coup_gamma_x[d2to1(m, m2)] * (S[d3to1 (n1, 0, m2)].x - S[d3to1 (n1, 0, m)].x);
        return cs + cg;
    }
    void events (double ct) {
        static thread_local std::mt19937 gen(std::random_device{}());
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            long idx = d3to1 (n1, n2, m);
            auto &p = pars[idx];
            if (p.poisson_x_f > 0.0) {
                std::exponential_distribution<double> exp(p.poisson_x_f);
                while (tx_event[idx] <= ct) {
                    poisson_events.push_back ( {
                        tx_event[idx],
                        n1
                    } );
                    state[idx].x += p.poisson_x_a;
                    hist[idx].correct(state[idx]);
                    n_input[idx]++;
                    tx_event[idx] += exp(gen);
                }
            }
        }
    }
    void step (double dt, double ct) {
        static thread_local std::mt19937 gen(std::random_device{}());
        static thread_local std::normal_distribution<double> normal(0.0, 1.0);
        events (ct);
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            long idx = d3to1 (n1, n2, m);
            auto &s = state[idx];
            auto &p = pars[idx];
            double cx = coupling_x(n1, m, state);
            dW[idx] = std::sqrt(dt) * normal(gen);
            FHNState in {cx, 0.0};
            auto f0 = rhs(s, p, in);
            tmp[idx].x = s.x + dt * f0.dx;
            tmp[idx].y = s.y + dt * f0.dy + p.gauss_y_s * dW[idx];
        }
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            long idx = d3to1 (n1, n2, m);
            auto &s = state[idx];
            auto &p = pars[idx];
            double cx0 = coupling_x(n1, m, state);
            double cx1 = coupling_x(n1, m, tmp);
            FHNState in0 {cx0, 0.0};
            FHNState in1 {cx1, 0.0};
            auto f0 = rhs(s, p, in0);
            auto f1 = rhs(tmp[idx], p, in1);
            s.x += 0.5 * dt * (f0.dx + f1.dx);
            s.y += 0.5 * dt * (f0.dy + f1.dy) + p.gauss_y_s * dW[idx];
            hist[idx].push(s, ct + dt);
        }
        update_spikes ();
    }
    void update_spikes () {
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            double x;
            double t;
            long idx = d3to1 (n1, n2, m);
            if (hist[idx].detect_max_x(x, t)) {
                n_spikes[idx]++;
                if (m == 0) {
                    layer1_events.push_back ( {
                        t,
                        n1
                    } );
                }
                if (m == 1) {
                    layer2_events.push_back ( {
                        t,
                        n1
                    } );
                }
            }
        }
    }
};

void runnet (double sigmaloc, double gammaloc) {
    FHNens net;
    net.resize (65, 1, 2);
    net.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
    net.set_coup_sigma_x (0.0, 0);
    net.set_coup_sigma_P_x (0, 0);
    net.set_coup_sigma_x (sigmaloc, 1);
    net.set_coup_sigma_P_x (10, 1);
    net.set_coup_gamma_x (0.0);
    net.set_coup_gamma_x (gammaloc, 1, 0);
    net.equilibriumIC ();
    net.set_poisson_pars (0.0, 0.0);
    net.set_poisson_pars (0.6, 1.0, 0);
    const double dt = 0.0001;
    const double T  = 10000.0;
    for (double t = 0.0; t < T; t += dt)
        net.step(dt, t);
    std::ostringstream prefix;
    prefix << std::fixed << std::setprecision(4)
           << "s" << sigmaloc << "-g" << gammaloc;
    std::ofstream fout0(prefix.str() + "-input.dat");
    for (const auto &e : net.poisson_events)
        fout0 << e.t << " " << e.neuron << "\n";
    std::ofstream fout1(prefix.str() + "-layer1.dat");
    for (const auto &e : net.layer1_events)
        fout1 << e.t << " " << e.neuron << "\n";
    std::ofstream fout2(prefix.str() + "-layer2.dat");
    for (const auto &e : net.layer2_events)
        fout2 << e.t << " " << e.neuron << "\n";
}

struct RunPars {
    double sigma;
    double gamma;
};

void run () {
    std::vector<RunPars> jobs;
    for (double s = 0.00; s < 1.00; s += 0.05)
        for (double g = 0.00; g < 1.00; g += 0.05)
            jobs.push_back({s, g});
    #pragma omp parallel for schedule(dynamic)
    for (long i = 0; i < (long)jobs.size(); i++) {
        runnet(
            jobs[i].sigma,
            jobs[i].gamma
        );
    }
}

#endif
