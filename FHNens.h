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

struct FHNState {
    double x, y;
};

struct FHNDeriv {
    double dx, dy;
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
        n_spikes.resize(TN);
        n_input.resize(TN);
        coup_sigma_x.resize(M);
        coup_sigma_P_x.resize(M);
        coup_gamma_x.resize(M*M);
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
                    double theta = dist(gen);
                    state[d3to1 (n1, n2, m)].x = rx * std::cos(theta);
                    state[d3to1 (n1, n2, m)].y = ry * std::sin(theta);
                    hist[d3to1 (n1, n2, m)].push(state[d3to1 (n1, n2, m)], 0.0);
                }
            }
        }
    }
    void equilibriumIC () {
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    state[d3to1 (n1, n2, m)].x = -pars[d3to1 (n1, n2, m)].alpha;
                    state[d3to1 (n1, n2, m)].y = pars[d3to1 (n1, n2, m)].alpha * (pars[d3to1 (n1, n2, m)].alpha*pars[d3to1 (n1, n2, m)].alpha/3.0 - 1.0);
                    hist[d3to1 (n1, n2, m)].push(state[d3to1 (n1, n2, m)], 0.0);
                }
            }
        }
    }
    void setFHNpars (double new_alpha, double new_beta, double new_epsilon, double new_I, double new_R) {
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    pars[d3to1 (n1, n2, m)].alpha = new_alpha;
                    pars[d3to1 (n1, n2, m)].beta = new_beta;
                    pars[d3to1 (n1, n2, m)].epsilon = new_epsilon;
                    pars[d3to1 (n1, n2, m)].I = new_I;
                    pars[d3to1 (n1, n2, m)].R = new_R;
                }
            }
        }
    }
    void set_poisson_pars (double new_f, double new_a) {
        std::mt19937 gen(std::random_device{}());
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    pars[d3to1 (n1, n2, m)].poisson_x_f = new_f;
                    pars[d3to1 (n1, n2, m)].poisson_x_a = new_a;
                    auto &p = pars[d3to1 (n1, n2, m)];
                    std::exponential_distribution<double> exp(p.poisson_x_f);
                    tx_event[d3to1 (n1, n2, m)] = exp(gen);
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
        double cg = 0.0;
        for (int m2 = 0; m2 < M; m2++)
            cg += coup_gamma_x[d2to1(m, m2)] * (S[d3to1 (n1, 0, m2)].x - S[d3to1 (n1, 0, m)].x);
        return cs + cg;
    }
    void events (double ct) {
        static thread_local std::mt19937 gen(std::random_device{}());
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            auto &p = pars[d3to1 (n1, n2, m)];
            if (p.poisson_x_f > 0.0) {
                std::exponential_distribution<double> exp(p.poisson_x_f);
                while (tx_event[d3to1 (n1, n2, m)] <= ct) {
                    state[d3to1 (n1, n2, m)].x += p.poisson_x_a;
                    hist[d3to1 (n1, n2, m)].correct(state[d3to1 (n1, n2, m)]);
                    n_input[d3to1 (n1, n2, m)]++;
                    tx_event[d3to1 (n1, n2, m)] += exp(gen);
                }
            }
        }
    }
    void step (double dt, double ct) {
        static thread_local std::mt19937 gen(std::random_device{}());
        static thread_local std::normal_distribution<double> normal(0.0, 1.0);
        events (ct);
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            auto &s = state[d3to1 (n1, n2, m)];
            auto &p = pars[d3to1 (n1, n2, m)];
            double cx = coupling_x(n1, m, state);
            dW[d3to1 (n1, n2, m)] = std::sqrt(dt) * normal(gen);
            FHNState in {cx, 0.0};
            auto f0 = rhs(s, p, in);
            tmp[d3to1 (n1, n2, m)].x = s.x + dt * f0.dx;
            tmp[d3to1 (n1, n2, m)].y = s.y + dt * f0.dy + p.gauss_y_s * dW[d3to1 (n1, n2, m)];
        }
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            auto &s = state[d3to1 (n1, n2, m)];
            auto &p = pars[d3to1 (n1, n2, m)];
            double cx0 = coupling_x(n1, m, state);
            double cx1 = coupling_x(n1, m, tmp);
            FHNState in0 {cx0, 0.0};
            FHNState in1 {cx1, 0.0};
            auto f0 = rhs(s, p, in0);
            auto f1 = rhs(tmp[d3to1 (n1, n2, m)], p, in1);
            s.x += 0.5 * dt * (f0.dx + f1.dx);
            s.y += 0.5 * dt * (f0.dy + f1.dy) + p.gauss_y_s * dW[d3to1 (n1, n2, m)];
            hist[d3to1 (n1, n2, m)].push(s, ct + dt);
        }
        update_spikes ();
    }
    void update_spikes () {
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            double x;
            double t;
            if (hist[d3to1 (n1, n2, m)].detect_max_x(x, t)) {
                n_spikes[d3to1 (n1, n2, m)]++;
            }
        }
    }
    void collect_max_x (std::ostream &out) {
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            double x;
            double t;
            if (hist[d3to1 (n1, n2, m)].detect_max_x(x, t)) {
                out << t << " " << n1 << " " << n2 << " " << m << " " << x << "\n";
            }
        }
    }
};

void Ring100std () {
    FHNens ring;
    ring.resize (100, 1, 1);
    ring.set_coup_sigma_x (0.0);
    ring.set_coup_gamma_x (0.0);
    ring.set_coup_sigma_P_x (1);
    ring.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
    ring.equilibriumIC ();
    double dt = 0.0001;
    ring.set_poisson_pars (0.5, 1.0);
    long steps = 200000;
    long trans = 0;
    long save_every = 2000;
    std::ofstream fout("spacetime.dat");
    std::ofstream sout("spikes.dat");
    for (long step_i = 0; step_i < steps; step_i++) {
        ring.step(dt, step_i * dt);
        if ((step_i > trans) && (step_i % save_every == 0))
            for (int n1 = 0; n1 < ring.N1; n1++) fout << step_i * dt << " " << n1 << " " << ring.state[ring.d3to1 (n1, 0, 0)].x << "\n";
        ring.collect_max_x(sout);
    }
    sout.close();
    fout.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << 'EOF'\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20';";
    plot_command << " set key tmargin center horizontal;";
    plot_command << " set xlabel 'i'; set ylabel 't'; set cblabel 'x';";
    plot_command << " set xrange [0.5:100.5]; set yrange [:]; set cbrange [-2.5:2.5];";
    plot_command << " set xtics 20 out; set ytics 10 out; set cbtics 1;";
    plot_command << " set mxtics 5; set mytics 5; set mcbtics 5;";
    plot_command << " set output 'spacetime.png';";
    plot_command << " plot 'spacetime.dat' u (1.0+$2):1:3 w image notitle,";
    plot_command << " 'spikes.dat' u (1.0+$2):1 w p pt 7 ps 0.9 lc rgb 'white' notitle;\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

void Ring100FoutOnFin () {
    long nofth = 200;
    long N1 = 100;
    double dt = 0.0001;
    long steps = 10000*5000;
    std::vector<double> nin;
    std::vector<double> nout;
    std::vector<double> freq;
    nin.resize(nofth);
    nout.resize(nofth);
    freq.resize(nofth);
    #pragma omp parallel for schedule(dynamic)
    for (long iter = 0; iter < nofth; iter++) {
        FHNens ring;
        freq[iter] = pow(10.0, iter*0.025-3.0);
        ring.resize (N1, 1, 1);
        ring.set_coup_sigma_x (0.0);
        ring.set_coup_gamma_x (0.0);
        ring.set_coup_sigma_P_x (1);
        ring.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
        ring.equilibriumIC ();
        ring.set_poisson_pars (freq[iter], 1.0);
        for (long step_i = 0; step_i < steps; step_i++) {
            ring.step(dt, step_i * dt);
        }
        nin[iter] = 0.0;
        nout[iter] = 0.0;
        for (int n1 = 0; n1 < N1; n1++) {
            nin[iter] += ring.n_input[ring.d3to1 (n1, 0, 0)];
            nout[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 0)];
        }
    }
    std::ofstream fout("FoutOnFin.dat");
    for (long iter = 0; iter < nofth; iter++) {
        fout << freq[iter] << " " << nin[iter]/(steps*dt*N1) << " " << nout[iter]/(steps*dt*N1) << "\n";
    }
    fout.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << 'EOF'\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20';";
    plot_command << " set key tmargin center horizontal;";
    plot_command << " set xlabel 'f(Poisson)'; set ylabel 'f(sp)';";
    plot_command << " set log; set xrange [1.0e-3:1.0]; set yrange [:];";
    plot_command << " set output 'FoutOnFin.png';";
    plot_command << " plot 'FoutOnFin.dat' u 1:3 w l lw 3 notitle;\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

void Ring100 () {
    //Ring100std ();
    Ring100FoutOnFin ();
}

#endif