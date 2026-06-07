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
    void set_poisson_pars (double new_f, double new_a, long layernum) {
        std::mt19937 gen(std::random_device{}());
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                int m = layernum; {
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
        return cs / (2.0*coup_sigma_P_x[m]) + cg;
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

// void Ring100std () {
//     FHNens ring;
//     ring.resize (100, 1, 1);
//     ring.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
//     ring.set_coup_sigma_x (0.0);
//     ring.set_coup_sigma_P_x (1);
//     ring.set_coup_gamma_x (0.0);
//     ring.equilibriumIC ();
//     double dt = 0.0001;
//     ring.set_poisson_pars (0.5, 1.0);
//     long steps = 200000;
//     long trans = 0;
//     long save_every = 2000;
//     std::ofstream fout("spacetime.dat");
//     std::ofstream sout("spikes.dat");
//     for (long step_i = 0; step_i < steps; step_i++) {
//         ring.step(dt, step_i * dt);
//         if ((step_i > trans) && (step_i % save_every == 0))
//             for (int n1 = 0; n1 < ring.N1; n1++) fout << step_i * dt << " " << n1 << " " << ring.state[ring.d3to1 (n1, 0, 0)].x << "\n";
//         ring.collect_max_x(sout);
//     }
//     sout.close();
//     fout.close();
//     std::stringstream plot_command;
//     plot_command << "gnuplot << 'EOF'\n";
//     plot_command << "set terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; set key tmargin center horizontal; set xlabel 'i'; set ylabel 't'; set cblabel 'x'; set xrange [0.5:100.5]; set yrange [:]; set cbrange [-2.5:2.5]; set xtics 20 out; set ytics 10 out; set cbtics 1; set mxtics 5; set mytics 5; set mcbtics 5; set output 'spacetime.png'; plot 'spacetime.dat' u (1.0+$2):1:3 w image notitle, 'spikes.dat' u (1.0+$2):1 w p pt 7 ps 0.9 lc rgb 'white' notitle;";
//     plot_command << "\nEOF";
//     std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
//     std::system (plot_command.str().c_str());
// }

// void SingleFHNPoissonTrace (long time, double freq) {
//     FHNens neu;
//     neu.resize(1, 1, 1);
//     neu.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
//     neu.equilibriumIC ();
//     neu.set_poisson_pars (freq, 1.0);
//     double dt = 0.0001;
//     long steps = time;
//     long save_every = 10;
//     std::ofstream fout("trace.dat");
//     std::ofstream fpulse("pulse.dat");
//     long prev_input = 0;
//     for (long step_i = 0; step_i < steps; step_i++) {
//         double t = step_i * dt;
//         neu.step(dt, t);
//         if (step_i % save_every == 0)
//             fout << t << " " << neu.state[0].x << "\n";
//         if (neu.n_input[0] > prev_input) {
//             prev_input = neu.n_input[0];
//             fpulse << t << " " << neu.state[0].x << "\n";
//         }
//     }
//     fout.close();
//     fpulse.close();
//     std::stringstream plot_command;
//     plot_command << "gnuplot << 'EOF'\n";
//     plot_command << "set terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; set key tmargin center horizontal; set xlabel 't'; set ylabel 'x'; set ytics 1; set output 'SingleFHNPoissonTrace.png'; plot 'trace.dat' u 1:2 w l lw 3 title 'x(t)', 'pulse.dat' u 1:(1) w p pt 7 ps 1.2 title 'Poisson pulses';";
//     plot_command << "\nEOF";
//     std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
//     std::system (plot_command.str().c_str());
// }

// void Ring100FoutOnFin () {
//     long nofth = 101;
//     long N1 = 100;
//     double dt = 0.0001;
//     long steps = 10000*5000;
//     std::vector<double> nin;
//     std::vector<double> nout;
//     std::vector<double> freq;
//     nin.resize(nofth);
//     nout.resize(nofth);
//     freq.resize(nofth);
//     #pragma omp parallel for schedule(dynamic)
//     for (long iter = 0; iter < nofth; iter++) {
//         FHNens ring;
//         freq[iter] = 0.0 + iter * 0.02;//pow(10.0, iter * 0.025 - 3.0);
//         ring.resize (N1, 1, 1);
//         ring.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
//         ring.set_coup_sigma_x (0.0);
//         ring.set_coup_sigma_P_x (1);
//         ring.set_coup_gamma_x (0.0);
//         ring.equilibriumIC ();
//         ring.set_poisson_pars (freq[iter], 1.0);
//         for (long step_i = 0; step_i < steps; step_i++) {
//             ring.step(dt, step_i * dt);
//         }
//         nin[iter] = 0.0;
//         nout[iter] = 0.0;
//         for (int n1 = 0; n1 < N1; n1++) {
//             nin[iter] += ring.n_input[ring.d3to1 (n1, 0, 0)];
//             nout[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 0)];
//         }
//     }
//     std::ofstream fout("FoutOnFin.dat");
//     for (long iter = 0; iter < nofth; iter++) {
//         fout << freq[iter] << " " << nin[iter]/(steps*dt*N1) << " " << nout[iter]/(steps*dt*N1) << "\n";
//     }
//     fout.close();
//     std::stringstream plot_command;
//     plot_command << "gnuplot << 'EOF'\n";
//     plot_command << "set terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; set key tmargin center horizontal; set xlabel 'f(Poisson)'; set ylabel 'f(sp)'; set xrange [1.0e-3:2.0]; set yrange [0:0.4]; set output 'FoutOnFin.png'; plot 'FoutOnFin.dat' u 1:3 w l lw 3 notitle, x w l lw 3 dt 2 notitle;";
//     plot_command << "\nEOF";
//     std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
//     std::system (plot_command.str().c_str());
// }

// void Ring100FoutOnFinMultiplex (double gamma) {
//     long nofth = 101;
//     long N1 = 100;
//     double dt = 0.0001;
//     long steps = 10000 * 5000;
//     std::vector<double> nout0;
//     std::vector<double> nout1;
//     std::vector<double> freq;
//     nout0.resize(nofth);
//     nout1.resize(nofth);
//     freq.resize(nofth);
//     #pragma omp parallel for schedule(dynamic)
//     for (long iter = 0; iter < nofth; iter++) {
//         FHNens ring;
//         freq[iter] = 0.0 + iter * 0.02;//pow(10.0, iter * 0.025 - 3.0);
//         ring.resize (N1, 1, 2);
//         ring.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
//         ring.set_coup_sigma_x (0.00, 0);
//         ring.set_coup_sigma_x (0.00, 1);
//         ring.set_coup_sigma_P_x (1, 0);
//         ring.set_coup_sigma_P_x (1, 1);
//         ring.set_coup_gamma_x (0.0);
//         ring.set_coup_gamma_x (gamma, 0, 1);
//         ring.set_coup_gamma_x (gamma, 1, 0);
//         ring.equilibriumIC ();
//         ring.set_poisson_pars (freq[iter], 1.0);
//         for (long step_i = 0; step_i < steps; step_i++) {
//             ring.step(dt, step_i * dt);
//         }
//         nout0[iter] = 0.0;
//         nout1[iter] = 0.0;
//         for (int n1 = 0; n1 < N1; n1++) {
//             nout0[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 0)];
//             nout1[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 1)];
//         }
//     }
//     std::stringstream filename;
//     std::stringstream filename_dat;
//     filename << "FoutOnFinMultiplex-g" << gamma;
//     filename_dat << filename.str() << ".dat";
//     std::ofstream fout(filename_dat.str());
//     for (long iter = 0; iter < nofth; iter++) {
//         fout << freq[iter] << " " << nout0[iter]/(steps*dt*N1) << " " << nout1[iter]/(steps*dt*N1) << "\n";
//     }
//     fout.close();
//     std::stringstream plot_command;
//     plot_command << "gnuplot << 'EOF'\n";
//     plot_command << "set terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; set key tmargin center horizontal; set xlabel 'f(Poisson)'; set ylabel 'f(sp)'; set xrange [1.0e-3:2.0]; set yrange [0:0.4]; set output '" << filename.str() << ".png'; plot '" << filename.str() << ".dat' u 1:(0.5*($2+$3)) w l lw 3 title 'σ=" << gamma << "', '../1layer/FoutOnFin.dat' u 1:3 w l lw 3 lc rgb 'black' title 'σ=0', x w l lw 3 lc rgb 'gray' dt 2 title 'f_{out}=f_{in}';";
//     plot_command << "\nEOF";
//     std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
//     std::system (plot_command.str().c_str());
// }

// void Ring10FoutOnFinMultiplex3 (double gamma) {
//     long nofth = 101;
//     long N1 = 10;
//     double dt = 0.0001;
//     long steps = 10000 * 5000;
//     std::vector<double> nout0;
//     std::vector<double> nout1;
//     std::vector<double> nout2;
//     std::vector<double> freq;
//     nout0.resize(nofth);
//     nout1.resize(nofth);
//     nout2.resize(nofth);
//     freq.resize(nofth);
//     #pragma omp parallel for schedule(dynamic)
//     for (long iter = 0; iter < nofth; iter++) {
//         FHNens ring;
//         freq[iter] = 0.0 + iter * 0.02;//pow(10.0, iter * 0.025 - 3.0);
//         ring.resize (N1, 1, 3);
//         ring.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
//         ring.set_coup_sigma_x (0.00, 0);
//         ring.set_coup_sigma_x (0.00, 1);
//         ring.set_coup_sigma_x (0.00, 2);
//         ring.set_coup_sigma_P_x (1, 0);
//         ring.set_coup_sigma_P_x (1, 1);
//         ring.set_coup_sigma_P_x (1, 2);
//         ring.set_coup_gamma_x (0.0);
//         ring.set_coup_gamma_x (gamma, 0, 1);
//         ring.set_coup_gamma_x (gamma, 1, 0);
//         ring.set_coup_gamma_x (gamma, 1, 2);
//         ring.set_coup_gamma_x (gamma, 2, 1);
//         ring.set_coup_gamma_x (gamma, 0, 2);
//         ring.set_coup_gamma_x (gamma, 2, 0);
//         ring.equilibriumIC ();
//         ring.set_poisson_pars (freq[iter], 1.0);
//         for (long step_i = 0; step_i < steps; step_i++) {
//             ring.step(dt, step_i * dt);
//         }
//         nout0[iter] = 0.0;
//         nout1[iter] = 0.0;
//         nout2[iter] = 0.0;
//         for (int n1 = 0; n1 < N1; n1++) {
//             nout0[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 0)];
//             nout1[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 1)];
//             nout2[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 2)];
//         }
//     }
//     std::stringstream filename;
//     std::stringstream filename_dat;
//     filename << "FoutOnFinMultiplex-g" << gamma;
//     filename_dat << filename.str() << ".dat";
//     std::ofstream fout(filename_dat.str());
//     for (long iter = 0; iter < nofth; iter++) {
//         fout << freq[iter] << " " << nout0[iter]/(steps*dt*N1) << " " << nout1[iter]/(steps*dt*N1) << " " << nout2[iter]/(steps*dt*N1) << "\n";
//     }
//     fout.close();
//     std::stringstream plot_command;
//     plot_command << "gnuplot << 'EOF'\n";
//     plot_command << "set terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; set key tmargin center horizontal; set xlabel 'f(Poisson)'; set ylabel 'f(sp)'; set xrange [1.0e-3:2.0]; set yrange [0:0.4]; set output '" << filename.str() << ".png'; plot '" << filename.str() << ".dat' u 1:(($2+$3+$4)/3.0) w l lw 3 title 'σ=" << gamma << "', '../1layer/FoutOnFin.dat' u 1:3 w l lw 3 lc rgb 'black' title 'σ=0', x w l lw 3 lc rgb 'gray' dt 2 title 'f_{out}=f_{in}';";
//     plot_command << "\nEOF";
//     std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
//     std::system (plot_command.str().c_str());
// }

void RingFoutOnFinSizesRanges (double sigmaloc, double freqloc) {
    long nofth = 0;
    for (long totsize = 5; totsize <= 125; totsize+=30) {
        for (long range = 1; range < (totsize/2); range++) {
            nofth++;
        }
    }
    double dt = 0.0001;
    long steps = 10000 * 20000;
    std::vector<double> nout0;
    std::vector<long> sizes;
    std::vector<long> ranges;
    nout0.resize(nofth);
    sizes.resize(nofth);
    ranges.resize(nofth);
    long nofthloc = 0;
    for (long totsize = 5; totsize <= 125; totsize+=30) {
        for (long range = 1; range < (totsize/2); range++) {
            sizes[nofthloc] = totsize;
            ranges[nofthloc] = range;
            nofthloc++;
        }
    }
    #pragma omp parallel for schedule(dynamic)
    for (long iter = 0; iter < nofth; iter++) {
        FHNens ring;
        ring.resize (sizes[iter], 1, 1);
        ring.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
        ring.set_coup_sigma_x (sigmaloc, 0);
        ring.set_coup_sigma_P_x (ranges[iter], 0);
        ring.set_coup_gamma_x (0.0);
        ring.equilibriumIC ();
        ring.set_poisson_pars (0.0, 0.0);
        ring.set_poisson_pars (freqloc, 1.0, 0);
        for (long step_i = 0; step_i < steps; step_i++) {
            ring.step(dt, step_i * dt);
        }
        nout0[iter] = 0.0;
        for (int n1 = 0; n1 < sizes[iter]; n1++) {
            nout0[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 0)];
        }
    }
    std::stringstream filename;
    std::stringstream filename_dat;
    filename << "FoutOnFinMultiplex-s" << sigmaloc << "-f" << freqloc;
    filename_dat << filename.str() << ".dat";
    std::ofstream fout(filename_dat.str());
    long currentN = sizes[0];
    for (long iter = 0; iter < nofth; iter++) {
        if (sizes[iter]!=currentN) {
            fout << "\n";
            currentN = sizes[iter];
        }
        fout << sizes[iter] << " " << ranges[iter] << " " << nout0[iter]/(steps*dt*sizes[iter]) << "\n";
    }
    fout.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << 'EOF'\n";
    plot_command << "set terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; set key tmargin center horizontal; set xlabel 'r'; set ylabel 'f(sp)'; set xrange [0.0:0.5]; set yrange [0:0.4]; set output '" << filename.str() << ".png'; plot '" << filename.str() << ".dat' u ($2/$1):($1==5?$3:1/0) w l lw 3 title 'N=5', '' u ($2/$1):($1==35?$3:1/0) w l lw 3 title 'N=35', '' u ($2/$1):($1==65?$3:1/0) w l lw 3 title 'N=65', '' u ($2/$1):($1==95?$3:1/0) w l lw 3 title 'N=95', '' u ($2/$1):($1==125?$3:1/0) w l lw 3 title 'N=125';";
    plot_command << "\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

void RingsFoutOnFinMultiplex3 (long Nloc, long Ploc, double sigmaloc, double gammaloc) {
    long nofth = 101;
    double dt = 0.0001;
    long steps = 10000 * 50000;
    std::vector<double> nout0;
    std::vector<double> nout1;
    std::vector<double> nout2;
    std::vector<double> freq;
    nout0.resize(nofth);
    nout1.resize(nofth);
    nout2.resize(nofth);
    freq.resize(nofth);
    #pragma omp parallel for schedule(dynamic)
    for (long iter = 0; iter < nofth; iter++) {
        FHNens ring;
        freq[iter] = 0.0 + iter * 0.02;//pow(10.0, iter * 0.025 - 3.0);
        ring.resize (Nloc, 1, 3);
        ring.setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
        ring.set_coup_sigma_x (sigmaloc, 0);
        ring.set_coup_sigma_x (sigmaloc, 1);
        ring.set_coup_sigma_x (sigmaloc, 2);
        ring.set_coup_sigma_P_x (Ploc, 0);
        ring.set_coup_sigma_P_x (Ploc, 1);
        ring.set_coup_sigma_P_x (Ploc, 2);
        ring.set_coup_gamma_x (0.0);
        ring.set_coup_gamma_x (gammaloc, 2, 1);
        ring.set_coup_gamma_x (gammaloc, 1, 0);
        ring.equilibriumIC ();
        ring.set_poisson_pars (0.0, 0.0);
        ring.set_poisson_pars (freq[iter], 1.0, 0);
        for (long step_i = 0; step_i < steps; step_i++) {
            ring.step(dt, step_i * dt);
        }
        nout0[iter] = 0.0;
        nout1[iter] = 0.0;
        nout2[iter] = 0.0;
        for (int n1 = 0; n1 < Nloc; n1++) {
            nout0[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 0)];
            nout1[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 1)];
            nout2[iter] += ring.n_spikes[ring.d3to1 (n1, 0, 2)];
        }
    }
    std::stringstream filename;
    std::stringstream filename_dat;
    filename << "FoutOnFinMultiplex-g" << sigmaloc;
    filename_dat << filename.str() << ".dat";
    std::ofstream fout(filename_dat.str());
    for (long iter = 0; iter < nofth; iter++) {
        fout << freq[iter] << " " << nout0[iter]/(steps*dt*Nloc) << " " << nout1[iter]/(steps*dt*Nloc) << " " << nout2[iter]/(steps*dt*Nloc) << "\n";
    }
    fout.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << 'EOF'\n";
    plot_command << "set terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; set key tmargin center horizontal; set xlabel 'f(Poisson)'; set ylabel 'f(sp)'; set xrange [1.0e-3:2.0]; set yrange [0:0.4]; set output '" << filename.str() << ".png'; plot '" << filename.str() << ".dat' u 1:2 w l lw 3 title 'σ=" << sigmaloc << "', '../1layer/FoutOnFin.dat' u 1:3 w l lw 3 lc rgb 'black' title 'σ=0', x w l lw 3 lc rgb 'gray' dt 2 title 'f_{out}=f_{in}';";
    plot_command << "\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
    std::stringstream plot_command2;
    plot_command2 << "gnuplot << 'EOF'\n";
    plot_command2 << "set terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; set key tmargin center horizontal; set xlabel 'f(Poisson)'; set ylabel 'f(sp)'; set xrange [1.0e-3:2.0]; set yrange [0:0.4]; set output '" << filename.str() << "full.png'; plot '" << filename.str() << ".dat' u 1:2 w l lw 3 title 'σ=" << sigmaloc << "', '' u 1:3 w l lw 3 title '2nd layer', '' u 1:4 w l lw 3 title '3rd layer', '../1layer/FoutOnFin.dat' u 1:3 w l lw 3 lc rgb 'black' title 'σ=0', x w l lw 3 lc rgb 'gray' dt 2 title 'f_{out}=f_{in}';";
    plot_command2 << "\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command2.str() << "\n\n";
    std::system (plot_command2.str().c_str());
}

void Ring100 () {
    // SingleFHNPoissonTrace(200000, 2.00);
    // Ring100std ();
    // Ring100FoutOnFin ();
    // for (double g = 0.0; g < 0.2; g += 0.01)
    //     Ring100FoutOnFinMultiplex (g);
    // for (double g = 0.0; g < 0.2; g += 0.01)
    //     Ring10FoutOnFinMultiplex3 (g);
    // for (double g = 0.25; g < 1.0; g += 0.05)
    //     Ring10FoutOnFinMultiplex3 (g);
    // for (double s = 0.0; s < 0.2; s += 0.01)
    //     RingsFoutOnFinMultiplex3 (5, 1, s, 0.1);
    // for (double s = 0.2; s < 1.0; s += 0.05)
    //     RingsFoutOnFinMultiplex3 (5, 1, s, 0.1);
    for (double s = 0.50; s < 0.56; s += 0.05)
        for (double f = 0.30; f < 0.31; f += 0.05)
            RingFoutOnFinSizesRanges (s, f);
}

#endif
