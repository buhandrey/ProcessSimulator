#ifndef FHN_ENS_H
#define FHN_ENS_H

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

struct FHNInput {
    double ix = 0.0;
    double iy = 0.0;
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

inline FHNDeriv rhs (const FHNState& s, const FHNParams& p, const FHNInput& in) {
    return {
        (s.x * (1.0 - s.x * s.x / 3.0) - s.y + p.R * p.I + in.ix) / p.epsilon,
        s.x + p.alpha - p.beta * s.y + in.iy
    };
}

class FHNens {
public:
    int N1 = 0, N2 = 0, M = 0;
    std::vector<std::vector<std::vector<FHNState>>> state;
    std::vector<std::vector<std::vector<FHNParams>>> pars;
    std::vector<std::vector<std::vector<FHNState>>> tmp;
    std::vector<std::vector<std::vector<double>>> tx_event;
    std::vector<std::vector<std::vector<double>>> dW;
    std::vector<int> coup_sigma_P_x;
    std::vector<double> coup_sigma_x;
    std::vector<std::vector<double>> coup_gamma_x;
    FHNens () {
        N1 = 0;
        N2 = 0;
        M = 0;
    }
    void resize (int newN1, int newN2, int newM) {
        N1 = newN1; N2 = newN2; M = newM;
        state.resize (N1);
        pars.resize (N1);
        tmp.resize (N1);
        tx_event.resize(N1);
        dW.resize(N1);
        for (int n1 = 0; n1 < N1; n1++) {
            state[n1].resize (N2);
            pars[n1].resize (N2);
            tmp[n1].resize (N2);
            tx_event[n1].resize(N2);
            dW[n1].resize(N2);
            for (int n2 = 0; n2 < N2; n2++) {
                state[n1][n2].resize (M);
                pars[n1][n2].resize (M);
                tmp[n1][n2].resize (M);
                tx_event[n1][n2].resize (M);
                dW[n1][n2].resize (M);
                for (int m = 0; m < M; m++) {
                    tx_event[n1][n2][m] = 0.0;
                }
            }
        }
        coup_sigma_x.resize (M);
        coup_sigma_P_x.resize (M);
        coup_gamma_x.resize (M, std::vector<double> (M, 0.0));
    }
    void randomIC (double rx, double ry) {
        std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<double> dist(0.0, 2.0 * M_PI);
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    double theta = dist(gen);
                    state[n1][n2][m].x = rx * std::cos(theta);
                    state[n1][n2][m].y = ry * std::sin(theta);
                }
            }
        }
    }
    void equilibriumIC () {
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    state[n1][n2][m].x = -pars[n1][n2][m].alpha;
                    state[n1][n2][m].y = pars[n1][n2][m].alpha * (pars[n1][n2][m].alpha*pars[n1][n2][m].alpha/3.0 - 1.0);
                }
            }
        }
    }
    void setFHNpars (double new_alpha, double new_beta, double new_epsilon, double new_I, double new_R) {
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    pars[n1][n2][m].alpha = new_alpha;
                    pars[n1][n2][m].beta = new_beta;
                    pars[n1][n2][m].epsilon = new_epsilon;
                    pars[n1][n2][m].I = new_I;
                    pars[n1][n2][m].R = new_R;
                }
            }
        }
    }
    void set_poisson_pars (double new_f, double new_a) {
        for (int n1 = 0; n1 < N1; n1++) {
            for (int n2 = 0; n2 < N2; n2++) {
                for (int m = 0; m < M; m++) {
                    pars[n1][n2][m].poisson_x_f = new_f;
                    pars[n1][n2][m].poisson_x_a = new_a;
                }
            }
        }
        init_events();
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
        coup_gamma_x[m1][m2] = new_gamma;
    }
    double coupling_x (int n1, int m, const std::vector<std::vector<std::vector<FHNState>>>& S) {
        double cs = 0.0;
        if (N2==1) {
            for (int n1s = n1 - coup_sigma_P_x[m]; n1s <= n1 + coup_sigma_P_x[m]; n1s++) {
                int j = (n1s + N1) % N1;
                if (j != n1) cs += coup_sigma_x[m] * (S[j][0][m].x - S[n1][0][m].x);
            }
        }
        double cg = 0.0;
        for (int m2 = 0; m2 < M; m2++)
            cg += coup_gamma_x[m][m2] * (S[n1][0][m2].x - S[n1][0][m].x);
        return cs + cg;
    }
    double next_event_time(double freq, std::mt19937& gen) {
        std::exponential_distribution<double> exp(freq);
        return exp(gen);
    }
    void init_events() {
        static thread_local std::mt19937 gen(std::random_device{}());
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            tx_event[n1][n2][m] = next_event_time(pars[n1][n2][m].poisson_x_f, gen);
        }
    }
    void step(double dt) {
        static thread_local std::mt19937 gen(std::random_device{}());
        static thread_local std::normal_distribution<double> normal(0.0, 1.0);
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
             auto &p = pars[n1][n2][m];
             tx_event[n1][n2][m] -= dt;
             if (tx_event[n1][n2][m] <= 0.0 && p.poisson_x_f > 0.0) {
                state[n1][n2][m].x += p.poisson_x_a;
                std::exponential_distribution<double> exp(p.poisson_x_f);
                tx_event[n1][n2][m] = exp(gen);
             }
        }
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            auto &s = state[n1][n2][m];
            auto &p = pars[n1][n2][m];
            double cx = coupling_x(n1, m, state);
            dW[n1][n2][m] = std::sqrt(dt) * normal(gen);
            FHNInput in {cx, 0.0};
            auto f0 = rhs(s, p, in);
            tmp[n1][n2][m].x = s.x + dt * f0.dx;
            tmp[n1][n2][m].y = s.y + dt * f0.dy + p.gauss_y_s * dW[n1][n2][m];
        }
        for (int n1 = 0; n1 < N1; n1++) for (int n2 = 0; n2 < N2; n2++) for (int m = 0; m < M; m++) {
            auto &s = state[n1][n2][m];
            auto &p = pars[n1][n2][m];
            double cx0 = coupling_x(n1, m, state);
            double cx1 = coupling_x(n1, m, tmp);
            FHNInput in0 {cx0, 0.0};
            FHNInput in1 {cx1, 0.0};
            auto f0 = rhs(s, p, in0);
            auto f1 = rhs(tmp[n1][n2][m], p, in1);
            s.x += 0.5 * dt * (f0.dx + f1.dx);
            s.y += 0.5 * dt * (f0.dy + f1.dy) + p.gauss_y_s * dW[n1][n2][m];
        }
    }
    void Ring100 () {
        resize (100, 1, 1);
        set_coup_sigma_x (0.0);
        set_coup_gamma_x (0.0);
        set_coup_sigma_P_x (1);
        setFHNpars (1.2, 0.0, 0.01, 0.0, 0.0);
        equilibriumIC ();
        set_poisson_pars (0.1, 1.0);
        double dt = 0.0001;
        long steps = 50000;
        long trans = 0;
        long save_every = 500;
        std::ofstream fout("spacetime.dat");
        for (long step_i = 0; step_i < steps; step_i++) {
            step(dt);
            double t = step_i * dt;
            if ((step_i > trans) && (step_i % save_every == 0))
                for (int n1 = 0; n1 < N1; n1++) fout << t << " " << n1 << " " << state[n1][0][0].x << "\n";
        }
        fout.close();
        std::stringstream plot_command;
        plot_command << "gnuplot << 'EOF'\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20';";
        plot_command << " unset warnings; set key tmargin center horizontal;";
        plot_command << " set xlabel 'i'; set ylabel 't'; set cblabel 'x';";
        plot_command << " set xrange [1:100]; set yrange [:]; set cbrange [-2.5:2.5];";
        plot_command << " set xtics 20 out; set ytics 1 out; set cbtics 1;";
        plot_command << " set mxtics 5; set mytics 5; set mcbtics 5;";
        plot_command << " set output 'spacetime.png';";
        plot_command << " plot 'spacetime.dat' using (1.0+$2):1:3 w image notitle;\nEOF";
        std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
        std::system (plot_command.str().c_str());
    }
};

#endif