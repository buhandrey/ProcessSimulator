#ifndef RESEARCHER_H
#define RESEARCHER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

template <typename specificmodel>
class researcher {
public:
    researcher () = default;
    ~researcher() = default;
    /*Методы для класса FHNsimple*/
    void TS0xy_x0_y0_alpha_epsilon_dt_t0_t1 (double x0, double y0, double alpha, double epsilon, double dt, long t0, long t1);
    void PP0xy_x0_y0_alpha_epsilon_dt_t0_t1 (double x0, double y0, double alpha, double epsilon, double dt, long t0, long t1);
    /*Методы для класса Rulkov_map*/
    void TS0xy_x0_y0_alpha_mu_sigma_t0_t1 (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1);
    void TS0xyLCE_x0_y0_alpha_mu_sigma_t0_t1 (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1);
    void PP0xy_x0_y0_alpha_mu_sigma_t0_t1 (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1);
    void PP0xyLCE_x0_y0_alpha_mu_sigma_t0_t1 (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1);
    void PDF0x_x0_y0_alpha_mu_sigma_t0_t1_winsize (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1, long PDF_win_size);
    void X0onS_s0_ds_s1_x0_y0_alpha_mu_t0_t1 (double s0, double ds, double s1, double x0, double y0, double alpha, double mu, long t0, long t1);
    void X0onA_a0_da_a1_x0_y0_mu_sigma_t0_t1 (double a0, double da, double a1, double x0, double y0, double mu, double sigma, long t0, long t1);
    void L12onS_s0_ds_s1_x0_y0_alpha_mu_t0_t1 (double s0, double ds, double s1, double x0, double y0, double alpha, double mu, long t0, long t1);
    void L12onA_a0_da_a1_x0_y0_mu_sigma_t0_t1 (double a0, double da, double a1, double x0, double y0, double mu, double sigma, long t0, long t1);
    void ExonS_s0_ds_s1_x0_y0_alpha_mu_t0_t1_winsize (double s0, double ds, double s1, double x0, double y0, double alpha, double mu, long t0, long t1, long ws);
    void ExonA_a0_da_a1_x0_y0_mu_sigma_t0_t1_winsize (double a0, double da, double a1, double x0, double y0, double mu, double sigma, long t0, long t1, long ws);
    void Min0x_x0_y0_alpha_mu_sigma_t0_t1(double x0, double y0, double alpha, double mu, double sigma, long t0, long t1);
    /*Методы для класса levy*/
    void LevySourceStat_alpha_beta_sigma_mu_idum_N (double alpha, double beta, double sigma, double mu, long idum, long N);
    /*Методы для класса poisson*/
    void PoissonSourceStat_alpha_beta_sigma_mu_idum_N (double freq, long idum, long N);
};

template <typename specificmodel>
void researcher<specificmodel>::TS0xy_x0_y0_alpha_epsilon_dt_t0_t1 (double x0, double y0, double alpha, double epsilon, double dt, long t0, long t1) {
    specificmodel model;
    model.set_par_a(alpha);
    model.set_par_e(epsilon);
    model.set_buffer_size(2);
    model.set_x0(x0);
    model.set_y0(y0);
    model.set_dt(dt);
    std::stringstream datfile, pngfile;
    datfile << "TS0xy" << model.get_basename() << ".dat";
    pngfile << "TS0xy" << model.get_basename() << ".png";
    std::ofstream datafile;
    datafile.open(datfile.str());
    long every = (t1-t0)/100000;
    if (every==0) every = 1;
    for (long t = 0; t <= t1; t++) {
        if ( (t>=t0) && ((t%every) == 0) ) {
            datafile << t*model.getdt();
            for (long n = 0; n < model.varnumber(); n++) {
                if ( (n == 0) || (n == 1) )
                    datafile << " " << model.get_vars(n);
            }
            datafile << "\n";
        }
        model.step();
    }
    datafile.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 't'; set ylabel 'x, y'; set output '" << pngfile.str() << "'; stat '" << datfile.str() << "' u 2:3 nooutput; maxval = (STATS_max_x>STATS_max_y)?STATS_max_x:STATS_max_y; minval = (STATS_min_x<STATS_min_y)?STATS_min_x:STATS_min_y; range = maxval-minval; set yrange [minval - 0.05 * range : maxval + 0.05 * range]; ";
    plot_command << "plot '" << datfile.str() << "' u 1:2 w l lw 3 lc rgb 0x008040 title 'x(t)', '' u 1:3 w l lw 3 lc rgb 0x8b0000 title 'y(t)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

template <typename specificmodel>
void researcher<specificmodel>::PP0xy_x0_y0_alpha_epsilon_dt_t0_t1 (double x0, double y0, double alpha, double epsilon, double dt, long t0, long t1) {
    specificmodel model;
    model.set_par_a(alpha);
    model.set_par_e(epsilon);
    model.set_buffer_size(2);
    model.set_x0(x0);
    model.set_y0(y0);
    model.set_dt(dt);
    std::stringstream datfile, pngfile;
    datfile << "PP0xy" << model.get_basename() << ".dat";
    pngfile << "PP0xy" << model.get_basename() << ".png";
    std::ofstream datafile;
    datafile.open(datfile.str());
    long every = (t1-t0)/100000;
    if (every==0) every = 1;
    for (long t = 0; t <= t1; t++) {
        if ( (t>=t0) && ((t%every) == 0) ) {
            datafile << t*model.getdt();
            for (long n = 0; n < model.varnumber(); n++) {
                if ( (n == 0) || (n == 1) )
                    datafile << " " << model.get_vars(n);
            }
            datafile << "\n";
        }
        model.step();
    }
    datafile.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 'x'; set ylabel 'y'; set output '" << pngfile.str() << "'; stat '" << datfile.str() << "' u 2:3 nooutput; xrange = STATS_max_x-STATS_min_x; yrange = STATS_max_y-STATS_min_y; set xrange [STATS_min_x-0.05*xrange:STATS_max_x+0.05*xrange]; set yrange [STATS_min_y-0.05*yrange:STATS_max_y+0.05*yrange]; ";
    plot_command << "plot '" << datfile.str() << "' u 2:3 w l lw 3 lc rgb 0x50191970 title 'y(x)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

template <typename specificmodel>
void researcher<specificmodel>::TS0xy_x0_y0_alpha_mu_sigma_t0_t1 (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1) {
    specificmodel model;
    model.set_par_a(alpha);
    model.set_par_m(mu);
    model.set_par_s(sigma);
    model.set_buffer_size(2);
    model.set_x0(x0);
    model.set_y0(y0);
    std::stringstream datfile, pngfile;
    datfile << "TS0xy" << model.get_basename() << ".dat";
    pngfile << "TS0xy" << model.get_basename() << ".png";
    std::ofstream datafile;
    datafile.open(datfile.str());
    long every = (t1-t0)/100000;
    if (every==0) every = 1;
    for (long t = 0; t <= t1; t++) {
        if ( (t>=t0) && ((t%every) == 0) ) {
            datafile << t*model.getdt();
            for (long n = 0; n < model.varnumber(); n++) {
                if ( (n == 0) || (n == 1) )
                    datafile << " " << model.get_vars(n);
            }
            datafile << "\n";
        }
        model.step();
    }
    datafile.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 't'; set ylabel 'x, y'; set output '" << pngfile.str() << "'; stat '" << datfile.str() << "' u 2:3 nooutput; maxval = (STATS_max_x>STATS_max_y)?STATS_max_x:STATS_max_y; minval = (STATS_min_x<STATS_min_y)?STATS_min_x:STATS_min_y; range = maxval-minval; set yrange [minval - 0.05 * range : maxval + 0.05 * range]; ";
    plot_command << "plot '" << datfile.str() << "' u 1:2 w l lw 3 lc rgb 0x008040 title 'x(t)', '' u 1:3 w l lw 3 lc rgb 0x8b0000 title 'y(t)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

template <typename specificmodel>
void researcher<specificmodel>::TS0xyLCE_x0_y0_alpha_mu_sigma_t0_t1 (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1) {
    specificmodel model;
    model.lyap_init();
    model.set_par_a(alpha);
    model.set_par_m(mu);
    model.set_par_s(sigma);
    model.set_buffer_size(2);
    model.set_x0(x0);
    model.set_y0(y0);
    std::stringstream datfile, pngfile;
    datfile << "TS0xyLCE" << model.get_basename() << ".dat";
    pngfile << "TS0xyLCE" << model.get_basename() << ".png";
    std::ofstream datafile;
    datafile.open(datfile.str());
    model.normalization_init();
    long every = (t1-t0)/100000;
    if (every==0) every = 1;
    for (long t = 0; t <= t1; t++) {
        if ( (t>=t0) && ((t%every) == 0) ) {
            datafile << t*model.getdt();
            for (long n = 0; n < model.varnumber(); n++) {
                if ( (n == 0) || (n == 1) )
                    datafile << " " << model.get_vars(n);
            }
            datafile << "\n";
        }
        model.step();
        if (t==t0) model.normalization_reset();
        else model.GramSmidt();
    }
    model.GramSmidt();
    datafile.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 't'; set ylabel 'x, y'; set output '" << pngfile.str() << "'; stat '" << datfile.str() << "' u 2:3 nooutput; maxval = (STATS_max_x>STATS_max_y)?STATS_max_x:STATS_max_y; minval = (STATS_min_x<STATS_min_y)?STATS_min_x:STATS_min_y; range = maxval-minval; set yrange [minval - 0.05 * range : maxval + 0.05 * range]; ";
    plot_command << "set label 'L0 = " << model.l(0) << "' at screen 0.2,0.95;  set label 'L1 = " << model.l(1) << "' at screen 0.7,0.95; ";
    plot_command << "plot '" << datfile.str() << "' u 1:2 w l lw 3 lc rgb 0x008040 title 'x(t)', '' u 1:3 w l lw 3 lc rgb 0x8b0000 title 'y(t)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

template <typename specificmodel>
void researcher<specificmodel>::PP0xy_x0_y0_alpha_mu_sigma_t0_t1 (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1) {
    specificmodel model;
    model.set_par_a(alpha);
    model.set_par_m(mu);
    model.set_par_s(sigma);
    model.set_buffer_size(2);
    model.set_x0(x0);
    model.set_y0(y0);
    long every = (t1-t0)/100000;
    if (every==0) every = 1;
    std::stringstream datfile, pngfile;
    datfile << "PP0xy" << model.get_basename() << ".dat";
    pngfile << "PP0xy" << model.get_basename() << ".png";
    std::ofstream datafile(datfile.str());
    for (long t = 0; t <= t1; t++) {
        if ( (t>=t0) && (t<=t1) && ((t%every) == 0) ) {
            datafile << t*model.getdt();
            for (long n = 0; n < model.varnumber(); n++) {
                if ( (n == 0) || (n == 1) )
                    datafile << " " << model.get_vars(n);
            }
            datafile << "\n";
        }
        model.step();
    }
    datafile.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 'x'; set ylabel 'y'; set output '" << pngfile.str() << "'; stat '" << datfile.str() << "' u 2:3 nooutput; xrange = STATS_max_x-STATS_min_x; yrange = STATS_max_y-STATS_min_y; set xrange [STATS_min_x-0.05*xrange:STATS_max_x+0.05*xrange]; set yrange [STATS_min_y-0.05*yrange:STATS_max_y+0.05*yrange]; ";
    plot_command << "plot '" << datfile.str() << "' u 2:3:(0.8) w p pt 7 ps variable lc rgb 0x50191970 title 'y(x)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

template <typename specificmodel>
void researcher<specificmodel>::PP0xyLCE_x0_y0_alpha_mu_sigma_t0_t1 (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1) {
    specificmodel model;
    model.lyap_init();
    model.set_par_a(alpha);
    model.set_par_m(mu);
    model.set_par_s(sigma);
    model.set_buffer_size(2);
    model.set_x0(x0);
    model.set_y0(y0);
    long every = (t1-t0)/100000;
    if (every==0) every = 1;
    std::stringstream datfile, pngfile;
    datfile << "PP0xyLCE" << model.get_basename() << ".dat";
    pngfile << "PP0xyLCE" << model.get_basename() << ".png";
    std::ofstream datafile(datfile.str());
    model.normalization_init();
    for (long t = 0; t <= t1; t++) {
        if ( (t>=t0) && (t<=t1) && ((t%every) == 0) ) {
            datafile << t*model.getdt();
            for (long n = 0; n < model.varnumber(); n++) {
                if ( (n == 0) || (n == 1) )
                    datafile << " " << model.get_vars(n);
            }
            datafile << "\n";
        }
        model.step();
        if (t==t0) model.normalization_reset();
        else model.GramSmidt();
    }
    model.GramSmidt();
    datafile.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 'x'; set ylabel 'y'; set output '" << pngfile.str() << "'; stat '" << datfile.str() << "' u 2:3 nooutput; xrange = STATS_max_x-STATS_min_x; yrange = STATS_max_y-STATS_min_y; set xrange [STATS_min_x-0.05*xrange:STATS_max_x+0.05*xrange]; set yrange [STATS_min_y-0.05*yrange:STATS_max_y+0.05*yrange]; ";
    plot_command << "set label 'L0 = " << model.l(0) << "' at screen 0.2,0.95;  set label 'L1 = " << model.l(1) << "' at screen 0.7,0.95; ";
    plot_command << "plot '" << datfile.str() << "' u 2:3:(0.8) w p pt 7 ps variable lc rgb 0x50191970 title 'y(x)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

template <typename specificmodel>
void researcher<specificmodel>::PDF0x_x0_y0_alpha_mu_sigma_t0_t1_winsize (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1, long PDF_win_size) {
    specificmodel model;
    model.set_par_a(alpha);
    model.set_par_m(mu);
    model.set_par_s(sigma);
    model.set_buffer_size(2);
    model.set_x0(x0);
    model.set_y0(y0);
    model.prob_init_vars(0, PDF_win_size);
    for (long t = 0; t <= t1; t++) {
        model.step();
        if (t>t0) model.update_stat_vars(0);
    }
    std::stringstream datfile, pngfile;
    datfile << "PDF0x" << model.get_basename() << ".dat";
    pngfile << "PDF0x" << model.get_basename() << ".png";
    std::ofstream datafile(datfile.str());
    for (long i = 0; i < model.get_bins_number(0); i++) {
        double tval = model.get_prob_dens_vars(0, i);
        if (tval>0.0) {
            double mid = 0.5 * (model.get_prob_right_vars(0, i) + model.get_prob_left_vars(0, i));
            double width = model.get_prob_right_vars(0, i) - model.get_prob_left_vars(0, i);
            datafile << mid << " " << tval << " " << width << "\n";
        }
    }
    datafile.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 'x'; set ylabel 'p(x)'; set output '" << pngfile.str() << "'; stat '" << datfile.str() << "' u 1:2 nooutput; maxval = STATS_max_y; set yrange[0:maxval*1.1]; minval = STATS_min_x; maxval = STATS_max_x; range = maxval-minval+0.01; set xrange[minval-0.05*range:maxval+0.05*range]; set label 'entropy = " << model.get_entropy2_vars(0) << " bits' at screen 0.7,0.95; plot '" << datfile.str() << "' u 1:2:3 w boxes lc rgb 0x00556b2f title 'p(x)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
}

template <typename specificmodel>
void researcher<specificmodel>::X0onS_s0_ds_s1_x0_y0_alpha_mu_t0_t1 (double s0, double ds, double s1, double x0, double y0, double alpha, double mu, long t0, long t1) {
    std::stringstream datfile_for, datfile_rev, pngfile_for, pngfile_rev, pngfile_both;
    double par_for = s0;
    double par_rev = s1;
    double x0_for = x0;
    double x0_rev = x0;
    double y0_for = y0;
    double y0_rev = y0;
    std::ofstream datafile_for, datafile_rev;
    for (long i = 0; par_for <= s1; i++) {
        specificmodel model_for, model_rev;
        model_for.set_par_a(alpha);
        model_rev.set_par_a(alpha);
        model_for.set_par_m(mu);
        model_rev.set_par_m(mu);
        model_for.set_par_s(par_for);
        model_rev.set_par_s(par_rev);
        model_for.set_buffer_size(2);
        model_rev.set_buffer_size(2);
        model_for.set_x0(x0_for);
        model_rev.set_x0(x0_rev);
        model_for.set_y0(y0_for);
        model_rev.set_y0(y0_rev);
        if (i==0) {
            datfile_for << "X0onS" << model_for.get_basename() << ".dat";
            datfile_rev << "X0onS" << model_for.get_basename() << "_rev.dat";
            pngfile_for << "X0onS" << model_for.get_basename() << ".png";
            pngfile_rev << "X0onS" << model_for.get_basename() << "_rev.png";
            pngfile_both << "X0onS" << model_for.get_basename() << "_both.png";
            datafile_for.open(datfile_for.str());
            datafile_rev.open(datfile_rev.str());
        }
        for (long t = 0; t <= t1; t++) {
            if ( (t>=t0) && (t<=t1) ) {
                datafile_for << par_for << " " << t*model_for.getdt();
                datafile_rev << par_rev << " " << t*model_rev.getdt();
                for (long n = 0; n < model_for.varnumber(); n++)
                    if ( (n == 0) ) {
                        datafile_for << " " << model_for.get_vars(n);
                        datafile_rev << " " << model_rev.get_vars(n);
                    }
                datafile_for << "\n";
                datafile_rev << "\n";
            }
            model_for.step();
            model_rev.step();
        }
        par_for+=ds;
        par_rev-=ds;
        x0_for = model_for.get_vars(0);
        x0_rev = model_rev.get_vars(0);
        y0_for = model_for.get_vars(1);
        y0_rev = model_rev.get_vars(1);
    }
    datafile_for.close();
    datafile_rev.close();
    std::string parname = "σ";
    for (long picn = 0; picn < 3; picn++) {
        std::stringstream plot_command;
        plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel '" << parname << "'; set ylabel 'x'; set output '";
        if (picn==0)
            plot_command << pngfile_both.str();
        if (picn==1)
            plot_command << pngfile_for.str();
        if (picn==2)
            plot_command << pngfile_rev.str();
        plot_command << "'; stat '" << datfile_for.str() << "' u 1:3 nooutput; miny = STATS_min_y; maxy = STATS_max_y; stat '" << datfile_rev.str() << "' u 1:3 nooutput; miny = STATS_min_y<miny?STATS_min_y:miny; maxy = STATS_max_y<maxy?STATS_max_y:maxy; range = maxy-miny; set yrange [miny-0.05*range:maxy+0.05*range]; ";
        if (picn==0)
            plot_command << "plot '" << datfile_for.str() << "' u 1:3:(0.2) w p pt 7 ps variable lc rgb 0x80ff8040 title 'x(increasing " << parname << ")', '" << datfile_rev.str() << "' u 1:3:(0.2) w p pt 7 ps variable lc rgb 0x802e8b57 title 'x(decreasing " << parname << ")';\nEOF";
        if (picn==1)
            plot_command << "plot '" << datfile_for.str() << "' u 1:3:(0.2) w p pt 7 ps variable lc rgb 0x80ff8040 title 'x(increasing " << parname << ")';\nEOF";
        if (picn==2)
            plot_command << "plot '" << datfile_rev.str() << "' u 1:3:(0.2) w p pt 7 ps variable lc rgb 0x802e8b57 title 'x(decreasing " << parname << ")';\nEOF";
        std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
        std::system (plot_command.str().c_str());
    }
}

template <typename specificmodel>
void researcher<specificmodel>::X0onA_a0_da_a1_x0_y0_mu_sigma_t0_t1 (double a0, double da, double a1, double x0, double y0, double mu, double sigma, long t0, long t1) {
    std::stringstream datfile_for, datfile_rev, pngfile_for, pngfile_rev, pngfile_both;
    double par_for = a0;
    double par_rev = a1;
    double x0_for = x0;
    double x0_rev = x0;
    double y0_for = y0;
    double y0_rev = y0;
    std::ofstream datafile_for, datafile_rev;
    for (long i = 0; par_for <= a1; i++) {
        specificmodel model_for, model_rev;
        model_for.set_par_a(par_for);
        model_rev.set_par_a(par_rev);
        model_for.set_par_m(mu);
        model_rev.set_par_m(mu);
        model_for.set_par_s(sigma);
        model_rev.set_par_s(sigma);
        model_for.set_buffer_size(2);
        model_rev.set_buffer_size(2);
        model_for.set_x0(x0_for);
        model_rev.set_x0(x0_rev);
        model_for.set_y0(y0_for);
        model_rev.set_y0(y0_rev);
        if (i==0) {
            datfile_for << "X0onA" << model_for.get_basename() << ".dat";
            datfile_rev << "X0onA" << model_for.get_basename() << "_rev.dat";
            pngfile_for << "X0onA" << model_for.get_basename() << ".png";
            pngfile_rev << "X0onA" << model_for.get_basename() << "_rev.png";
            pngfile_both << "X0onA" << model_for.get_basename() << "_both.png";
            datafile_for.open(datfile_for.str());
            datafile_rev.open(datfile_rev.str());
        }
        for (long t = 0; t <= t1; t++) {
            if ( (t>=t0) && (t<=t1) ) {
                datafile_for << par_for << " " << t*model_for.getdt();
                datafile_rev << par_rev << " " << t*model_rev.getdt();
                for (long n = 0; n < model_for.varnumber(); n++)
                    if ( (n == 0) ) {
                        datafile_for << " " << model_for.get_vars(n);
                        datafile_rev << " " << model_rev.get_vars(n);
                    }
                    datafile_for << "\n";
                datafile_rev << "\n";
            }
            model_for.step();
            model_rev.step();
        }
        par_for+=da;
        par_rev-=da;
        x0_for = model_for.get_vars(0);
        x0_rev = model_rev.get_vars(0);
        y0_for = model_for.get_vars(1);
        y0_rev = model_rev.get_vars(1);
    }
    datafile_for.close();
    datafile_rev.close();
    std::string parname = "α";
    for (long picn = 0; picn < 3; picn++) {
        std::stringstream plot_command;
        plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel '" << parname << "'; set ylabel 'x'; set output '";
        if (picn==0)
            plot_command << pngfile_both.str();
        if (picn==1)
            plot_command << pngfile_for.str();
        if (picn==2)
            plot_command << pngfile_rev.str();
        plot_command << "'; stat '" << datfile_for.str() << "' u 1:3 nooutput; miny = STATS_min_y; maxy = STATS_max_y; stat '" << datfile_rev.str() << "' u 1:3 nooutput; miny = STATS_min_y<miny?STATS_min_y:miny; maxy = STATS_max_y<maxy?STATS_max_y:maxy; range = maxy-miny; set yrange [miny-0.05*range:maxy+0.05*range]; ";
        if (picn==0)
            plot_command << "plot '" << datfile_for.str() << "' u 1:3:(0.2) w p pt 7 ps variable lc rgb 0x80ff8040 title 'x(increasing " << parname << ")', '" << datfile_rev.str() << "' u 1:3:(0.2) w p pt 7 ps variable lc rgb 0x802e8b57 title 'x(decreasing " << parname << ")';\nEOF";
        if (picn==1)
            plot_command << "plot '" << datfile_for.str() << "' u 1:3:(0.2) w p pt 7 ps variable lc rgb 0x80ff8040 title 'x(increasing " << parname << ")';\nEOF";
        if (picn==2)
            plot_command << "plot '" << datfile_rev.str() << "' u 1:3:(0.2) w p pt 7 ps variable lc rgb 0x802e8b57 title 'x(decreasing " << parname << ")';\nEOF";
        std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
        std::system (plot_command.str().c_str());
    }
}

template <typename specificmodel>
void researcher<specificmodel>::L12onS_s0_ds_s1_x0_y0_alpha_mu_t0_t1 (double s0, double ds, double s1, double x0, double y0, double alpha, double mu, long t0, long t1) {
    std::stringstream datfile_for, datfile_rev, pngfile_for, pngfile_rev, pngfile_both;
    double par_for = s0;
    double par_rev = s1;
    double x0_for = x0;
    double x0_rev = x0;
    double y0_for = y0;
    double y0_rev = y0;
    std::ofstream datafile_for, datafile_rev;
    for (long i = 0; par_for <= s1; i++) {
        specificmodel model_for, model_rev;
        model_for.lyap_init();
        model_rev.lyap_init();
        model_for.set_par_a(alpha);
        model_rev.set_par_a(alpha);
        model_for.set_par_m(mu);
        model_rev.set_par_m(mu);
        model_for.set_par_s(par_for);
        model_rev.set_par_s(par_rev);
        model_for.set_buffer_size(2);
        model_rev.set_buffer_size(2);
        model_for.set_x0(x0_for);
        model_rev.set_x0(x0_rev);
        model_for.set_y0(y0_for);
        model_rev.set_y0(y0_rev);
        if (i==0) {
            datfile_for << "L12onS" << model_for.get_basename() << ".dat";
            datfile_rev << "L12onS" << model_for.get_basename() << "_rev.dat";
            pngfile_for << "L12onS" << model_for.get_basename() << ".png";
            pngfile_rev << "L12onS" << model_for.get_basename() << "_rev.png";
            pngfile_both << "L12onS" << model_for.get_basename() << "_both.png";
            datafile_for.open(datfile_for.str());
            datafile_rev.open(datfile_rev.str());
        }
        model_for.normalization_init();
        model_rev.normalization_init();
        for (long t = 0; t <= t1; t++) {
            model_for.step();
            model_rev.step();
            if (t==t0) {
                model_for.normalization_reset();
                model_rev.normalization_reset();
            }
            else {
                model_for.GramSmidt();
                model_rev.GramSmidt();
            }
        }
        model_for.GramSmidt();
        model_rev.GramSmidt();
        datafile_for << par_for;
        datafile_rev << par_rev;
        for (long l = 0; l < model_for.lyapnumber(); l++) {
            if ( ( l == 0 ) || ( l == 1) ) {
                datafile_for << " " << model_for.l(l);
                datafile_rev << " " << model_rev.l(l);
            }
        }
        datafile_for << "\n";
        datafile_rev << "\n";
        par_for+=ds;
        par_rev-=ds;
        x0_for = model_for.get_vars(0);
        x0_rev = model_rev.get_vars(0);
        y0_for = model_for.get_vars(1);
        y0_rev = model_rev.get_vars(1);
    }
    datafile_for.close();
    datafile_rev.close();
    std::string parname = "σ";
    for (long picn = 0; picn < 3; picn++) {
        std::stringstream plot_command;
        plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel '" << parname << "'; set ylabel 'λ'; set output '";
        if (picn==0)
            plot_command << pngfile_both.str();
        if (picn==1)
            plot_command << pngfile_for.str();
        if (picn==2)
            plot_command << pngfile_rev.str();
        plot_command << "'; ";
        if (picn==0)
            plot_command << "plot '" << datfile_for.str() << "' u 1:2 w l lw 3 lc rgb 0x801400 title 'λ1(increasing " << parname << ")', '' u 1:3 w l lw 3 lc rgb 0xffa500 title 'λ2(increasing " << parname << ")', '' u 1:4 w l lw 3 lc rgb 0xffffff title 'λ3(increasing " << parname << ")', '" << datfile_rev.str() << "' u 1:2 w l dt 2 lw 3 lc rgb 0x801400 title 'λ1(decreasing " << parname << ")', '' u 1:3 w l dt 2 lw 3 lc rgb 0xffa500 title 'λ2(decreasing " << parname << ")', '' u 1:4 w l dt 2 lw 3 lc rgb 0xffffff title 'λ3(decreasing " << parname << ")', 0 w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        if (picn==1)
            plot_command << "plot '" << datfile_for.str() << "' u 1:2 w l lw 3 lc rgb 0x801400 title 'λ1(increasing " << parname << ")', '' u 1:3 w l lw 3 lc rgb 0xffa500 title 'λ2(increasing " << parname << ")', '' u 1:4 w l lw 3 lc rgb 0xffffff title 'λ3(increasing " << parname << ")', 0 w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        if (picn==2)
            plot_command << "plot '" << datfile_rev.str() << "' u 1:2 w l lw 3 lc rgb 0x801400 title 'λ1(decreasing " << parname << ")', '' u 1:3 w l lw 3 lc rgb 0xffa500 title 'λ2(decreasing " << parname << ")', '' u 1:4 w l lw 3 lc rgb 0xffffff title 'λ3(decreasing " << parname << ")', 0 w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
        std::system (plot_command.str().c_str());
    }
}

template <typename specificmodel>
void researcher<specificmodel>::L12onA_a0_da_a1_x0_y0_mu_sigma_t0_t1 (double a0, double da, double a1, double x0, double y0, double mu, double sigma, long t0, long t1) {
    std::stringstream datfile_for, datfile_rev, pngfile_for, pngfile_rev, pngfile_both;
    double par_for = a0;
    double par_rev = a1;
    double x0_for = x0;
    double x0_rev = x0;
    double y0_for = y0;
    double y0_rev = y0;
    std::ofstream datafile_for, datafile_rev;
    for (long i = 0; par_for <= a1; i++) {
        specificmodel model_for, model_rev;
        model_for.lyap_init();
        model_rev.lyap_init();
        model_for.set_par_a(par_for);
        model_rev.set_par_a(par_rev);
        model_for.set_par_m(mu);
        model_rev.set_par_m(mu);
        model_for.set_par_s(sigma);
        model_rev.set_par_s(sigma);
        model_for.set_buffer_size(2);
        model_rev.set_buffer_size(2);
        model_for.set_x0(x0_for);
        model_rev.set_x0(x0_rev);
        model_for.set_y0(y0_for);
        model_rev.set_y0(y0_rev);
        if (i==0) {
            datfile_for << "L12onA" << model_for.get_basename() << ".dat";
            datfile_rev << "L12onA" << model_for.get_basename() << "_rev.dat";
            pngfile_for << "L12onA" << model_for.get_basename() << ".png";
            pngfile_rev << "L12onA" << model_for.get_basename() << "_rev.png";
            pngfile_both << "L12onA" << model_for.get_basename() << "_both.png";
            datafile_for.open(datfile_for.str());
            datafile_rev.open(datfile_rev.str());
        }
        model_for.normalization_init();
        model_rev.normalization_init();
        for (long t = 0; t <= t1; t++) {
            model_for.step();
            model_rev.step();
            if (t==t0) {
                model_for.normalization_reset();
                model_rev.normalization_reset();
            }
            else {
                model_for.GramSmidt();
                model_rev.GramSmidt();
            }
        }
        model_for.GramSmidt();
        model_rev.GramSmidt();
        datafile_for << par_for;
        datafile_rev << par_rev;
        for (long l = 0; l < model_for.lyapnumber(); l++) {
            if ( ( l == 0 ) || ( l == 1) ) {
                datafile_for << " " << model_for.l(l);
                datafile_rev << " " << model_rev.l(l);
            }
        }
        datafile_for << "\n";
        datafile_rev << "\n";
        par_for+=da;
        par_rev-=da;
        x0_for = model_for.get_vars(0);
        x0_rev = model_rev.get_vars(0);
        y0_for = model_for.get_vars(1);
        y0_rev = model_rev.get_vars(1);
    }
    datafile_for.close();
    datafile_rev.close();
    std::string parname = "α";
    for (long picn = 0; picn < 3; picn++) {
        std::stringstream plot_command;
        plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel '" << parname << "'; set ylabel 'λ'; set output '";
        if (picn==0)
            plot_command << pngfile_both.str();
        if (picn==1)
            plot_command << pngfile_for.str();
        if (picn==2)
            plot_command << pngfile_rev.str();
        plot_command << "'; ";
        if (picn==0)
            plot_command << "plot '" << datfile_for.str() << "' u 1:2 w l lw 3 lc rgb 0x801400 title 'λ1(increasing " << parname << ")', '' u 1:3 w l lw 3 lc rgb 0xffa500 title 'λ2(increasing " << parname << ")', '' u 1:4 w l lw 3 lc rgb 0xffffff title 'λ3(increasing " << parname << ")', '" << datfile_rev.str() << "' u 1:2 w l dt 2 lw 3 lc rgb 0x801400 title 'λ1(decreasing " << parname << ")', '' u 1:3 w l dt 2 lw 3 lc rgb 0xffa500 title 'λ2(decreasing " << parname << ")', '' u 1:4 w l dt 2 lw 3 lc rgb 0xffffff title 'λ3(decreasing " << parname << ")', 0 w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        if (picn==1)
            plot_command << "plot '" << datfile_for.str() << "' u 1:2 w l lw 3 lc rgb 0x801400 title 'λ1(increasing " << parname << ")', '' u 1:3 w l lw 3 lc rgb 0xffa500 title 'λ2(increasing " << parname << ")', '' u 1:4 w l lw 3 lc rgb 0xffffff title 'λ3(increasing " << parname << ")', 0 w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        if (picn==2)
            plot_command << "plot '" << datfile_rev.str() << "' u 1:2 w l lw 3 lc rgb 0x801400 title 'λ1(decreasing " << parname << ")', '' u 1:3 w l lw 3 lc rgb 0xffa500 title 'λ2(decreasing " << parname << ")', '' u 1:4 w l lw 3 lc rgb 0xffffff title 'λ3(decreasing " << parname << ")', 0 w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
        std::system (plot_command.str().c_str());
    }
}

template <typename specificmodel>
void researcher<specificmodel>::ExonS_s0_ds_s1_x0_y0_alpha_mu_t0_t1_winsize (double s0, double ds, double s1, double x0, double y0, double alpha, double mu, long t0, long t1, long ws) {
    std::stringstream datfile_for, datfile_rev, pngfile_for, pngfile_rev, pngfile_both;
    double par_for = s0;
    double par_rev = s1;
    double x0_for = x0;
    double x0_rev = x0;
    double y0_for = y0;
    double y0_rev = y0;
    std::ofstream datafile_for, datafile_rev;
    for (long i = 0; par_for <= s1; i++) {
        specificmodel model_for, model_rev;
        model_for.set_par_a(alpha);
        model_rev.set_par_a(alpha);
        model_for.set_par_m(mu);
        model_rev.set_par_m(mu);
        model_for.set_par_s(par_for);
        model_rev.set_par_s(par_rev);
        model_for.set_buffer_size(2);
        model_rev.set_buffer_size(2);
        model_for.set_x0(x0_for);
        model_rev.set_x0(x0_rev);
        model_for.set_y0(y0_for);
        model_rev.set_y0(y0_rev);
        if (i==0) {
            datfile_for << "ExonSv" << model_for.get_basename() << ".dat";
            datfile_rev << "ExonSv" << model_for.get_basename() << "_rev.dat";
            pngfile_for << "ExonSv" << model_for.get_basename() << ".png";
            pngfile_rev << "ExonSv" << model_for.get_basename() << "_rev.png";
            pngfile_both << "ExonSv" << model_for.get_basename() << "_both.png";
            datafile_for.open(datfile_for.str());
            datafile_rev.open(datfile_rev.str());
        }
        for (long n = 0; n < model_for.varnumber(); n++) {
            if ( ( n == 0 ) ) {
                model_for.prob_init_vars(n, 501);
                model_rev.prob_init_vars(n, 501);
            }
        }
        for (long t = 0; t <= t1; t++) {
            model_for.step();
            model_rev.step();
            for (long n = 0; n < model_for.varnumber(); n++) {
                if ( ( n == 0 ) ) {
                    model_for.update_stat_vars(n);
                    model_rev.update_stat_vars(n);
                }
            }
        }
        datafile_for << par_for;
        datafile_rev << par_rev;
        for (long n = 0; n < model_for.varnumber(); n++) {
            if ( ( n == 0 ) ) {
                datafile_for << " " << model_for.get_entropy_vars(n);
                datafile_rev << " " << model_rev.get_entropy_vars(n);
            }
        }
        datafile_for << "\n";
        datafile_rev << "\n";
        par_for+=ds;
        par_rev-=ds;
        x0_for = model_for.get_vars(0);
        x0_rev = model_rev.get_vars(0);
        y0_for = model_for.get_vars(1);
        y0_rev = model_rev.get_vars(1);
    }
    datafile_for.close();
    datafile_rev.close();
    std::string parname = "σ";
    for (long picn = 0; picn < 3; picn++) {
        std::stringstream plot_command;
        plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel '" << parname << "'; set ylabel 'Entropy power'; set output '";
        if (picn==0)
            plot_command << pngfile_both.str();
        if (picn==1)
            plot_command << pngfile_for.str();
        if (picn==2)
            plot_command << pngfile_rev.str();
        plot_command << "'; ";
        if (picn==0)
            plot_command << "plot '" << datfile_for.str() << "' u 1:(exp(\\$2)) w l lw 3 lc rgb 0x408000 title 'Entropy power(incr. " << parname << ")', '" << datfile_rev.str() << "' u 1:(exp(\\$2)) w l dt 2 lw 3 lc rgb 0x408000 title 'Entropy power(decr. " << parname << ")', exp(1.4189) w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        if (picn==1)
            plot_command << "plot '" << datfile_for.str() << "' u 1:(exp(\\$2)) w l lw 3 lc rgb 0x408000 title 'Entropy power(incr. " << parname << ")', exp(1.4189) w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        if (picn==2)
            plot_command << "plot '" << datfile_rev.str() << "' u 1:(exp(\\$2)) w l lw 3 lc rgb 0x408000 title 'Entropy power(decr. " << parname << ")', exp(1.4189) w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
        std::system (plot_command.str().c_str());
    }
}

template <typename specificmodel>
void researcher<specificmodel>::ExonA_a0_da_a1_x0_y0_mu_sigma_t0_t1_winsize (double a0, double da, double a1, double x0, double y0, double mu, double sigma, long t0, long t1, long ws) {
    std::stringstream datfile_for, datfile_rev, pngfile_for, pngfile_rev, pngfile_both;
    double par_for = a0;
    double par_rev = a1;
    double x0_for = x0;
    double x0_rev = x0;
    double y0_for = y0;
    double y0_rev = y0;
    std::ofstream datafile_for, datafile_rev;
    for (long i = 0; par_for <= a1; i++) {
        specificmodel model_for, model_rev;
        model_for.set_par_a(par_for);
        model_rev.set_par_a(par_rev);
        model_for.set_par_m(mu);
        model_rev.set_par_m(mu);
        model_for.set_par_s(sigma);
        model_rev.set_par_s(sigma);
        model_for.set_buffer_size(2);
        model_rev.set_buffer_size(2);
        model_for.set_x0(x0_for);
        model_rev.set_x0(x0_rev);
        model_for.set_y0(y0_for);
        model_rev.set_y0(y0_rev);
        if (i==0) {
            datfile_for << "ExonAv" << model_for.get_basename() << ".dat";
            datfile_rev << "ExonAv" << model_for.get_basename() << "_rev.dat";
            pngfile_for << "ExonAv" << model_for.get_basename() << ".png";
            pngfile_rev << "ExonAv" << model_for.get_basename() << "_rev.png";
            pngfile_both << "ExonAv" << model_for.get_basename() << "_both.png";
            datafile_for.open(datfile_for.str());
            datafile_rev.open(datfile_rev.str());
        }
        for (long n = 0; n < model_for.varnumber(); n++) {
            if ( ( n == 0 ) ) {
                model_for.prob_init_vars(n, 501);
                model_rev.prob_init_vars(n, 501);
            }
        }
        for (long t = 0; t <= t1; t++) {
            model_for.step();
            model_rev.step();
            for (long n = 0; n < model_for.varnumber(); n++) {
                if ( ( n == 0 ) ) {
                    model_for.update_stat_vars(n);
                    model_rev.update_stat_vars(n);
                }
            }
        }
        datafile_for << par_for;
        datafile_rev << par_rev;
        for (long n = 0; n < model_for.varnumber(); n++) {
            if ( ( n == 0 ) ) {
                datafile_for << " " << model_for.get_entropy_vars(n);
                datafile_rev << " " << model_rev.get_entropy_vars(n);
            }
        }
        datafile_for << "\n";
        datafile_rev << "\n";
        par_for+=da;
        par_rev-=da;
        x0_for = model_for.get_vars(0);
        x0_rev = model_rev.get_vars(0);
        y0_for = model_for.get_vars(1);
        y0_rev = model_rev.get_vars(1);
    }
    datafile_for.close();
    datafile_rev.close();
    std::string parname = "α";
    for (long picn = 0; picn < 3; picn++) {
        std::stringstream plot_command;
        plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel '" << parname << "'; set ylabel 'Entropy power'; set output '";
        if (picn==0)
            plot_command << pngfile_both.str();
        if (picn==1)
            plot_command << pngfile_for.str();
        if (picn==2)
            plot_command << pngfile_rev.str();
        plot_command << "'; ";
        if (picn==0)
            plot_command << "plot '" << datfile_for.str() << "' u 1:(exp(\\$2)) w l lw 3 lc rgb 0x408000 title 'Entropy power(incr. " << parname << ")', '" << datfile_rev.str() << "' u 1:(exp(\\$2)) w l dt 2 lw 3 lc rgb 0x408000 title 'Entropy power(decr. " << parname << ")', exp(1.4189) w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        if (picn==1)
            plot_command << "plot '" << datfile_for.str() << "' u 1:(exp(\\$2)) w l lw 3 lc rgb 0x408000 title 'Entropy power(incr. " << parname << ")', exp(1.4189) w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        if (picn==2)
            plot_command << "plot '" << datfile_rev.str() << "' u 1:(exp(\\$2)) w l lw 3 lc rgb 0x408000 title 'Entropy power(decr. " << parname << ")', exp(1.4189) w l dt 2 lw 1 lc rgb 0x000000 title '';\nEOF";
        std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
        std::system (plot_command.str().c_str());
    }
}

template <typename specificmodel>
void researcher<specificmodel>::Min0x_x0_y0_alpha_mu_sigma_t0_t1 (double x0, double y0, double alpha, double mu, double sigma, long t0, long t1) {
    specificmodel model;
    model.set_par_a(alpha);
    model.set_par_m(mu);
    model.set_par_s(sigma);
    model.set_buffer_size(2);
    model.set_x0(x0);
    model.set_y0(y0);
    std::stringstream datfile, sttfile, datfileMIN, datfileMINt, pngfile;
    datfile << "TS0x" << model.get_basename() << ".dat";
    sttfile << "STAT0x" << model.get_basename() << ".dat";
    datfileMIN << "TS0xMIN" << model.get_basename() << ".dat";
    datfileMINt << "TS0xMINt" << model.get_basename() << ".dat";
    pngfile << "TS0x" << model.get_basename() << ".png";
    std::ofstream datafile(datfile.str());
    std::ofstream datafileMIN(datfileMIN.str());
    std::ofstream datafileMINt(datfileMINt.str());
    long every = (t1-t0)/100000;
    if (every==0) every = 1;
    for (long i = 0; i <= t1; i++) {
        model.step();
        if (i==t0) {
            model.init_minima(0, -1.5);
            model.init_minimaRR(0, -1.5);
        }
        if (i>=t0) {
            long var = 0; bool cs; double ct; double cx;
            model.get_minima(0, cs, ct, cx);
            if (cs) datafileMIN << ct + t0 << " " << cx << "\n";
            model.get_minimaRR(0, cs, ct, cx);
            if (cs) datafileMINt << ct + t0 << " " << cx << "\n";
            if ((i % every) == 0) datafile << model.get_time(0) << " " << model.get_vars(0) << "\n";
        }
    }
    datafile.close();
    datafileMIN.close();
    datafileMINt.close();
    std::ofstream statfile(sttfile.str());
    statfile << "meanA    " << model.stat_get_mean_minima(0) << "\n";
    statfile << "sigmaA   " << model.stat_get_disp_root_minima(0) << "\n";
    statfile << "skewA    " << model.stat_get_skew_minima(0) << "\n";
    statfile << "kurtA    " << model.stat_get_kurt_minima(0) << "\n";
    statfile << "entropyA " << model.get_entropy_minima(0) << "\n";
    statfile << "outof3sA " << model.extreme_part_outside_minima(0, 3.0) << "\n";
    statfile << "outof5sA " << model.extreme_part_outside_minima(0, 5.0) << "\n";
    statfile << "outof7sA " << model.extreme_part_outside_minima(0, 7.0) << "\n";
    statfile << "meanT    " << model.stat_get_mean_minimaRR(0) << "\n";
    statfile << "sigmaT   " << model.stat_get_disp_root_minimaRR(0) << "\n";
    statfile << "skewT    " << model.stat_get_skew_minimaRR(0) << "\n";
    statfile << "kurtT    " << model.stat_get_kurt_minimaRR(0) << "\n";
    statfile << "entropyT " << model.get_entropy_minimaRR(0) << "\n";
    statfile << "outof3sT " << model.extreme_part_outside_minimaRR(0, 3.0) << "\n";
    statfile << "outof5sT " << model.extreme_part_outside_minimaRR(0, 5.0) << "\n";
    statfile << "outof7sT " << model.extreme_part_outside_minimaRR(0, 7.0) << "\n";
    std::stringstream plot_command;
    plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 't'; set ylabel 'x, y'; set output '" << pngfile.str() << "'; stat '" << datfile.str() << "' u 1:2 nooutput; range = STATS_max_y-STATS_min_y; set yrange [STATS_min_y - 0.05 * range : STATS_max_y + 0.05 * range]; ";
    plot_command << "plot '" << datfile.str() << "' u 1:2 w l lw 3 lc rgb 0x008040 title 'x(t)', '" << datfileMIN.str() << "' u 1:2 w p pt 5 ps 1 lc rgb 'black' title 'minima';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
    std::stringstream probAdat, probApng;
    probAdat << "Prob0xAmin" << model.get_basename() << ".dat";
    probApng << "Prob0xAmin" << model.get_basename() << ".png";
    std::ofstream datafileA(probAdat.str());
    for (long i = 0; i < model.get_bins_number_minima(0); i++) {
        double tval = model.get_prob_dens_minima(0, i);
        if (tval>0.0) {
            double mid = 0.5 * (model.get_prob_right_minima(0, i) + model.get_prob_left_minima(0, i));
            double width = model.get_prob_right_minima(0, i) - model.get_prob_left_minima(0, i);
            datafileA << mid << " " << tval << " " << width << "\n";
        }
    }
    datafileA.close();
    std::stringstream plot_command2;
    plot_command2 << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 'x'; set ylabel 'p(x)'; set output '" << probApng.str() << "'; stat '" << probAdat.str() << "' u 1:2 nooutput; maxval = STATS_max_y; set yrange[0:maxval*1.1]; minval = STATS_min_x; maxval = STATS_max_x; range = maxval-minval+0.01; set xrange[minval-0.05*range:maxval+0.05*range]; set label 'entropy = " << model.get_entropy2_minima(0) << " bits' at screen 0.7,0.95; plot '" << probAdat.str() << "' u 1:2:3 w boxes lc rgb 0x00556b2f title 'p(x)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command2.str() << "\n\n";
    std::system (plot_command2.str().c_str());
    std::stringstream probTdat, probTpng;
    probTdat << "Prob0xTmin" << model.get_basename() << ".dat";
    probTpng << "Prob0xTmin" << model.get_basename() << ".png";
    std::ofstream datafileT(probTdat.str());
    for (long i = 0; i < model.get_bins_number_minimaRR(0); i++) {
        double tval = model.get_prob_dens_minimaRR(0, i);
        if (tval>0.0) {
            double mid = 0.5 * (model.get_prob_right_minimaRR(0, i) + model.get_prob_left_minimaRR(0, i));
            double width = model.get_prob_right_minimaRR(0, i) - model.get_prob_left_minimaRR(0, i);
            datafileT << mid << " " << tval << " " << width << "\n";
        }
    }
    datafileT.close();
    std::stringstream plot_command3;
    plot_command3 << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal; set xlabel 'x'; set ylabel 'p(x)'; set output '" << probTpng.str() << "'; stat '" << probTdat.str() << "' u 1:2 nooutput; maxval = STATS_max_y; set yrange[0:maxval*1.1]; minval = STATS_min_x; maxval = STATS_max_x; range = maxval-minval+0.01; set xrange[minval-0.05*range:maxval+0.05*range]; set label 'entropy = " << model.get_entropy2_minimaRR(0) << " bits' at screen 0.7,0.95; plot '" << probTdat.str() << "' u 1:2:3 w boxes lc rgb 0x00556b2f title 'p(x)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command3.str() << "\n\n";
    std::system (plot_command3.str().c_str());
}

template <typename specificmodel>
void researcher<specificmodel>::LevySourceStat_alpha_beta_sigma_mu_idum_N (double alpha, double beta, double sigma, double mu, long idum, long N) {
    specificmodel model;
    model.set_pars(2, mu, alpha, beta, sigma, idum);
    model.stat_init();
    model.corr_init();
    model.prob_init(501);
    std::stringstream sttfile, pdtfile, pngfile;
    pdtfile << "STAT_Levy_prob" << model.get_basename() << ".dat";
    sttfile << "STAT_Levy_stat" << model.get_basename() << ".dat";
    pngfile << "STAT_Levy_prob" << model.get_basename() << ".png";
    for (long s=0;s<N;s++) {
        model.step();
        model.update_stat();
    }
    std::ofstream statfile(sttfile.str());
    statfile << "meanA    " << model.stat_get_mean() << "\n";
    statfile << "sigmaA   " << model.stat_get_disp_root() << "\n";
    statfile << "skewA    " << model.stat_get_skew() << "\n";
    statfile << "kurtA    " << model.stat_get_kurt() << "\n";
    statfile << "entropyA " << model.get_entropy() << "\n";
    statfile << "outof3sA " << model.extreme_part_outside(3.0) << "\n";
    statfile << "outof5sA " << model.extreme_part_outside(5.0) << "\n";
    statfile << "outof7sA " << model.extreme_part_outside(7.0) << "\n";
    std::ofstream pdatafile(pdtfile.str());
    for (long i = 0; i < model.get_bins_number(); i++) {
        double tval = model.get_prob_dens(i);
        if (tval>0.0) {
            double mid = 0.5 * (model.get_prob_right(i) + model.get_prob_left(i));
            double width = model.get_prob_right(i) - model.get_prob_left(i);
            pdatafile << mid << " " << tval << " " << width << "\n";
        }
    }
    pdatafile.close();
    std::stringstream plot_command2;
    plot_command2 << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal;"
                  << " set xlabel 'η'; set ylabel 'p(η)'; set output '" << pngfile.str() << "'; stat '" << pdtfile.str() << "' u 1:2 nooutput; maxval = STATS_max_y; set yrange[0:maxval*1.1]; minval = STATS_min_x; maxval = STATS_max_x; range = maxval-minval+0.01; set xrange[minval-0.05*range:maxval+0.05*range]; set label 'entropy = " << model.get_entropy2() << " bits' at screen 0.7,0.95; plot '" << pdtfile.str() << "' u 1:2:3 w boxes lc rgb 0x00556b2f title 'p(x)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command2.str() << "\n\n";
    std::system (plot_command2.str().c_str());
}

template <typename specificmodel>
void researcher<specificmodel>::PoissonSourceStat_alpha_beta_sigma_mu_idum_N (double freq, long idum, long N) {
    specificmodel model;
    model.set_pars(2, freq, idum);
    model.init_events();
    std::stringstream sttfile, pdtfile, pngfile, datTSfile, pngTSfile;
    pdtfile << "STAT_Poisson_prob" << model.get_basename() << ".dat";
    sttfile << "STAT_Poisson_stat" << model.get_basename() << ".dat";
    pngfile << "STAT_Poisson_prob" << model.get_basename() << ".png";
    datTSfile << "TS_Poisson" << model.get_basename() << ".dat";
    pngTSfile << "TS_Poisson" << model.get_basename() << ".png";
    std::ofstream dataTSfile(datTSfile.str());
    double start_time = -1.0e6;
    double finit_time = -1.0e6;
    long points_number = 0;
    for (long s=0;s<N;s++) {
        model.step();
        bool cs; double ct; double cx;
        model.get_events(cs, ct, cx);
        if (cs) {
            if (points_number<1e4)
                dataTSfile << ct << " " << 0.0 << "\n" << ct << " " << 1.0 << "\n" << ct << " " << 0.0 << "\n";
            if (start_time<0) start_time = ct;
            points_number++;
            if (points_number==80) finit_time = ct;
        }
    }
    dataTSfile.close();
    std::stringstream plot_command;
    plot_command << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal;"
                 << " set xlabel 't'; set ylabel 'φ(t)'; set output '" << pngTSfile.str() << "'; set yrange [-0.1 : 1.1]; set xrange [" << start_time*0.5 << " : " << finit_time << "]; ";
    plot_command << "plot '" << datTSfile.str() << "' u 1:2 w l lw 3 lc rgb 'black' title 'Poisson process';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command.str() << "\n\n";
    std::system (plot_command.str().c_str());
    std::ofstream statfile(sttfile.str());
    statfile << "meanA    " << model.stat_get_mean_events() << "\n";
    statfile << "sigmaA   " << model.stat_get_disp_root_events() << "\n";
    statfile << "skewA    " << model.stat_get_skew_events() << "\n";
    statfile << "kurtA    " << model.stat_get_kurt_events() << "\n";
    statfile << "entropyA " << model.get_entropy_events() << "\n";
    statfile << "outof3sA " << model.extreme_part_outside_events(3.0) << "\n";
    statfile << "outof5sA " << model.extreme_part_outside_events(5.0) << "\n";
    statfile << "outof7sA " << model.extreme_part_outside_events(7.0) << "\n";
    std::ofstream pdatafile(pdtfile.str());
    for (long i = 0; i < model.get_bins_number_events(); i++) {
        double tval = model.get_prob_dens_events(i);
        if (tval>0.0) {
            double mid = 0.5 * (model.get_prob_right_events(i) + model.get_prob_left_events(i));
            double width = model.get_prob_right_events(i) - model.get_prob_left_events(i);
            pdatafile << mid << " " << tval << " " << width << "\n";
        }
    }
    pdatafile.close();
    std::stringstream plot_command2;
    plot_command2 << "gnuplot << EOF\nset terminal pngcairo size 1200,600 enhanced font 'Verdana,20'; unset warnings; set key tmargin center horizontal;"
                  << " set xlabel 'φ'; set ylabel 'p(φ)'; set output '" << pngfile.str() << "'; stat '" << pdtfile.str() << "' u 1:2 nooutput; maxval = STATS_max_y; set yrange[0:maxval*1.1]; minval = STATS_min_x; maxval = STATS_max_x; range = maxval-minval+0.01; set xrange[minval-0.05*range:maxval+0.05*range]; set label 'entropy = " << model.get_entropy2_events() << " bits' at screen 0.7,0.95; plot '" << pdtfile.str() << "' u 1:2:3 w boxes lc rgb 0x00556b2f title 'p(x)';\nEOF";
    std::cout << "\nGnuplot:\n" << plot_command2.str() << "\n\n";
    std::system (plot_command2.str().c_str());
}

#endif