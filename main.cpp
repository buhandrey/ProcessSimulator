#include "researcher.h"
#include "Rulkov_map.h"
#include "levy.h"
#include "poisson.h"

int main() {
    // researcher<Rulkov_map> rkm;
    // rkm.TS0xy_x0_y0_alpha_mu_sigma_t0_t1 (1.0, -2.9, 4.0, 0.001, 1.0, 99e4, 1e6);
    // rkm.TS0xyLCE_x0_y0_alpha_mu_sigma_t0_t1 (1.0, -2.9, 4.0, 0.001, 1.0, 99e4, 1e6);
    // rkm.PP0xy_x0_y0_alpha_mu_sigma_t0_t1 (1.0, -2.9, 4.0, 0.001, 1.0, 99e4, 1e6);
    // rkm.PP0xyLCE_x0_y0_alpha_mu_sigma_t0_t1 (1.0, -2.9, 4.0, 0.001, 1.0, 99e4, 1e6);
    // rkm.PDF0x_x0_y0_alpha_mu_sigma_t0_t1_winsize (1.0, -2.9, 4.0, 0.001, 1.0, 3e5, 1e6, 501);
    // rkm.X0onS_s0_ds_s1_x0_y0_alpha_mu_t0_t1 (-2.0, 0.004, 2.0, 1.0, -2.9, 4.0, 0.001, 999800, 1e6);
    // rkm.X0onA_a0_da_a1_x0_y0_mu_sigma_t0_t1 (0.0, 0.006, 6.0, 1.0, -2.9, 0.001, 1.0, 999800, 1e6);
    // rkm.L12onS_s0_ds_s1_x0_y0_alpha_mu_t0_t1 (-2.0, 0.004, 2.0, 1.0, -2.9, 4.0, 0.001, 3e5, 1e6);
    // rkm.L12onA_a0_da_a1_x0_y0_mu_sigma_t0_t1 (0.0, 0.006, 6.0, 1.0, -2.9, 0.001, 1.0, 3e5, 1e6);
    // rkm.ExonS_s0_ds_s1_x0_y0_alpha_mu_t0_t1_winsize (-2.0, 0.004, 2.0, 1.0, -2.9, 4.0, 0.001, 3e5, 1e6, 501);
    // rkm.ExonA_a0_da_a1_x0_y0_mu_sigma_t0_t1_winsize (0.0, 0.006, 6.0, 1.0, -2.9, 0.001, 1.0, 3e5, 1e6, 501);
    // rkm.Min0x_x0_y0_alpha_mu_sigma_t0_t1(1.0, -2.9, 4.0, 0.001, 1.0, 3e5, 1e6);
    // researcher<levy> levynoise;
    // levynoise.LevySourceStat_alpha_beta_sigma_mu_idum_N (1.9999, -0.5, 1.0, 0.0, -341, 1e7);
    researcher<poisson> poissonnoise;
    poissonnoise.PoissonSourceStat_alpha_beta_sigma_mu_idum_N (0.001, -341, 1e5);
    return 0;
}