#pragma once
#include "boxes.h"
#include "galaxy.h"
#include "hist.h"
#include "io.h"

#ifdef _USE_INV_WEIGHTS

void auto_pair_count_monopole_invpweights(const GalaxyBoxes *cat, JKHists *hist, Params params);

void cross_pair_count_monopole_invpweights(const GalaxyBoxes *cat_1, const GalaxyBoxes *cat_2, JKHists *hist, Params params);

void auto_correlate_monopole_invpweights(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params);

void bin_pair_sigma_pi_invpweights(const Galaxy *gal_1, const Galaxy *gal_2, JKHists2D *hist, const float min_r_sq, const float max_r_sq, JKHists *ang_upweight, Params params);

void auto_pair_count_sigma_pi_invpweights(const GalaxyBoxes *cat, JKHists2D *hist, Params params);

void cross_pair_count_sigma_pi_invpweights(const GalaxyBoxes *cat_1, const GalaxyBoxes *cat_2, JKHists2D *hist, Params params);

void auto_correlate_sigma_pi_invpweights(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params);

void bin_pair_s_mu_invpweights(const Galaxy *gal_1, const Galaxy *gal_2, JKHists2D *hist, const float min_r_sq, const float max_r_sq, JKHists *ang_upweight, Params params);

void auto_pair_count_s_mu_invpweights(const GalaxyBoxes *cat, JKHists2D *hist, Params params);

void cross_pair_count_s_mu_invpweights(const GalaxyBoxes *cat_1, const GalaxyBoxes *cat_2, JKHists2D *hist, Params params);

void auto_correlate_s_mu_invpweights(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params);

void auto_pair_count_angular_invpweights(const GalaxyPixels *cat, JKHists *hist, Params params);

void cross_pair_count_angular_invpweights(const GalaxyPixels *cat_1, const GalaxyPixels *cat_2, JKHists *hist);

void calc_angular_upweight(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params);

#endif