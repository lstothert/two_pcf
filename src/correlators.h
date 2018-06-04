#pragma once
#include "boxes.h"
#include "galaxy.h"
#include "hist.h"
#include "io.h"

void auto_pair_count_monopole(const GalaxyBoxes *cat, JKHists *hist);

void cross_pair_count_monopole(const GalaxyBoxes *cat_1, const GalaxyBoxes *cat_2, JKHists *hist);

void auto_correlate_monopole(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params);

void bin_pair_sigma_pi(const Galaxy *gal_1, const Galaxy *gal_2, JKHists2D *hist, const float min_r_sq, const float max_r_sq);

void auto_pair_count_sigma_pi(const GalaxyBoxes *cat, JKHists2D *hist);

void cross_pair_count_sigma_pi(const GalaxyBoxes *cat_1, const GalaxyBoxes *cat_2, JKHists2D *hist);

void auto_correlate_sigma_pi(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params);

void bin_pair_s_mu(const Galaxy *gal_1, const Galaxy *gal_2, JKHists2D *hist, const float min_r_sq, const float max_r_sq);

void auto_pair_count_s_mu(const GalaxyBoxes *cat, JKHists2D *hist);

void cross_pair_count_s_mu(const GalaxyBoxes *cat_1, const GalaxyBoxes *cat_2, JKHists2D *hist);

void auto_correlate_s_mu(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params);

void auto_pair_count_angular(const GalaxyPixels *cat, JKHists *hist);

void cross_pair_count_angular(const GalaxyPixels *cat_1, const GalaxyPixels *cat_2, JKHists *hist);