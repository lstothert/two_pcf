#include "hist.h"
#include "galaxy.h"
#include "boxes.h"
#include "io.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "healpix/pointing.h"
#include "healpix/rangeset.h"

void auto_pair_count_monopole(const GalaxyBoxes *cat, JKHists *hist) {

	float max_r_sq = hist->full_hist.ax.max * hist->full_hist.ax.max;
	float min_r_sq = hist->full_hist.ax.min * hist->full_hist.ax.min;

#pragma omp parallel
	{

		JKHists priv_hist(hist->n_jk_regions, hist->full_hist);

#pragma omp for schedule(dynamic)
		for (int ii = 0; ii < cat->n_boxes_used; ++ii) {
			for (int adj_id = 0; adj_id < 27; ++adj_id) {

				int jj = adj_id_to_box_index(adj_id, cat->boxes[ii].loc3d, cat->id_to_box_index, cat->n_boxes_3d);
				if (ii == jj) {

					for (int kk = 0; kk < cat->boxes[ii].numgals; ++kk) {
						for (int ll = kk + 1; ll < cat->boxes[jj].numgals; ++ll) {

							Galaxy *gal_1 = &cat->boxes[ii].gal[kk], *gal_2 = &cat->boxes[jj].gal[ll];

							float diff_vec[3];
							for (int mm = 0; mm < 3; ++mm) {
								diff_vec[mm] = gal_1->x[mm] - gal_2->x[mm];
							}
							float s = diff_vec[0] * diff_vec[0] + diff_vec[1] * diff_vec[1] + diff_vec[2] * diff_vec[2];

							if (s < max_r_sq && s > min_r_sq) {
								priv_hist.add_pair(gal_1->jk_region, gal_2->jk_region, sqrtf(s), gal_1->weight*gal_2->weight);
							}

						}
					}

				}
				else if (jj > ii) {

					for (int kk = 0; kk < cat->boxes[ii].numgals; ++kk) {
						for (int ll = 0; ll < cat->boxes[jj].numgals; ++ll) {

							Galaxy *gal_1 = &cat->boxes[ii].gal[kk], *gal_2 = &cat->boxes[jj].gal[ll];

							float diff_vec[3];
							for (int mm = 0; mm < 3; ++mm) {
								diff_vec[mm] = gal_1->x[mm] - gal_2->x[mm];
							}
							float s = diff_vec[0] * diff_vec[0] + diff_vec[1] * diff_vec[1] + diff_vec[2] * diff_vec[2];

							if (s < max_r_sq && s > min_r_sq) {
								priv_hist.add_pair(gal_1->jk_region, gal_2->jk_region, sqrtf(s), gal_1->weight*gal_2->weight);
							}

						}

					}
				}

			}

		}

		hist->reduce(priv_hist);
	}

	hist->invert_jk_hists();
	hist->multiply_hists(2.0);
}

void cross_pair_count_monopole(const GalaxyBoxes *cat_1, const GalaxyBoxes *cat_2, JKHists *hist) {

	float max_r_sq = hist->full_hist.ax.max * hist->full_hist.ax.max;
	float min_r_sq = hist->full_hist.ax.min * hist->full_hist.ax.min;

#pragma omp parallel 
	{

		JKHists priv_hist(hist->n_jk_regions, hist->full_hist);

#pragma omp for schedule(dynamic)
		for (int ii = 0; ii < cat_1->n_boxes_used; ++ii) {
			for (int adj_id = 0; adj_id < 27; ++adj_id) {

				int jj = adj_id_to_box_index(adj_id, cat_1->boxes[ii].loc3d, cat_2->id_to_box_index, cat_2->n_boxes_3d);
				if (jj != -1) {

					for (int kk = 0; kk < cat_1->boxes[ii].numgals; ++kk) {
						for (int ll = 0; ll < cat_2->boxes[jj].numgals; ++ll) {

							Galaxy *gal_1 = &cat_1->boxes[ii].gal[kk], *gal_2 = &cat_2->boxes[jj].gal[ll];

							float diff_vec[3];
							for (int mm = 0; mm < 3; ++mm) {
								diff_vec[mm] = gal_1->x[mm] - gal_2->x[mm];
							}
							float s = diff_vec[0] * diff_vec[0] + diff_vec[1] * diff_vec[1] + diff_vec[2] * diff_vec[2];

							if (s < max_r_sq && s > min_r_sq) {
								priv_hist.add_pair(gal_1->jk_region, gal_2->jk_region, sqrtf(s), gal_1->weight*gal_2->weight);
							}

						}
					}

				}

			}
		}

		(*hist).reduce(priv_hist);

	}

	hist->invert_jk_hists();
}

void auto_correlate_monopole(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params) {

	std::cout << "Calculating autocorrelation monopole\n";

	Hist temp_hist(params.monopole_min, params.monopole_max, params.monopole_n_bins, params.monopole_log_base);
	JKHists DD(params.n_jk_regions, temp_hist);
	JKHists DR(params.n_jk_regions, temp_hist);
	JKHists RR(params.n_jk_regions, temp_hist);

	//Place in boxes
	GalaxyBoxes data_boxes(data->gals, data->n_gals, randoms->gals, randoms->n_gals, params.monopole_max);
	GalaxyBoxes random_boxes(randoms->gals, randoms->n_gals, data->gals, data->n_gals, params.monopole_max);

	std::cout << "\nCalculating Monopole DDs\n";
	auto_pair_count_monopole(&data_boxes, &DD);
	DD.normalise(data, data_boxes.n_gals);

	std::cout << "Calculating Monopole DRs\n";
	cross_pair_count_monopole(&data_boxes, &random_boxes, &DR);
	DR.normalise(data, data_boxes.n_gals, randoms, random_boxes.n_gals);

	std::cout << "Calculating Monopole RRs\n";
	auto_pair_count_monopole(&random_boxes, &RR);
	RR.normalise(randoms, random_boxes.n_gals);

	output_xi_1D(params.monopole_filename, params.monopole_output_type, &DD, &DR, &RR, data, randoms);

}

inline void bin_pair_sigma_pi(const Galaxy *gal_1, const Galaxy *gal_2, JKHists2D *hist, const float min_r_sq, const float max_r_sq) {

	float diff_vec[3];
	for (int mm = 0; mm < 3; ++mm) {
		diff_vec[mm] = gal_1->x[mm] - gal_2->x[mm];
	}
	float s = diff_vec[0] * diff_vec[0] + diff_vec[1] * diff_vec[1] + diff_vec[2] * diff_vec[2];

	if (s < max_r_sq && s > min_r_sq) {

#ifdef _ANGULAR_PI_DEFN
		float pi = diff_vec[0] * (gal_1->x_norm[0] + gal_2->x_norm[0]);
		pi += diff_vec[1] * (gal_1->x_norm[1] + gal_2->x_norm[1]);
		pi += diff_vec[2] * (gal_1->x_norm[2] + gal_2->x_norm[2]);
		pi *= 0.5;

#else
		float sum_vec[3], inverse_avg_r = 0.0;
		for (int mm = 0; mm < 3; ++mm) {
			sum_vec[mm] = gal_1->x[mm] + gal_2->x[mm];
			inverse_avg_r += sum_vec[mm] * sum_vec[mm];
		}
		inverse_avg_r = 1.0 / sqrtf(inverse_avg_r);

		float pi = diff_vec[0] * sum_vec[0] * inverse_avg_r;
		pi += diff_vec[1] * sum_vec[1] * inverse_avg_r;
		pi += diff_vec[2] * sum_vec[2] * inverse_avg_r;
#endif

		hist->add_pair(gal_1->jk_region, gal_2->jk_region, sqrtf(s - pi*pi), fabs(pi), gal_1->weight*gal_2->weight);

	}
}

void auto_pair_count_sigma_pi(const GalaxyBoxes *cat, JKHists2D *hist) {

	float max_r_sq = hist->full_hist.ax1.max*hist->full_hist.ax1.max + hist->full_hist.ax2.max*hist->full_hist.ax2.max;
	float min_r_sq = hist->full_hist.ax1.min*hist->full_hist.ax1.min;
	if (hist->full_hist.ax1.min > hist->full_hist.ax2.min)
		min_r_sq = hist->full_hist.ax2.min*hist->full_hist.ax2.min;

#pragma omp parallel
	{

		JKHists2D priv_hist(hist->n_jk_regions, hist->full_hist);

#pragma omp for schedule(dynamic)
		for (int ii = 0; ii < cat->n_boxes_used; ++ii) {
			for (int adj_id = 0; adj_id < 27; ++adj_id) {

				int jj = adj_id_to_box_index(adj_id, cat->boxes[ii].loc3d, cat->id_to_box_index, cat->n_boxes_3d);
				if (ii == jj) {

					for (int kk = 0; kk < cat->boxes[ii].numgals; ++kk) {
						for (int ll = kk + 1; ll < cat->boxes[jj].numgals; ++ll) {

							bin_pair_sigma_pi(&cat->boxes[ii].gal[kk], &cat->boxes[jj].gal[ll], &priv_hist, min_r_sq, max_r_sq);

						}
					}

				}
				else if (jj > ii) {

					for (int kk = 0; kk < cat->boxes[ii].numgals; ++kk) {
						for (int ll = 0; ll < cat->boxes[jj].numgals; ++ll) {

							bin_pair_sigma_pi(&cat->boxes[ii].gal[kk], &cat->boxes[jj].gal[ll], &priv_hist, min_r_sq, max_r_sq);

						}
					}

				}

			}
		}

		(*hist).reduce(priv_hist);

	}

	hist->invert_jk_hists();
	hist->multiply_hists(2.0);
}

void cross_pair_count_sigma_pi(const GalaxyBoxes *cat_1, const GalaxyBoxes *cat_2, JKHists2D *hist) {

	float max_r_sq = hist->full_hist.ax1.max*hist->full_hist.ax1.max + hist->full_hist.ax2.max*hist->full_hist.ax2.max;
	float min_r_sq = hist->full_hist.ax1.min*hist->full_hist.ax1.min;
	if (hist->full_hist.ax1.min > hist->full_hist.ax2.min)
		min_r_sq = hist->full_hist.ax2.min*hist->full_hist.ax2.min;

#pragma omp parallel 
	{

		JKHists2D priv_hist(hist->n_jk_regions, hist->full_hist);

#pragma omp for schedule(dynamic)
		for (int ii = 0; ii < cat_1->n_boxes_used; ++ii) {
			for (int adj_id = 0; adj_id < 27; ++adj_id) {

				int jj = adj_id_to_box_index(adj_id, cat_1->boxes[ii].loc3d, cat_2->id_to_box_index, cat_2->n_boxes_3d);
				if (jj != -1) {

					for (int kk = 0; kk < cat_1->boxes[ii].numgals; ++kk) {
						for (int ll = 0; ll < cat_2->boxes[jj].numgals; ++ll) {

							bin_pair_sigma_pi(&cat_1->boxes[ii].gal[kk], &cat_2->boxes[jj].gal[ll], &priv_hist, min_r_sq, max_r_sq);

						}
					}

				}

			}
		}

		(*hist).reduce(priv_hist);

	}

	hist->invert_jk_hists();
}

void auto_correlate_sigma_pi(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params) {

	std::cout << "Calculating autocorrelation sigma-pi decomposition\n";

	Hist2D temp_hist(params.sigma_min, params.sigma_max, params.sigma_n_bins, params.sigma_log_base, \
		params.pi_min, params.pi_max, params.pi_n_bins, params.pi_log_base);
	JKHists2D DD(params.n_jk_regions, temp_hist);
	JKHists2D DR(params.n_jk_regions, temp_hist);
	JKHists2D RR(params.n_jk_regions, temp_hist);

	//Place in boxes
	GalaxyBoxes data_boxes(data->gals, data->n_gals, randoms->gals, randoms->n_gals, sqrt(params.sigma_max*params.sigma_max + params.pi_max*params.pi_max));
	GalaxyBoxes random_boxes(randoms->gals, randoms->n_gals, data->gals, data->n_gals, sqrt(params.sigma_max*params.sigma_max + params.pi_max*params.pi_max));

	//Run the correlators
	std::cout << "\nCalculating sigma-pi DDs\n";
	auto_pair_count_sigma_pi(&data_boxes, &DD);
	DD.normalise(data, data_boxes.n_gals);
	std::cout << "Calculating sigma-pi DRs\n";
	cross_pair_count_sigma_pi(&data_boxes, &random_boxes, &DR);
	DR.normalise(data, data_boxes.n_gals, randoms, random_boxes.n_gals);
	std::cout << "Calculating sigma-pi RRs\n";
	auto_pair_count_sigma_pi(&random_boxes, &RR);
	RR.normalise(randoms, random_boxes.n_gals);

	output_xi_2D(params.sigma_pi_filename, params.sigma_pi_output_type, "sigma", "pi", &DD, &DR, &RR, data, randoms);

}

inline void bin_pair_s_mu(const Galaxy *gal_1, const Galaxy *gal_2, JKHists2D *hist, const float min_r_sq, const float max_r_sq) {

	float diff_vec[3];
	for (int mm = 0; mm < 3; ++mm) {
		diff_vec[mm] = gal_1->x[mm] - gal_2->x[mm];
	}
	float s = diff_vec[0] * diff_vec[0] + diff_vec[1] * diff_vec[1] + diff_vec[2] * diff_vec[2];

	if (s < max_r_sq && s > min_r_sq) {

#ifdef _ANGULAR_PI_DEFN
		float pi = diff_vec[0] * (gal_1->x_norm[0] + gal_2->x_norm[0]);
		pi += diff_vec[1] * (gal_1->x_norm[1] + gal_2->x_norm[1]);
		pi += diff_vec[2] * (gal_1->x_norm[2] + gal_2->x_norm[2]);
		pi *= 0.5;

#else
		float sum_vec[3], inverse_avg_r = 0.0;
		for (int mm = 0; mm < 3; ++mm) {
			sum_vec[mm] = gal_1->x[mm] + gal_2->x[mm];
			inverse_avg_r += sum_vec[mm] * sum_vec[mm];
		}
		inverse_avg_r = 1.0 / sqrtf(inverse_avg_r);

		float pi = diff_vec[0] * sum_vec[0] * inverse_avg_r;
		pi += diff_vec[1] * sum_vec[1] * inverse_avg_r;
		pi += diff_vec[2] * sum_vec[2] * inverse_avg_r;
#endif

		s = sqrtf(s);

		hist->add_pair(gal_1->jk_region, gal_2->jk_region, s, fabs(pi) / s, gal_1->weight*gal_2->weight);

	}
}

void auto_pair_count_s_mu(const GalaxyBoxes *cat, JKHists2D *hist) {

	float max_r_sq = hist->full_hist.ax1.max * hist->full_hist.ax1.max;
	float min_r_sq = hist->full_hist.ax1.min * hist->full_hist.ax1.min;

#pragma omp parallel
	{

		JKHists2D priv_hist(hist->n_jk_regions, hist->full_hist);

#pragma omp for schedule(dynamic)
		for (int ii = 0; ii < cat->n_boxes_used; ++ii) {
			for (int adj_id = 0; adj_id < 27; ++adj_id) {

				int jj = adj_id_to_box_index(adj_id, cat->boxes[ii].loc3d, cat->id_to_box_index, cat->n_boxes_3d);
				if (ii == jj) {

					for (int kk = 0; kk < cat->boxes[ii].numgals; ++kk) {
						for (int ll = kk + 1; ll < cat->boxes[jj].numgals; ++ll) {

							bin_pair_s_mu(&cat->boxes[ii].gal[kk], &cat->boxes[jj].gal[ll], &priv_hist, min_r_sq, max_r_sq);

						}
					}

				}
				else if (jj > ii) {

					for (int kk = 0; kk < cat->boxes[ii].numgals; ++kk) {
						for (int ll = 0; ll < cat->boxes[jj].numgals; ++ll) {

							bin_pair_s_mu(&cat->boxes[ii].gal[kk], &cat->boxes[jj].gal[ll], &priv_hist, min_r_sq, max_r_sq);

						}
					}

				}

			}
		}


		(*hist).reduce(priv_hist);

	}

	hist->invert_jk_hists();
	hist->multiply_hists(2.0);
}

void cross_pair_count_s_mu(const GalaxyBoxes *cat_1, const GalaxyBoxes *cat_2, JKHists2D *hist) {

	float max_r_sq = hist->full_hist.ax1.max * hist->full_hist.ax1.max;
	float min_r_sq = hist->full_hist.ax1.min * hist->full_hist.ax1.min;

#pragma omp parallel 
	{

		JKHists2D priv_hist(hist->n_jk_regions, hist->full_hist);

#pragma omp for schedule(dynamic)
		for (int ii = 0; ii < cat_1->n_boxes_used; ++ii) {
			for (int adj_id = 0; adj_id < 27; ++adj_id) {

				int jj = adj_id_to_box_index(adj_id, cat_1->boxes[ii].loc3d, cat_2->id_to_box_index, cat_2->n_boxes_3d);
				if (jj != -1) {

					for (int kk = 0; kk < cat_1->boxes[ii].numgals; ++kk) {
						for (int ll = 0; ll < cat_2->boxes[jj].numgals; ++ll) {

							bin_pair_s_mu(&cat_1->boxes[ii].gal[kk], &cat_2->boxes[jj].gal[ll], &priv_hist, min_r_sq, max_r_sq);

						}
					}

				}

			}
		}

		(*hist).reduce(priv_hist);

	}

	hist->invert_jk_hists();
}

void auto_correlate_s_mu(GalaxyCatalogue *data, GalaxyCatalogue *randoms, Params params) {

	std::cout << "Calculating autocorrelation s-mu decomposition\n";

	Hist2D temp_hist(params.s_min, params.s_max, params.s_n_bins, params.s_log_base, \
		0.0, 1.0, params.mu_n_bins, 0);
	JKHists2D DD(params.n_jk_regions, temp_hist);
	JKHists2D DR(params.n_jk_regions, temp_hist);
	JKHists2D RR(params.n_jk_regions, temp_hist);

	//Place in boxes
	GalaxyBoxes data_boxes(data->gals, data->n_gals, randoms->gals, randoms->n_gals, params.s_max);
	GalaxyBoxes random_boxes(randoms->gals, randoms->n_gals, data->gals, data->n_gals, params.s_max);

	//Run the correlators
	std::cout << "\nCalculating s-mu DDs\n";
	auto_pair_count_s_mu(&data_boxes, &DD);
	DD.normalise(data, data_boxes.n_gals);
	std::cout << "Calculating s-mu DRs\n";
	cross_pair_count_s_mu(&data_boxes, &random_boxes, &DR);
	DR.normalise(data, data_boxes.n_gals, randoms, random_boxes.n_gals);
	std::cout << "Calculating s-mu RRs\n";
	auto_pair_count_s_mu(&random_boxes, &RR);
	RR.normalise(randoms, random_boxes.n_gals);

	output_xi_2D(params.s_mu_filename, params.s_mu_output_type, "s", "mu", &DD, &DR, &RR, data, randoms);

}