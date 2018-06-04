#include "hist.h"
#include "galaxy.h"
#include "io.h"
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _HAVE_HDF5 
#include "hdf5.h"
#endif

Axis::Axis(float input_min, float input_max, int input_n_bins, float log_base) {

	min = input_min;
	minus_min = -input_min;
	max = input_max;
	n_bins = input_n_bins;
	step = (max - min) / (float)n_bins;
	inv_step = 1.0 / step;
	base = log_base;

	if (base > 1.02) {
		inv_log_of_base = 1.0 / log(log_base);
		inv_log_step = (pow(log_base, input_n_bins) - 1) / (input_max - input_min);
		log_binning = true;
	}
	else {
		inv_log_of_base = 0.0;
		inv_log_step = 0.0;
		log_binning = false;
	}

};

double Axis::bin_centre(const int bin_num) {
	return (value_from_bin(bin_num + 1.0) + value_from_bin(bin_num))/2.0;
}

double Axis::bin_width(const int bin_num) {
	return value_from_bin(bin_num + 1.0) - value_from_bin(bin_num);
}

float Axis::bin_num(const float value) {
	if (log_binning) {
		return inv_log_of_base * log1pf(inv_log_step * (value - min));
	}
	else {
		return (value + minus_min)*inv_step;
	}
}

double Axis::value_from_bin(const double bin_num) {
	if (log_binning) {
		return (exp(bin_num / inv_log_of_base) - 1.0) / inv_log_step + min;
	}
	else {
		return min + bin_num*step;
	}
}

Hist::Hist(float input_min, float input_max, int input_n_bins, float log_base) : ax(input_min, input_max, input_n_bins, log_base) {

	bins = (double*)calloc(ax.n_bins, sizeof(double));

};

Hist::Hist(const Hist &copy_obj) : ax(copy_obj.ax.min, copy_obj.ax.max, copy_obj.ax.n_bins, copy_obj.ax.base) {

	bins = (double*)calloc(ax.n_bins, sizeof(double));

	for (int ii = 0; ii < ax.n_bins; ++ii) {
		bins[ii] = copy_obj.bins[ii];
	}

};

Hist::~Hist() {
	free(bins);
}

void Hist::add(const int bin_num, const float weight) {

	if (bin_num < ax.n_bins && bin_num >= 0) {
		bins[bin_num] += weight;
	}

}

void Hist::add(const float value, const float weight) {

	int bin_num = (int)ax.bin_num(value);

	if (bin_num < ax.n_bins && bin_num >= 0 ) {
		bins[bin_num] += weight;
	}

}

void Hist::reduce(const Hist &hist) {

	for (int ii = 0; ii < hist.ax.n_bins; ++ii) {
#pragma omp atomic
		bins[ii] += hist.bins[ii];
	}
}

JKHists::JKHists(int n_jk_regions, const Hist &template_hist) : full_hist(template_hist), n_jk_regions(n_jk_regions) {
	sub_hist = (double*)calloc(n_jk_regions * full_hist.ax.n_bins, sizeof(double));
}

#ifdef _HAVE_HDF5 
void read_hdf5_hist(std::string filename, JKHists **hist) {

	std::cout << "\nReading file: " << filename << "\n";

	hid_t file_id, dataset_id, group_id, attr;
	herr_t status;
	file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	if (file_id < 0) {
		std::cout << "File cannot be opened\n";
		return;
	}

	int n_jk_regions, n_bins;
	double min, max, base;

	group_id = H5Gopen(file_id, "full_result", H5P_DEFAULT);

	attr = H5Aopen(group_id, "n_jk_regions", H5P_DEFAULT);
	status = H5Aread(attr, H5T_NATIVE_INT, &n_jk_regions);
	status = H5Aclose(attr);

	attr = H5Aopen(group_id, "n_bins", H5P_DEFAULT);
	status = H5Aread(attr, H5T_NATIVE_INT, &n_bins);
	status = H5Aclose(attr);

	attr = H5Aopen(group_id, "min", H5P_DEFAULT);
	status = H5Aread(attr, H5T_NATIVE_DOUBLE, &min);
	status = H5Aclose(attr);

	attr = H5Aopen(group_id, "max", H5P_DEFAULT);
	status = H5Aread(attr, H5T_NATIVE_DOUBLE, &max);
	status = H5Aclose(attr);

	attr = H5Aopen(group_id, "base", H5P_DEFAULT);
	status = H5Aread(attr, H5T_NATIVE_DOUBLE, &base);
	status = H5Aclose(attr);

	Hist temp_hist(min, max, n_bins, base);
	(*hist) = new JKHists(n_jk_regions, temp_hist);

	dataset_id = H5Dopen2(group_id, "hist", H5P_DEFAULT);
	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (*hist)->full_hist.bins);
	H5Dclose(dataset_id);

	H5Gclose(group_id);

	// Output jack-knife results if needed
	for (int jk_reg = 0; jk_reg < n_jk_regions; ++jk_reg) {

		std::string group_name = "jk_reg";
		group_name.append(std::to_string((long long)jk_reg));

		group_id = H5Gopen(file_id, group_name.c_str(), H5P_DEFAULT);
		dataset_id = H5Dopen2(group_id, "hist", H5P_DEFAULT);
		H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &((*hist)->sub_hist[jk_reg*(*hist)->full_hist.ax.n_bins]));
		H5Dclose(dataset_id);
		H5Gclose(group_id);

	}

	H5Fclose(file_id);

}

void JKHists::output_hdf5(std::string filename) {

	hid_t file_id, group_id;
	file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	group_id = H5Gcreate2(file_id, "full_result", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	write_hdf5_int_attribute(&(n_jk_regions), "n_jk_regions", group_id);
	write_hdf5_int_attribute(&(full_hist.ax.n_bins), "n_bins", group_id);
	write_hdf5_float_attribute(&(full_hist.ax.min), "min", group_id);
	write_hdf5_float_attribute(&(full_hist.ax.max), "max", group_id);
	write_hdf5_float_attribute(&(full_hist.ax.base), "base", group_id);

	write_hdf5_dataset_1d(full_hist.bins, full_hist.ax.n_bins, "hist", group_id);

	double *axis_data = (double*)malloc(full_hist.ax.n_bins * sizeof(double));
	for (int ii = 0; ii < full_hist.ax.n_bins; ++ii)
		axis_data[ii] = full_hist.ax.bin_centre(ii);
	write_hdf5_dataset_1d(axis_data, full_hist.ax.n_bins, "Bin_centre", group_id);
	for (int ii = 0; ii < full_hist.ax.n_bins; ++ii)
		axis_data[ii] = full_hist.ax.bin_width(ii);
	write_hdf5_dataset_1d(axis_data, full_hist.ax.n_bins, "Bin_width", group_id);
	free(axis_data);

	H5Gclose(group_id);

	// Output jack-knife results if needed
	for (int jk_reg = 0; jk_reg < n_jk_regions; ++jk_reg) {

		std::string group_name = "jk_reg";
		group_name.append(std::to_string((long long)jk_reg));

		group_id = H5Gcreate2(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		write_hdf5_dataset_1d(&(sub_hist[jk_reg*full_hist.ax.n_bins]), full_hist.ax.n_bins, "hist", group_id);
		H5Gclose(group_id);

	}

	H5Fclose(file_id);
}

#endif

void JKHists::divide_by_hist(const JKHists &input_hist) {

	for (int ii = 0; ii < full_hist.ax.n_bins; ++ii) 
		full_hist.bins[ii] /= input_hist.full_hist.bins[ii];
	for (int ii = 0; ii < n_jk_regions * full_hist.ax.n_bins; ++ii)
		sub_hist[ii] /= input_hist.sub_hist[ii];

}

JKHists::~JKHists() {
	free(sub_hist);
}

void JKHists::reduce(const JKHists &input_hist) {

	full_hist.reduce(input_hist.full_hist);
	for (int ii = 0; ii < n_jk_regions*full_hist.ax.n_bins; ++ii) {
#pragma omp atomic
		sub_hist[ii] += input_hist.sub_hist[ii];
	}
}

void JKHists::add_pair(const int jk_region_1, const int jk_region_2, const float value, const float weight) {

	int bin_num = full_hist.ax.bin_num(value);

	if (bin_num < full_hist.ax.n_bins && bin_num >= 0) {

		full_hist.bins[bin_num] += weight;
		sub_hist[jk_region_1*full_hist.ax.n_bins + bin_num] += weight;

		if (jk_region_1 != jk_region_2) {
			sub_hist[jk_region_2*full_hist.ax.n_bins + bin_num] += weight;
		}

	}
}

double JKHists::get(const int bin_num) {
	return full_hist.bins[bin_num];
}

double JKHists::get(const int jk_region, const int bin_num) {
	return sub_hist[jk_region*full_hist.ax.n_bins + bin_num];
}

void JKHists::invert_jk_hists() {

	for (int ii = 0; ii < n_jk_regions; ++ii) {
		for (int jj = 0; jj < full_hist.ax.n_bins; ++jj) {
			sub_hist[ii*full_hist.ax.n_bins + jj] = full_hist.bins[jj] - sub_hist[ii*full_hist.ax.n_bins + jj];
		}
	}
}

void JKHists::multiply_hists(double factor) {

	for (int ii = 0; ii < full_hist.ax.n_bins; ++ii) {
		full_hist.bins[ii] *= factor;
	}

	for (int ii = 0; ii < n_jk_regions*full_hist.ax.n_bins; ++ii) {
		sub_hist[ii] *= factor;
	}
}

void JKHists::normalise(const GalaxyCatalogue *cat, int n_gals) {

	// Calculate sum of weights total and in all jk regions
	double sum_weights = 0.0;
	double *sum_weights_jk = (double*)calloc(n_jk_regions, sizeof(double));
	int *n_gals_jk = (int*)calloc(n_jk_regions, sizeof(int));
	for (int ii = 0; ii < n_gals; ++ii) {
		sum_weights += cat->gals[ii].weight;
		sum_weights_jk[cat->gals[ii].jk_region] += cat->gals[ii].weight;
		++n_gals_jk[cat->gals[ii].jk_region];
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk[ii] = sum_weights - sum_weights_jk[ii];
	}

	// Calculate (N-1)N equivilant with weights
	sum_weights = sum_weights * sum_weights * (1.0 - 1.0 / (double)n_gals);
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk[ii] = sum_weights_jk[ii] * sum_weights_jk[ii] * (1.0 - 1.0 / n_gals_jk[ii]);
	}

	// Normalise histograms
	for (int ii = 0; ii < full_hist.ax.n_bins; ++ii) {
		full_hist.bins[ii] /= sum_weights;
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		for (int jj = 0; jj < full_hist.ax.n_bins; ++jj) {
			sub_hist[ii*full_hist.ax.n_bins + jj] /= sum_weights_jk[ii];
		}
	}

}

void JKHists::normalise(const GalaxyCatalogue *cat_1, int n_gals_1, const GalaxyCatalogue *cat_2, int n_gals_2) {

	// Calculate sum of weights total and in all jk regions
	double sum_weights_1 = 0.0;
	double *sum_weights_jk_1 = (double*)calloc(n_jk_regions, sizeof(double));
	int *n_gals_jk_1 = (int*)calloc(n_jk_regions, sizeof(int));
	for (int ii = 0; ii < n_gals_1; ++ii) {
		sum_weights_1 += cat_1->gals[ii].weight;
		sum_weights_jk_1[cat_1->gals[ii].jk_region] += cat_1->gals[ii].weight;
		++n_gals_jk_1[cat_1->gals[ii].jk_region];
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk_1[ii] = sum_weights_1 - sum_weights_jk_1[ii];
	}

	// Calculate sum of weights total and in all jk regions
	double sum_weights_2 = 0.0;
	double *sum_weights_jk_2 = (double*)calloc(n_jk_regions, sizeof(double));
	int *n_gals_jk_2 = (int*)calloc(n_jk_regions, sizeof(int));
	for (int ii = 0; ii < n_gals_2; ++ii) {
		sum_weights_2 += cat_2->gals[ii].weight;
		sum_weights_jk_2[cat_2->gals[ii].jk_region] += cat_2->gals[ii].weight;
		++n_gals_jk_2[cat_2->gals[ii].jk_region];
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk_2[ii] = sum_weights_2 - sum_weights_jk_2[ii];
	}

	// Calculate (N-1)N equivilant with weights
	sum_weights_1 *= sum_weights_2;
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk_1[ii] *= sum_weights_jk_2[ii];
	}

	// Normalise histograms
	for (int ii = 0; ii < full_hist.ax.n_bins; ++ii) {
		full_hist.bins[ii] /= sum_weights_1;
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		for (int jj = 0; jj < full_hist.ax.n_bins; ++jj) {
			sub_hist[ii*full_hist.ax.n_bins + jj] /= sum_weights_jk_1[ii];
		}
	}

}

double JKHists::interp(const double value) {

        if (value < full_hist.ax.min){
            return full_hist.bins[0];
        }
        else if(value > full_hist.ax.max) {
	    return full_hist.bins[full_hist.ax.n_bins - 1];
	}
	else {
		float bin_num = full_hist.ax.bin_num(value) - 0.5;
		if (bin_num < 0) {
			return full_hist.bins[0];
		}
		else if (bin_num > full_hist.ax.n_bins - 1){
			return full_hist.bins[full_hist.ax.n_bins - 1];
		}
		else {
			return full_hist.bins[(int)bin_num] + (full_hist.bins[(int)bin_num + 1] - full_hist.bins[(int)bin_num]) * \
				(value - full_hist.ax.bin_centre((int)bin_num)) / (full_hist.ax.bin_centre((int)bin_num + 1) - full_hist.ax.bin_centre((int)bin_num));
		}
	}

}

Hist2D::Hist2D(float input_min_1, float input_max_1, int input_n_bins_1, float axis_1_log_base,
	float input_min_2, float input_max_2, int input_n_bins_2, float axis_2_log_base) :
	ax1(input_min_1, input_max_1, input_n_bins_1, axis_1_log_base) , 
	ax2(input_min_2, input_max_2, input_n_bins_2, axis_2_log_base) {

	total_bins = ax1.n_bins*ax2.n_bins;
	bins = (double*)calloc(total_bins, sizeof(double));

}

Hist2D::Hist2D(const Hist2D &copy_obj) : 
	ax1(copy_obj.ax1.min, copy_obj.ax1.max, copy_obj.ax1.n_bins, copy_obj.ax1.base), 
	ax2(copy_obj.ax2.min, copy_obj.ax2.max, copy_obj.ax2.n_bins, copy_obj.ax2.base) {

	total_bins = ax1.n_bins*ax2.n_bins;
	bins = (double*)malloc(total_bins * sizeof(double));

	for (int ii = 0; ii < total_bins; ++ii) {
		bins[ii] = copy_obj.bins[ii];
	}

}

Hist2D::~Hist2D() {
	free(bins);
}

void Hist2D::add(const int bin_num_1, const int bin_num_2, const float weight) {

	if (bin_num_1 < ax1.n_bins && bin_num_1 >= 0 && bin_num_2 < ax2.n_bins && bin_num_2 >= 0) {
		bins[bin_num_1*ax2.n_bins + bin_num_2] += weight;
	}

}

void Hist2D::add(const float value_1, const float value_2, const float weight) {

	int bin_num[2];
	bin_num[0]= (int)ax1.bin_num(value_1);
	bin_num[1] = (int)ax2.bin_num(value_2);
	if (bin_num[0] >= 0 && bin_num[0] < ax1.n_bins && bin_num[1] >= 0 && bin_num[1] < ax2.n_bins) {
		bins[ax2.n_bins*bin_num[0] + bin_num[1]] += weight;
	}
}

void Hist2D::reduce(const Hist2D &hist) {

	for (int ii = 0; ii < hist.ax1.n_bins*hist.ax2.n_bins; ++ii) {
#pragma omp atomic
		bins[ii] += hist.bins[ii];
	}
}

JKHists2D::JKHists2D(int n_jk_regions, const Hist2D &template_hist) : full_hist(template_hist), n_jk_regions(n_jk_regions) {
	sub_hist = (double*)calloc(n_jk_regions * full_hist.total_bins, sizeof(double));
}

JKHists2D::~JKHists2D() {
	free(sub_hist);
}

void JKHists2D::reduce(const JKHists2D &input_hist) {

	full_hist.reduce(input_hist.full_hist);
	for (int ii = 0; ii < n_jk_regions*full_hist.total_bins; ++ii) {
#pragma omp atomic
		sub_hist[ii] += input_hist.sub_hist[ii];
	}

}

void JKHists2D::invert_jk_hists(){

	for (int ii = 0; ii < n_jk_regions; ++ii) {
		for (int jj = 0; jj < full_hist.total_bins; ++jj) {
			sub_hist[ii*full_hist.total_bins + jj] = full_hist.bins[jj] - sub_hist[ii*full_hist.total_bins + jj];
		}
	}
}

void JKHists2D::multiply_hists(double factor) {

	for (int ii = 0; ii < full_hist.total_bins; ++ii) {
		full_hist.bins[ii] *= factor;
	}

	for (int ii = 0; ii < n_jk_regions*full_hist.total_bins; ++ii) {
		sub_hist[ii] *= factor;
	}
}

void JKHists2D::add_pair(const int jk_region_1, const int jk_region_2, const float value_1, const float value_2, const float weight) {

	int bin_num_1 = full_hist.ax1.bin_num(value_1);
	if (bin_num_1 < full_hist.ax1.n_bins && bin_num_1 >= 0) {

		int bin_num_2 = full_hist.ax2.bin_num(value_2);
		if (bin_num_2 < full_hist.ax2.n_bins && bin_num_2 >= 0) {

			int bin_num = bin_num_1 * full_hist.ax2.n_bins + bin_num_2;
			full_hist.bins[bin_num] += weight;
			sub_hist[jk_region_1*full_hist.total_bins + bin_num] += weight;

			if (jk_region_1 != jk_region_2) {
				sub_hist[jk_region_2*full_hist.total_bins + bin_num] += weight;
			}

		}
	}
}

double JKHists2D::get(const int bin_num_1, const int bin_num_2) {
	return full_hist.bins[bin_num_1 * full_hist.ax2.n_bins + bin_num_2];
}

double JKHists2D::get(const int jk_region, const int bin_num_1, const int bin_num_2) {
	return sub_hist[jk_region*full_hist.total_bins + bin_num_1 * full_hist.ax2.n_bins + bin_num_2];
}

void JKHists2D::normalise(const GalaxyCatalogue *cat, int n_gals) {

	// Calculate sum of weights total and in all jk regions
	double sum_weights = 0.0;
	double *sum_weights_jk = (double*)calloc(n_jk_regions, sizeof(double));
	int *n_gals_jk = (int*)calloc(n_jk_regions, sizeof(int));
	for (int ii = 0; ii < n_gals; ++ii) {
		sum_weights += cat->gals[ii].weight;
		sum_weights_jk[cat->gals[ii].jk_region] += cat->gals[ii].weight;
		++n_gals_jk[cat->gals[ii].jk_region];
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk[ii] = sum_weights - sum_weights_jk[ii];
	}

	// Calculate (N-1)N equivilant with weights
	sum_weights = sum_weights * sum_weights * (1.0 - 1.0 / (double)n_gals);
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk[ii] = sum_weights_jk[ii] * sum_weights_jk[ii] * (1.0 - 1.0 / n_gals_jk[ii]);
	}

	// Normalise histograms
	for (int ii = 0; ii < full_hist.total_bins; ++ii) {
		full_hist.bins[ii] /= sum_weights;
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		for (int jj = 0; jj < full_hist.total_bins; ++jj) {
			sub_hist[ii*full_hist.total_bins + jj] /= sum_weights_jk[ii];
		}
	}

}

void JKHists2D::normalise(const GalaxyCatalogue *cat_1, int n_gals_1, const GalaxyCatalogue *cat_2, int n_gals_2) {

	// Calculate sum of weights total and in all jk regions
	double sum_weights_1 = 0.0;
	double *sum_weights_jk_1 = (double*)calloc(n_jk_regions, sizeof(double));
	int *n_gals_jk_1 = (int*)calloc(n_jk_regions, sizeof(int));
	for (int ii = 0; ii < n_gals_1; ++ii) {
		sum_weights_1 += cat_1->gals[ii].weight;
		sum_weights_jk_1[cat_1->gals[ii].jk_region] += cat_1->gals[ii].weight;
		++n_gals_jk_1[cat_1->gals[ii].jk_region];
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk_1[ii] = sum_weights_1 - sum_weights_jk_1[ii];
	}

	// Calculate sum of weights total and in all jk regions
	double sum_weights_2 = 0.0;
	double *sum_weights_jk_2 = (double*)calloc(n_jk_regions, sizeof(double));
	int *n_gals_jk_2 = (int*)calloc(n_jk_regions, sizeof(int));
	for (int ii = 0; ii < n_gals_2; ++ii) {
		sum_weights_2 += cat_2->gals[ii].weight;
		sum_weights_jk_2[cat_2->gals[ii].jk_region] += cat_2->gals[ii].weight;
		++n_gals_jk_2[cat_2->gals[ii].jk_region];
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk_2[ii] = sum_weights_2 - sum_weights_jk_2[ii];
	}

	// Calculate (N-1)N equivilant with weights
	sum_weights_1 *= sum_weights_2;
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		sum_weights_jk_1[ii] *= sum_weights_jk_2[ii];
	}

	// Normalise histograms
	for (int ii = 0; ii < full_hist.total_bins; ++ii) {
		full_hist.bins[ii] /= sum_weights_1;
	}
	for (int ii = 0; ii < n_jk_regions; ++ii) {
		for (int jj = 0; jj < full_hist.total_bins; ++jj) {
			sub_hist[ii*full_hist.total_bins + jj] /= sum_weights_jk_1[ii];
		}
	}

}
