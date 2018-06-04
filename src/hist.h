#pragma once
#include <string>
#include "galaxy.h"

class Axis {

public:
	int n_bins;
	bool log_binning;
	float min, minus_min, max, step, inv_step, base, inv_log_of_base, inv_log_step;
	Axis(float input_min, float input_max, int input_n_bins, float log_base);
	~Axis() {};

	double bin_centre(const int bin_num);
	double bin_width(const int bin_num);
	double value_from_bin(const double bin_num);
	float bin_num(const float value);

};

class Hist {

public:

	Axis ax;
	double *bins;

	Hist(float input_min, float input_max, int input_n_bins, float log_base);
	Hist(const Hist &copy_obj);
	~Hist();

	void add(const int bin_num, const float weight = 1.0);
	void add(const float value, const float weight = 1.0);
	void reduce(const Hist &hist);

};

class JKHists {

public:

	int n_jk_regions;
	Hist full_hist;
	double *sub_hist;

	JKHists(int jk_regions, const Hist &template_hist);
	~JKHists();

	void add_pair(const int jk_region_1, const int jk_region_2, const float value, const float weight = 1.0);
	void reduce(const JKHists &input_hist);
	void divide_by_hist(const JKHists &input_hist);
	void invert_jk_hists();
	void multiply_hists(double factor);
	double get(const int bin_num);
	double get(const int jk_region, const int bin_num);
	void normalise(const GalaxyCatalogue *cat, int n_gals);
	void normalise(const GalaxyCatalogue *cat_1, int n_gals_1, const GalaxyCatalogue *cat_2, int n_gals_2);
	double interp(const double value);

#ifdef _HAVE_HDF5 
	void output_hdf5(std::string filename);
#endif

};

#ifdef _HAVE_HDF5 
void read_hdf5_hist(std::string filename, JKHists **hist);
#endif

class Hist2D {

public:
	Axis ax1, ax2;
	double *bins;
	int total_bins;

	Hist2D(float input_min_1, float input_max_1, int input_n_bins_1, float axis_1_log_base, 
		float input_min_2, float input_max_2, int input_n_bins_2, float axis_2_log_base);
	Hist2D(const Hist2D &copy_obj);
	~Hist2D();

	void add(const int bin_num_1, const int bin_num_2, const float weight = 1.0);
	void add(const float value_1, const float value_2, const float weight = 1.0);
	void reduce(const Hist2D &hist);

};

class JKHists2D {

public:

	int n_jk_regions;
	Hist2D full_hist;
	double *sub_hist;

	JKHists2D(int jk_regions, const Hist2D &template_hist);
	~JKHists2D();
	void add_pair(const int jk_region_1, const int jk_region_2, const float value_1, const float value_2, const float weight = 1.0);
	void reduce(const JKHists2D &input_hist);
	void invert_jk_hists();
	void multiply_hists(double factor);
	double get(const int bin_num_1, const int bin_num_2);
	double get(const int jk_region, const int bin_num_1, const int bin_num_2);
	void normalise(const GalaxyCatalogue *cat, int n_gals);
	void normalise(const GalaxyCatalogue *cat_1, int n_gals_1, const GalaxyCatalogue *cat_2, int n_gals_2);

};