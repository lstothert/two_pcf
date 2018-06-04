#pragma once
#include <string>
#include "galaxy.h"
#include "hist.h"

#ifdef _HAVE_HDF5
#include "hdf5.h"
#endif

class Params {

public:

	std::string param_filename;

	std::string data_filename;
	std::string data_file_type;
	std::string random_filename;
	std::string random_file_type;

	std::string coord_system;

	std::string ra_x_dataset_name;
	std::string dec_y_dataset_name;
	std::string z_z_dataset_name;
	std::string weight_dataset_name;
	std::string jk_dataset_name;

	bool use_weights;
	int n_threads;

	int n_jk_regions;

	double omega_m;
	double h;
	double z_min;
	double z_max;

	bool plot_monopole;
	std::string monopole_filename;
	std::string monopole_output_type;
	double monopole_log_base;
	double monopole_min;
	double monopole_max;
	int monopole_n_bins;

	bool plot_sigma_pi;
	std::string sigma_pi_filename;
	std::string sigma_pi_output_type;
	double sigma_log_base;
	double sigma_min;
	double sigma_max;
	int sigma_n_bins;
	double pi_log_base;
	double pi_min;
	double pi_max;
	int pi_n_bins;

	bool plot_s_mu;
	std::string s_mu_filename;
	std::string s_mu_output_type;
	double s_log_base;
	double s_min;
	double s_max;
	int s_n_bins;
	int mu_n_bins;

	bool calculate_angular_dd;
	bool calculate_angular_dd_invpweights;
	bool calculate_angular_dr;
	bool calculate_angular_dr_invpweights;
	std::string angular_dd_filename;
	std::string angular_dd_invpweights_filename;
	std::string angular_dr_filename;
	std::string angular_dr_invpweights_filename;
	double theta_max;
	int theta_n_bins;
	double theta_log_base;
	int healpix_order;
	std::string bitwise_weight_dataset_name;
	std::string dither_weight_dataset_name;
	int n_bitwise_runs;

	Params(std::string input_filename);

};

GalaxyCatalogue * read_data_ascii(std::string filename, Cosmology *cosmo, Params params);
#ifdef _HAVE_HDF5
void write_hdf5_double_attribute(double *attr, std::string attr_name, hid_t loc_id);
void write_hdf5_float_attribute(float *attr, std::string attr_name, hid_t loc_id);
void write_hdf5_int_attribute(int *attr, std::string attr_name, hid_t loc_id);
void write_hdf5_string_attribute(std::string attr, std::string attr_name, hid_t loc_id);
void write_hdf5_dataset_1d(double *attr, hsize_t n_data, std::string attr_name, hid_t loc_id);
GalaxyCatalogue * read_data_hdf5(std::string filename, Cosmology *cosmo, Params params);
#endif
#ifdef _HAVE_CFITSIO
GalaxyCatalogue * read_data_fits(std::string filename, Cosmology *cosmo, Params params);
#endif
GalaxyCatalogue * read_file(std::string filename, std::string file_type, Cosmology *cosmo, Params params);

double ls_estimator(double DD, double DR, double RR);
void output_single_xi_1D_ascii(std::string filename, Axis *ax, double *DD, double *DR, double *RR);
void output_xi_1D(std::string filename, std::string file_type, JKHists *DD, JKHists *DR, JKHists *RR, GalaxyCatalogue *data_cat, GalaxyCatalogue *random_cat);
void output_single_xi_2D_ascii(std::string filename, Axis *ax_1, Axis *ax_2, double *DD, double *DR, double *RR);
void output_xi_2D(std::string filename, std::string file_type, std::string axis_1_name, std::string axis_2_name, \
	JKHists2D *DD, JKHists2D *DR, JKHists2D *RR, GalaxyCatalogue *data_cat, GalaxyCatalogue *random_cat);
