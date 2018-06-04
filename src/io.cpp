#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include "io.h"
#include "galaxy.h"
#include "cosmology.h"

#ifdef _HAVE_HDF5 
#include "hdf5.h"
#endif

// Params constructor, read param file and save params, there must be a cleaner way to do this! (Dictionary?)
Params::Params(std::string input_filename) {

	std::cout << "Reading parameter file: " << input_filename << "\n";

	std::ifstream file(input_filename);
	std::string line, key, value;

	n_jk_regions = 1;

	while (std::getline(file, line)) {

		if (line[0] != '#' && line.find("=") != std::string::npos) {

			std::stringstream not_trimmed(line);
			std::getline(not_trimmed, line, '#');

			std::stringstream ss(line);
			std::getline(ss, key, '=');
			std::getline(ss, value);
			key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
			value.erase(std::remove(value.begin(), value.end(), ' '), value.end());

			std::cout << "Key: " << key << ",\t Value: " << value << ",\t Result: ";

			bool key_found = true;
			if (key == "data_filename") {
				data_filename = value;
			}
			else if (key == "data_file_type") {
				data_file_type = value;
			}
			else if (key == "random_filename") {
				random_filename = value;
			}
			else if (key == "random_file_type") {
				random_file_type = value;
			}
			else if (key == "coord_system") {
				coord_system = value;
			}
			else if (key == "ra_x_dataset_name") {
				ra_x_dataset_name = value;
			}
			else if (key == "dec_y_dataset_name") {
				dec_y_dataset_name = value;
			}
			else if (key == "z_z_dataset_name") {
				z_z_dataset_name = value;
			}
			else if (key == "weight_dataset_name") {
				weight_dataset_name = value;
			}
			else if (key == "jk_dataset_name") {
				jk_dataset_name = value;
			}
			else if (key == "use_weights") {
				use_weights = std::stoi(value);
			}
			else if (key == "n_threads") {
				n_threads = std::stoi(value);
			}
			else if (key == "n_jk_regions") {
				n_jk_regions = std::stoi(value);
			}
			else if (key == "omega_m") {
				omega_m = std::stod(value);
			}
			else if (key == "h") {
				h = std::stod(value);
			}
			else if (key == "z_min") {
				z_min = std::stod(value);
			}
			else if (key == "z_max") {
				z_max = std::stod(value);
			}

			else if (key == "plot_monopole") {
				plot_monopole = std::stoi(value);
			}
			else if (key == "monopole_filename") {
				monopole_filename = value;
			}
			else if (key == "monopole_output_type") {
				monopole_output_type = value;
			}
			else if (key == "monopole_log_base") {
				monopole_log_base = std::stod(value);
			}
			else if (key == "monopole_min") {
				monopole_min = std::stod(value);
			}
			else if (key == "monopole_max") {
				monopole_max = std::stod(value);
			}
			else if (key == "monopole_n_bins") {
				monopole_n_bins = std::stoi(value);
			}

			else if (key == "plot_sigma_pi") {
				plot_sigma_pi = std::stoi(value);
			}
			else if (key == "sigma_pi_filename") {
				sigma_pi_filename = value;
			}
			else if (key == "sigma_pi_output_type") {
				sigma_pi_output_type = value;
			}
			else if (key == "sigma_log_base") {
				sigma_log_base = std::stod(value);
			}
			else if (key == "sigma_min") {
				sigma_min = std::stod(value);
			}
			else if (key == "sigma_max") {
				sigma_max = std::stod(value);
			}
			else if (key == "sigma_n_bins") {
				sigma_n_bins = std::stoi(value);
			}
			else if (key == "pi_log_base") {
				pi_log_base = std::stod(value);
			}
			else if (key == "pi_min") {
				pi_min = std::stod(value);
			}
			else if (key == "pi_max") {
				pi_max = std::stod(value);
			}
			else if (key == "pi_n_bins") {
				pi_n_bins = std::stoi(value);
			}

			else if (key == "plot_s_mu") {
				plot_s_mu = std::stoi(value);
			}
			else if (key == "s_mu_filename") {
				s_mu_filename = value;
			}
			else if (key == "s_mu_output_type") {
				s_mu_output_type = value;
			}
			else if (key == "s_log_base") {
				s_log_base = std::stod(value);
			}
			else if (key == "s_min") {
				s_min = std::stod(value);
			}
			else if (key == "s_max") {
				s_max = std::stod(value);
			}
			else if (key == "s_n_bins") {
				s_n_bins = std::stoi(value);
			}
			else if (key == "mu_n_bins") {
				mu_n_bins = std::stoi(value);
			}

			else if (key == "calculate_angular_dd") {
				calculate_angular_dd = std::stoi(value);
			}
			else if (key == "angular_dd_filename") {
				angular_dd_filename = value;
			}
			else if (key == "calculate_angular_dd_invpweights") {
				calculate_angular_dd_invpweights = std::stoi(value);
			}
			else if (key == "angular_dd_invpweights_filename") {
				angular_dd_invpweights_filename = value;
			}
			else if (key == "calculate_angular_dr") {
				calculate_angular_dr = std::stoi(value);
			}
			else if (key == "angular_dr_filename") {
				angular_dr_filename = value;
			}
			else if (key == "calculate_angular_dr_invpweights") {
				calculate_angular_dr_invpweights = std::stoi(value);
			}
			else if (key == "angular_dr_invpweights_filename") {
				angular_dr_invpweights_filename = value;
			}
			else if (key == "theta_max") {
				theta_max = std::stod(value);
			}
			else if (key == "theta_n_bins") {
				theta_n_bins = std::stoi(value);
			}
			else if (key == "theta_log_base") {
				theta_log_base = std::stod(value);
			}
			else if (key == "healpix_order") {
				healpix_order = std::stoi(value);
			}
			else if (key == "bitwise_weight_dataset_name") {
				bitwise_weight_dataset_name = value;
			}
			else if (key == "dither_weight_dataset_name") {
				dither_weight_dataset_name = value;
			}
			else if (key == "n_bitwise_runs") {
				n_bitwise_runs = std::stoi(value);
			}

			else {
				key_found = false;
				std::cout << "KEY NOT FOUND\n";
			}

			if (key_found)
				std::cout << "KEY FOUND\n";

		}
	}

	file.close();
	std::cout << "\n";

	if (n_jk_regions <= 1)
		n_jk_regions = 1;

}

// Reads ascii data file, comments starting with # ignored
GalaxyCatalogue * read_data_ascii(std::string filename, Cosmology *cosmo, Params params) {

	std::ifstream file(filename);

	if (file.is_open()) {
		std::cout << "File successfully opened\n";
	}
	else {
		std::cout << "File not open\n";
	}

	std::string line, part;

	//Count the number of lines in file
	int num_lines = 0;
	while (std::getline(file, line)) {
		if (line[0] != '#' && line != "")
			++num_lines;
	}

	std::cout << "Number of lines found: " << num_lines << "\n";

	file.clear();
	file.seekg(0);

	GalaxyCatalogue *cat = new GalaxyCatalogue(num_lines, params.n_jk_regions);

	if (params.coord_system == "equatorial") {

		std::cout << "Reading equatorial data\n";

		GalaxyEq *eqinput = new GalaxyEq[num_lines];

		num_lines = 0;
		while (std::getline(file, line)) {

			if (line[0] != '#' && line != "") {

				// Allow different separation characters, spaces will all be trimmed.
				char sepchar = ' ';
				if (line.find('\t') != std::string::npos)
					sepchar = '\t';
				if (line.find(',') != std::string::npos)
					sepchar = ',';

				std::stringstream ss(line);

				do {
					std::getline(ss, part, sepchar);
				} while (part == "");
				part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
				eqinput[num_lines].ra = stod(part);

				do {
					std::getline(ss, part, sepchar);
				} while (part == "");
				part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
				eqinput[num_lines].dec = stod(part);

				do {
					std::getline(ss, part, sepchar);
				} while (part == "");
				part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
				eqinput[num_lines].z = stod(part);

				if (params.use_weights) {
					do {
						std::getline(ss, part, sepchar);
					} while (part == "");
					part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
					cat->gals[num_lines].weight = stof(part);
				}
				else {
					cat->gals[num_lines].weight = 1.0;
				}

				if (params.n_jk_regions > 1) {
					do {
						std::getline(ss, part, sepchar);
					} while (part == "");
					part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
					cat->gals[num_lines].jk_region = stoi(part);
				}
				else {
					cat->gals[num_lines].jk_region = 0;
				}

				cat->gals[num_lines].d_comov = cosmo->comov_dist(eqinput[num_lines].z);
				cat->gals[num_lines].set(eqinput[num_lines].ra, eqinput[num_lines].dec);

				++num_lines;

			}
		}

		delete[] eqinput;
	}
	else if (params.coord_system == "cartesian") {

		std::cout << "Reading cartesian data\n";

		num_lines = 0;
		while (std::getline(file, line)) {
			if (line[0] != '#' && line != "") {

				// Allow different separation characters, spaces will all be trimmed.
				char sepchar = ' ';
				if (line.find('\t') != std::string::npos)
					sepchar = '\t';
				if (line.find(',') != std::string::npos)
					sepchar = ',';

				std::stringstream ss(line);

				do {
					std::getline(ss, part, sepchar);
				} while (part == "");
				part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
				cat->gals[num_lines].x[0] = stod(part);

				do {
					std::getline(ss, part, sepchar);
				} while (part == "");
				part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
				cat->gals[num_lines].x[1] = stod(part);

				do {
					std::getline(ss, part, sepchar);
				} while (part == "");
				part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
				cat->gals[num_lines].x[2] = stod(part);

				if (params.use_weights) {
					do {
						std::getline(ss, part, sepchar);
					} while (part == "");
					part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
					cat->gals[num_lines].weight = stof(part);
				}
				else {
					cat->gals[num_lines].weight = 1.0;
				}
				if (params.n_jk_regions > 1) {
					do {
						std::getline(ss, part, sepchar);
					} while (part == "");
					part.erase(std::remove(part.begin(), part.end(), ' '), part.end());
					cat->gals[num_lines].jk_region = stoi(part);
				}
				else {
					cat->gals[num_lines].jk_region = 0;
				}

				cat->gals[num_lines].d_comov = sqrt(cat->gals[num_lines].x[0] * cat->gals[num_lines].x[0] + \
					cat->gals[num_lines].x[1] * cat->gals[num_lines].x[1] + \
					cat->gals[num_lines].x[2] * cat->gals[num_lines].x[2]);

				cat->gals[num_lines].x_norm[0] = cat->gals[num_lines].x[0] / cat->gals[num_lines].d_comov;
				cat->gals[num_lines].x_norm[1] = cat->gals[num_lines].x[1] / cat->gals[num_lines].d_comov;
				cat->gals[num_lines].x_norm[2] = cat->gals[num_lines].x[2] / cat->gals[num_lines].d_comov;

				++num_lines;

			}
		}

	}
	else {
		std::cout << "Coord system not found\n";
	}

	file.close();

	return(cat);

}

#ifdef _HAVE_HDF5
GalaxyCatalogue * read_data_hdf5(std::string filename, Cosmology *cosmo, Params params) {

	hid_t file_id, dataset_id, dspace;
	herr_t status;
	file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	if (file_id < 0) {
		std::cout << "File cannot be opened\n";
	}
	else {
		std::cout << "File successfully opened\n";
	}

	dataset_id = H5Dopen2(file_id, params.ra_x_dataset_name.c_str(), H5P_DEFAULT);
	dspace = H5Dget_space(dataset_id);
	int num_lines = H5Sget_simple_extent_npoints(dspace);
	H5Dclose(dataset_id);

	std::cout << "N points in file: " << num_lines << "\n";

	GalaxyCatalogue *cat = new GalaxyCatalogue(num_lines, params.n_jk_regions);
	double *col1 = (double*)malloc(num_lines * sizeof(double));
	double *col2 = (double*)malloc(num_lines * sizeof(double));
	double *col3 = (double*)malloc(num_lines * sizeof(double));

	dataset_id = H5Dopen2(file_id, params.ra_x_dataset_name.c_str(), H5P_DEFAULT);
	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, col1);
	H5Dclose(dataset_id);

	dataset_id = H5Dopen2(file_id, params.dec_y_dataset_name.c_str(), H5P_DEFAULT);
	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, col2);
	H5Dclose(dataset_id);

	dataset_id = H5Dopen2(file_id, params.z_z_dataset_name.c_str(), H5P_DEFAULT);
	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, col3);
	H5Dclose(dataset_id);

	if (params.coord_system == "equatorial") {
		std::cout << "Reading equatorial data\n";
		for (int ii = 0; ii < num_lines; ++ii) {
			if (col3[ii] >= 0.0) {
				cat->gals[ii].d_comov = cosmo->comov_dist(col3[ii]);
			}
			else {
				cat->gals[ii].d_comov = -1.0;
			}
			cat->gals[ii].set(col1[ii], col2[ii]);
		}
	}
	else if (params.coord_system == "cartesian") {
		std::cout << "Reading cartesian data\n";

		for (int ii = 0; ii < num_lines; ++ii) {

			cat->gals[ii].x[0] = col1[ii];
			cat->gals[ii].x[1] = col2[ii];
			cat->gals[ii].x[2] = col3[ii];

			cat->gals[ii].d_comov = sqrt(col1[ii] * col1[ii] + col2[ii] * col2[ii] + col3[ii] * col3[ii]);

			cat->gals[ii].x_norm[0] = col1[ii] / cat->gals[ii].d_comov;
			cat->gals[ii].x_norm[1] = col2[ii] / cat->gals[ii].d_comov;
			cat->gals[ii].x_norm[2] = col3[ii] / cat->gals[ii].d_comov;

		}
	}
	else {
		std::cout << "Cannot find coordinate system\n";
	}
	free(col1); free(col2); free(col3);

	if (params.use_weights) {
		std::cout << "Reading input weights\n";
		double *weights = (double*)malloc(num_lines * sizeof(double));
		dataset_id = H5Dopen2(file_id, params.weight_dataset_name.c_str(), H5P_DEFAULT);
		H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, weights);
		H5Dclose(dataset_id);
		for (int ii = 0; ii < num_lines; ++ii) {
			cat->gals[ii].weight = weights[ii];
		}
		free(weights);
	}
	else {
		for (int ii = 0; ii < num_lines; ++ii)
			cat->gals[ii].weight = 1.0;
	}

	if (params.n_jk_regions > 1) {
		std::cout << "Reading jackknife regions\n";
		int *jk_regions = (int*)malloc(num_lines * sizeof(int));
		dataset_id = H5Dopen2(file_id, params.jk_dataset_name.c_str(), H5P_DEFAULT);
		H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, jk_regions);
		H5Dclose(dataset_id);
		for (int ii = 0; ii < num_lines; ++ii) {
			cat->gals[ii].jk_region = jk_regions[ii];
		}
		free(jk_regions);
	}
	else {
		for (int ii = 0; ii < num_lines; ++ii)
			cat->gals[ii].jk_region = 0;
	}

#ifdef _USE_INV_WEIGHTS

	std::cout << "Reading bitwise weights\n";
	long long *bitwise_weights = (long long*)malloc(num_lines * sizeof(long long));
	int *total = (int*)calloc(num_lines, sizeof(int));
	std::string bitwise_weight_dataset_name_ii;
	std::stringstream sstm;
	for (int jj = 0; jj < _N_MASK_INTS; ++jj) {
		sstm << params.bitwise_weight_dataset_name << jj;
		bitwise_weight_dataset_name_ii = sstm.str();
		dataset_id = H5Dopen2(file_id, bitwise_weight_dataset_name_ii.c_str(), H5P_DEFAULT);
		H5Dread(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, bitwise_weights);
		H5Dclose(dataset_id);
		for (int ii = 0; ii < num_lines; ++ii) {
			cat->gals[ii].bitmask[jj] = bitwise_weights[ii];
			total[ii] += __builtin_popcountll(bitwise_weights[ii]);
		}
		sstm.str("");
	}
	free(bitwise_weights);

	for (int ii = 0; ii < num_lines; ++ii) {
		cat->gals[ii].input_weight = cat->gals[ii].weight;
		if (total[ii] == 0) {
			total[ii] = 1;
		}
	}

#ifdef _USE_DITHER_WEIGHTS

	std::cout << "Reading dither weights\n";
	long long *dither_weights = (long long*)malloc(num_lines * sizeof(long long));
	int *dither_total = (int*)calloc(num_lines, sizeof(int));
	std::string dither_weight_dataset_name_ii;
	for (int jj = 0; jj < _N_MASK_INTS; ++jj) {
		sstm << params.dither_weight_dataset_name << jj;
		dither_weight_dataset_name_ii = sstm.str();
		dataset_id = H5Dopen2(file_id, dither_weight_dataset_name_ii.c_str(), H5P_DEFAULT);
		H5Dread(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, dither_weights);
		H5Dclose(dataset_id);
		for (int ii = 0; ii < num_lines; ++ii) {
			cat->gals[ii].dither_mask[jj] = dither_weights[ii];
			dither_total[ii] += __builtin_popcountll(dither_weights[ii]);
		}
		sstm.str("");
	}
	free(dither_weights);

	for (int ii = 0; ii < num_lines; ++ii) {
		if (total[ii] != params.n_bitwise_runs) {
			cat->gals[ii].weight *= (double)dither_total[ii] / (double)total[ii];
		}
	}
	free(dither_total);

#else

	for (int ii = 0; ii < num_lines; ++ii) {
		if (total[ii] != params.n_bitwise_runs) {
			cat->gals[ii].weight *= (double)params.n_bitwise_runs / (double)total[ii];
		}
	}

#endif

	free(total);

#endif

	return(cat);

}
#endif

#ifdef _HAVE_CFITSIO
GalaxyCatalogue * read_data_fits(std::string filename, Cosmology *cosmo, Params params) {
	return(NULL);
}
#endif

GalaxyCatalogue * read_file(std::string filename, std::string file_type, Cosmology *cosmo, Params params) {

	std::cout << "\nReading file: " << filename << "\n";

	GalaxyCatalogue *cat = NULL;
	if (file_type == "ascii") {
		std::cout << "File type ASCII\n";
		cat = read_data_ascii(filename, cosmo, params);
	}
#ifdef _HAVE_HDF5
	else if (file_type == "hdf5") {
		std::cout << "File type hdf5\n";
		cat = read_data_hdf5(filename, cosmo, params);
	}
#endif
#ifdef _HAVE_CFITSIO
	else if (file_type = "fits") {
		cat = read_data_fits(filename, cosmo, params);
	}
#endif
	else {
		std::cout << "Can't find file type for " << filename << "\n";
	}
	std::cout << "\n";

	return(cat);
}

//Takes a filename and appends a number before the file ending
std::string generate_sub_filename(std::string filename, long long sub_number) {

	size_t dot_index = filename.find_last_of('.');
	std::string file_ext = filename.substr(dot_index);

	std::string result = filename;
	result.erase(dot_index);
	result.append(std::to_string(sub_number));
	result.append(file_ext);

	return result;
}

double ls_estimator(double DD, double DR, double RR) {

	if (RR == 0)
		return 0;
	else {
		return DD / RR - 2.0*DR / RR + 1.0;
	}

}

#ifdef _HAVE_HDF5
void write_hdf5_double_attribute(double *attr, std::string attr_name, hid_t loc_id) {

	hid_t attr_id, space_id;
	space_id = H5Screate(H5S_SCALAR);
	attr_id = H5Acreate2(loc_id, attr_name.c_str(), H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, H5T_NATIVE_DOUBLE, attr);
	H5Aclose(attr_id);
	H5Sclose(space_id);

}

void write_hdf5_float_attribute(float *attr, std::string attr_name, hid_t loc_id) {

	hid_t attr_id, space_id;
	space_id = H5Screate(H5S_SCALAR);
	attr_id = H5Acreate2(loc_id, attr_name.c_str(), H5T_NATIVE_FLOAT, space_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, H5T_NATIVE_FLOAT, attr);
	H5Aclose(attr_id);
	H5Sclose(space_id);

}

void write_hdf5_int_attribute(int *attr, std::string attr_name, hid_t loc_id) {

	hid_t attr_id, space_id;
	space_id = H5Screate(H5S_SCALAR);
	attr_id = H5Acreate2(loc_id, attr_name.c_str(), H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, H5T_NATIVE_INT, attr);
	H5Aclose(attr_id);
	H5Sclose(space_id);

}

void write_hdf5_string_attribute(std::string attr, std::string attr_name, hid_t loc_id) {

	hid_t attr_id, space_id;
	space_id = H5Screate(H5S_SCALAR);
	attr_id = H5Acreate2(loc_id, attr_name.c_str(), H5T_C_S1, space_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, H5T_C_S1, attr.c_str());
	H5Aclose(attr_id);
	H5Sclose(space_id);

}

// Should be written to work in n-d
void write_hdf5_dataset_1d(double *attr, hsize_t n_data, std::string attr_name, hid_t loc_id) {

	hid_t data_id, space_id;
	space_id = H5Screate_simple(1, &n_data, NULL);
	data_id = H5Dcreate2(loc_id, attr_name.c_str(), H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, attr);
	H5Dclose(data_id);
	H5Sclose(space_id);

}

void write_hdf5_dataset_2d(double *attr, hsize_t *n_data, std::string attr_name, hid_t loc_id) {

	hid_t data_id, space_id, data_type;
	data_type = H5Tcopy(H5T_NATIVE_DOUBLE);
	H5Tset_order(data_type, H5T_ORDER_LE);
	space_id = H5Screate_simple(2, n_data, NULL);
	data_id = H5Dcreate2(loc_id, attr_name.c_str(), data_type, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, attr);
	H5Dclose(data_id);
	H5Sclose(space_id);

}
#endif

void output_single_xi_1D_ascii(std::string filename, Axis *ax, double *DD, double *DR, double *RR) {

	std::ofstream file(filename);
	for (int ii = 0; ii < ax->n_bins; ++ii) {

		file << ii << "\t";
		file << ax->bin_centre(ii) << "\t";
		file << ax->bin_width(ii) << "\t";
		file << ls_estimator(DD[ii], DR[ii], RR[ii]) << "\t";
		file << DD[ii] << "\t";
		file << DR[ii] << "\t";
		file << RR[ii] << "\n";
	}
	file.close();

}

void output_xi_1D(std::string filename, std::string file_type, JKHists *DD, JKHists *DR, JKHists *RR, GalaxyCatalogue *data_cat, GalaxyCatalogue *rand_cat) {

	if (file_type == "ascii") {

		// Output full result 
		output_single_xi_1D_ascii(filename, &(DD->full_hist.ax), DD->full_hist.bins, DR->full_hist.bins, RR->full_hist.bins);

		// Output jack-knife results if needed
		for (int jk_reg = 0; jk_reg < DD->n_jk_regions; ++jk_reg) {

			output_single_xi_1D_ascii(generate_sub_filename(filename, jk_reg), &(DD->full_hist.ax), \
				&(DD->sub_hist[jk_reg*DD->full_hist.ax.n_bins]), \
				&(DR->sub_hist[jk_reg*DD->full_hist.ax.n_bins]), \
				&(RR->sub_hist[jk_reg*DD->full_hist.ax.n_bins]));

		}

	}

#ifdef _HAVE_HDF5
	else if (file_type == "hdf5") {

		hid_t file_id, group_id;
		file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		group_id = H5Gcreate2(file_id, "full_result", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		write_hdf5_string_attribute("s", "Axis_name", group_id);

		double *axis_data = (double*)malloc(DD->full_hist.ax.n_bins * sizeof(double));
		for (int ii = 0; ii < DD->full_hist.ax.n_bins; ++ii)
			axis_data[ii] = DD->full_hist.ax.bin_centre(ii);
		write_hdf5_dataset_1d(axis_data, DD->full_hist.ax.n_bins, "Bin_centre", group_id);
		for (int ii = 0; ii < DD->full_hist.ax.n_bins; ++ii)
			axis_data[ii] = DD->full_hist.ax.bin_width(ii);
		write_hdf5_dataset_1d(axis_data, DD->full_hist.ax.n_bins, "Bin_width", group_id);
		free(axis_data);

		write_hdf5_int_attribute(&(DD->n_jk_regions), "n_jk_regions", group_id);

		//write_hdf5_int_attribute(&(data_cat->n_gals), "n_data", group_id);
		//write_hdf5_int_attribute(&(rand_cat->n_gals), "n_randoms", group_id);
		//write_hdf5_double_attribute(&(data_cat->sum_weights), "sum_data_weights", group_id);
		//write_hdf5_double_attribute(&(rand_cat->sum_weights), "sum_random_weights", group_id);

		write_hdf5_dataset_1d(DD->full_hist.bins, DD->full_hist.ax.n_bins, "DD", group_id);
		write_hdf5_dataset_1d(DR->full_hist.bins, DR->full_hist.ax.n_bins, "DR", group_id);
		write_hdf5_dataset_1d(RR->full_hist.bins, RR->full_hist.ax.n_bins, "RR", group_id);

		double *xi = (double*)malloc(DD->full_hist.ax.n_bins * sizeof(double));
		for (int ii = 0; ii < DD->full_hist.ax.n_bins; ++ii)
			xi[ii] = ls_estimator(DD->full_hist.bins[ii], DR->full_hist.bins[ii], RR->full_hist.bins[ii]);
		write_hdf5_dataset_1d(xi, DD->full_hist.ax.n_bins, "xi", group_id);

		H5Gclose(group_id);

		// Output jack-knife results if needed
		for (int jk_reg = 0; jk_reg < DD->n_jk_regions; ++jk_reg) {

			std::string group_name = "jk_reg";
			group_name.append(std::to_string((long long)jk_reg));

			group_id = H5Gcreate2(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			//write_hdf5_int_attribute(&(data_cat->n_gals_jk[jk_reg]), "n_data", group_id);
			//write_hdf5_int_attribute(&(rand_cat->n_gals_jk[jk_reg]), "n_randoms", group_id);
			//write_hdf5_double_attribute(&(data_cat->sum_weights_jk[jk_reg]), "sum_data_weights", group_id);
			//write_hdf5_double_attribute(&(rand_cat->sum_weights_jk[jk_reg]), "sum_random_weights", group_id);

			write_hdf5_dataset_1d(&(DD->sub_hist[jk_reg*DD->full_hist.ax.n_bins]), DD->full_hist.ax.n_bins, "DD", group_id);
			write_hdf5_dataset_1d(&(DR->sub_hist[jk_reg*DD->full_hist.ax.n_bins]), DR->full_hist.ax.n_bins, "DR", group_id);
			write_hdf5_dataset_1d(&(RR->sub_hist[jk_reg*DD->full_hist.ax.n_bins]), RR->full_hist.ax.n_bins, "RR", group_id);

			for (int ii = 0; ii < DD->full_hist.ax.n_bins; ++ii)
				xi[ii] = ls_estimator(DD->sub_hist[jk_reg*DD->full_hist.ax.n_bins + ii], \
					DR->sub_hist[jk_reg*DD->full_hist.ax.n_bins + ii], \
					RR->sub_hist[jk_reg*DD->full_hist.ax.n_bins + ii]);

			write_hdf5_dataset_1d(xi, DD->full_hist.ax.n_bins, "xi", group_id);

			H5Gclose(group_id);

		}

		free(xi);
		H5Fclose(file_id);

	}
#endif
	else {
		std::cout << "Output file type not found\n";
	}

}

void output_single_xi_2D_ascii(std::string filename, Axis *ax_1, Axis *ax_2, double *DD, double *DR, double *RR) {

	std::ofstream file(filename);
	for (int ii = 0; ii < ax_1->n_bins; ++ii) {
		for (int jj = 0; jj < ax_2->n_bins; ++jj) {

			int bin = ii*ax_2->n_bins + jj;

			file << ii << "\t";
			file << jj << "\t";
			file << ax_1->bin_centre(ii) << "\t";
			file << ax_2->bin_centre(jj) << "\t";
			file << ax_1->bin_width(ii) << "\t";
			file << ax_2->bin_width(jj) << "\t";
			file << ls_estimator(DD[bin], DR[bin], RR[bin]) << "\t";
			file << DD[bin] << "\t";
			file << DR[bin] << "\t";
			file << RR[bin] << "\n";

		}
	}
	file.close();

}

void output_xi_2D(std::string filename, std::string file_type, std::string axis_1_name, std::string axis_2_name, \
	JKHists2D *DD, JKHists2D *DR, JKHists2D *RR, GalaxyCatalogue *data_cat, GalaxyCatalogue *rand_cat) {

	if (file_type == "ascii") {

		output_single_xi_2D_ascii(filename, &(DD->full_hist.ax1), &(DD->full_hist.ax2), \
			DD->full_hist.bins, DR->full_hist.bins, RR->full_hist.bins);

		for (int jk_reg = 0; jk_reg < DD->n_jk_regions; ++jk_reg) {

			output_single_xi_2D_ascii(generate_sub_filename(filename, jk_reg), &(DD->full_hist.ax1), &(DD->full_hist.ax2),
				&(DD->sub_hist[jk_reg*DD->full_hist.total_bins]), \
				&(DR->sub_hist[jk_reg*DD->full_hist.total_bins]), \
				&(RR->sub_hist[jk_reg*DD->full_hist.total_bins]));

		}

	}

#ifdef _HAVE_HDF5
	else if (file_type == "hdf5") {

		hid_t file_id, group_id;
		file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		hsize_t n_bins[2];
		n_bins[0] = DD->full_hist.ax1.n_bins;
		n_bins[1] = DD->full_hist.ax2.n_bins;

		group_id = H5Gcreate2(file_id, "full_result", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		write_hdf5_string_attribute(axis_1_name, "Axis_1_name", group_id);
		write_hdf5_string_attribute(axis_2_name, "Axis_2_name", group_id);

		double *axis_data = (double*)malloc(DD->full_hist.ax1.n_bins * sizeof(double));
		for (int ii = 0; ii < DD->full_hist.ax1.n_bins; ++ii)
			axis_data[ii] = DD->full_hist.ax1.bin_centre(ii);
		write_hdf5_dataset_1d(axis_data, DD->full_hist.ax1.n_bins, "Axis_1_bin_centre", group_id);
		for (int ii = 0; ii < DD->full_hist.ax1.n_bins; ++ii)
			axis_data[ii] = DD->full_hist.ax1.bin_width(ii);
		write_hdf5_dataset_1d(axis_data, DD->full_hist.ax1.n_bins, "Axis_1_bin_width", group_id);
		free(axis_data);

		axis_data = (double*)malloc(DD->full_hist.ax2.n_bins * sizeof(double));
		for (int ii = 0; ii < DD->full_hist.ax2.n_bins; ++ii)
			axis_data[ii] = DD->full_hist.ax2.bin_centre(ii);
		write_hdf5_dataset_1d(axis_data, DD->full_hist.ax2.n_bins, "Axis_2_bin_centre", group_id);
		for (int ii = 0; ii < DD->full_hist.ax2.n_bins; ++ii)
			axis_data[ii] = DD->full_hist.ax2.bin_width(ii);
		write_hdf5_dataset_1d(axis_data, DD->full_hist.ax2.n_bins, "Axis_2_bin_width", group_id);
		free(axis_data);

		write_hdf5_int_attribute(&(DD->n_jk_regions), "n_jk_regions", group_id);

		//write_hdf5_int_attribute(&(data_cat->n_gals), "n_data", group_id);
		//write_hdf5_int_attribute(&(rand_cat->n_gals), "n_randoms", group_id);
		//write_hdf5_double_attribute(&(data_cat->sum_weights), "sum_data_weights", group_id);
		//write_hdf5_double_attribute(&(rand_cat->sum_weights), "sum_random_weights", group_id);

		write_hdf5_dataset_2d(DD->full_hist.bins, n_bins, "DD", group_id);
		write_hdf5_dataset_2d(DR->full_hist.bins, n_bins, "DR", group_id);
		write_hdf5_dataset_2d(RR->full_hist.bins, n_bins, "RR", group_id);

		double *xi = (double*)malloc(DD->full_hist.total_bins * sizeof(double));
		for (int ii = 0; ii < DD->full_hist.ax1.n_bins; ++ii) {
			for (int jj = 0; jj < DD->full_hist.ax2.n_bins; ++jj) {
				int bin = ii*DD->full_hist.ax2.n_bins + jj;
				xi[bin] = ls_estimator(DD->full_hist.bins[bin], DR->full_hist.bins[bin], RR->full_hist.bins[bin]);
			}
		}
		write_hdf5_dataset_2d(xi, n_bins, "xi", group_id);

		H5Gclose(group_id);

		// Output jack-knife results if needed
		for (int jk_reg = 0; jk_reg < DD->n_jk_regions; ++jk_reg) {

			std::string group_name = "jk_reg";
			group_name.append(std::to_string((long long)jk_reg));

			group_id = H5Gcreate2(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			//write_hdf5_int_attribute(&(data_cat->n_gals_jk[jk_reg]), "n_data", group_id);
			//write_hdf5_int_attribute(&(rand_cat->n_gals_jk[jk_reg]), "n_randoms", group_id);
			//write_hdf5_double_attribute(&(data_cat->sum_weights_jk[jk_reg]), "sum_data_weights", group_id);
			//write_hdf5_double_attribute(&(rand_cat->sum_weights_jk[jk_reg]), "sum_random_weights", group_id);

			write_hdf5_dataset_2d(&(DD->sub_hist[jk_reg*DD->full_hist.total_bins]), n_bins, "DD", group_id);
			write_hdf5_dataset_2d(&(DR->sub_hist[jk_reg*DR->full_hist.total_bins]), n_bins, "DR", group_id);
			write_hdf5_dataset_2d(&(RR->sub_hist[jk_reg*RR->full_hist.total_bins]), n_bins, "RR", group_id);

			for (int ii = 0; ii < DD->full_hist.ax1.n_bins; ++ii) {
				for (int jj = 0; jj < DD->full_hist.ax2.n_bins; ++jj) {
					int jk_bin = jk_reg*DD->full_hist.total_bins + ii*DD->full_hist.ax2.n_bins + jj;
					int bin = ii*DD->full_hist.ax2.n_bins + jj;
					xi[bin] = ls_estimator(DD->sub_hist[jk_bin], DR->sub_hist[jk_bin], RR->sub_hist[jk_bin]);
				}
			}
			write_hdf5_dataset_2d(xi, n_bins, "xi", group_id);

			H5Gclose(group_id);

		}

		free(xi);
		H5Fclose(file_id);

	}
#endif

	else {
		std::cout << "File type not found\n";
	}

}
