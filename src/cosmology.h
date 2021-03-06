#pragma once
#include "dataset.h"

class Cosmology {

public:
	double omega_m;
	double omega_l;
	double h;
	double c;

	double z_min;
	double z_max;

	Cosmology(double omega_m, double h, double z_min, double z_max, double z_save_step, int z_oversample);
	double comov_dist(double z);
	double lum_dist(double z);
	double angular_diameter_dist(double z);
	double comov_volume(double z, double area);
	double comov_dist_to_z(double comov_dist);
	double lum_dist_to_z(double lum_dist);
	double app_mag(double absmag, double lumdist, double kcorr);
	double abs_mag(double appmag, double lumdist, double kcorr);

private:
	double codist_integrand(double omega_m, double omega_l, double one_plus_z);

	DatasetFixedXY z_codist;
	DatasetFixedXY codist_z;
	DatasetFixedXY z_lumdist;
	DatasetFixedXY lumdist_z;

};
