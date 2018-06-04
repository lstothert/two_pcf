#include "cosmology.h"
#include "dataset.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

double Cosmology::codist_integrand(double omega_m, double omega_l, double one_plus_z) {
	return(1.0 / (sqrt(omega_m*one_plus_z*one_plus_z*one_plus_z + omega_l)));
}

Cosmology::Cosmology(double input_omega_m, double input_h, double input_z_min, double input_z_max, \
	double z_save_step, int z_oversample) :
	omega_m(input_omega_m),
	omega_l(1.0 - input_omega_m),
	h(input_h),
	z_min(input_z_min),
	z_max(input_z_max),
	z_codist( (int)(((z_max - z_min) / z_save_step) + 1), 1, z_save_step, input_z_min),
	z_lumdist((int)(((z_max - z_min) / z_save_step) + 1), 1, z_save_step, input_z_min),
	codist_z((int)(((z_max - z_min) / z_save_step) + 1)),
	lumdist_z((int)(((z_max - z_min) / z_save_step) + 1))
{
  c = 299792.458;

	std::cout << "\n--Initialising Cosmology--\n";
	std::cout << "Omega m: " << omega_m << "\n";
	std::cout << "Omega l: " << omega_l << "\n";
	std::cout << "h: " << h << "\n\n";

	double hubble_dist = c / (h * 100.0);

	double currentz = 0.0, currentcodist = 0.0, prevcodist = 0.0, currentlumdist = 0.0;
	double finezstep = z_save_step / (double)z_oversample;

	int numbins = (int)(((z_max - z_min) / z_save_step) + 1);

	// Loop until reach z_min
	while (currentz <= z_min) {

		prevcodist = currentcodist;
		currentcodist += codist_integrand(omega_m, omega_l, 1.0 + currentz) + \
			codist_integrand(omega_l, omega_m, 1.0 + currentz + finezstep);
		currentz += finezstep;

	}

	currentcodist -= ((currentz - z_min) / finezstep)*(currentcodist - prevcodist);
	currentz = z_min;

	for (int ii = 0; ii < numbins; ++ii) {

		// Save codist value
		z_codist.y[ii] = 0.5*finezstep*hubble_dist*currentcodist;
		z_lumdist.y[ii] = (1 + currentz)*z_codist.y[ii];

		// Loop over the over samplenumber to reach next save value
		for (int jj = 0; jj < z_oversample; ++jj) {

			currentcodist += codist_integrand(omega_m, omega_l, 1.0 + currentz) + \
				codist_integrand(omega_m, omega_l, 1 + currentz + finezstep);
			currentz += finezstep;

		}

	}

	//PRODUCE THE INVERTED COMOVDISTARRAY AND LUMDISTARRAY
	codist_z.x_step = (z_codist.y[numbins-1] - z_codist.y[0])/ (double)(numbins - 1);
	codist_z.x_0 = z_codist.y[0];
	codist_z.x_n = codist_z.x_0 + codist_z.x_step*(numbins - 1);

	lumdist_z.x_step = (z_lumdist.y[numbins - 1] - z_lumdist.y[0]) / (double)(numbins - 1);
	lumdist_z.x_0 = z_lumdist.y[0];
	lumdist_z.x_n = lumdist_z.x_0 + lumdist_z.x_step*(numbins - 1);
	
	currentcodist = z_codist.y[0], currentlumdist = z_lumdist.y[0]; 
	codist_z.y[0] = z_min, lumdist_z.y[0] = z_min;
	int jj = 0, kk = 0;
	for (int ii = 1; ii < numbins - 1; ++ii) {

		currentcodist = ii*codist_z.x_step;
		while (currentcodist > z_codist.y[jj]) {
			++jj;
		}
		codist_z.y[ii] = ((double)jj - 1.0 + (currentcodist - z_codist.y[jj - 1]) / (z_codist.y[jj] - z_codist.y[jj - 1]))*z_save_step;

		currentlumdist = ii*lumdist_z.x_step;
		while (currentlumdist > z_lumdist.y[kk]) {
			++kk;
		}
		lumdist_z.y[ii] = ((double)kk - 1.0 + (currentlumdist - z_lumdist.y[kk - 1]) / (z_lumdist.y[jj] - z_lumdist.y[kk - 1]))*z_save_step;

	}
	codist_z.y[numbins-1] = z_max, lumdist_z.y[numbins-1] = z_max;

}

double Cosmology::comov_dist(double z) {
	return(z_codist.lin_interp(z));
}


double Cosmology::lum_dist(double z) {
	return(z_lumdist.lin_interp(z));
}


double Cosmology::angular_diameter_dist(double z) {
	return(comov_dist(z) / (1.0 + z));
}


double Cosmology::comov_volume(double z, double area) {

	double comovingvolume;
	double skyfrac = area*M_PI / 129600.0;
	double codist = comov_dist(z);

	comovingvolume = skyfrac*4.0*M_PI*codist*codist*codist / 3.0;

	return(comovingvolume);
}


double Cosmology::comov_dist_to_z(double comov_dist) {
	return(codist_z.lin_interp(comov_dist));
}


double Cosmology::lum_dist_to_z(double lum_dist) {
	return(lumdist_z.lin_interp(lum_dist));
}


double Cosmology::app_mag(double absmag, double lumdist, double kcorr) {
	return( absmag + 5.0*log10(lumdist) + 25.0 + kcorr );
}


double Cosmology::abs_mag(double appmag, double lumdist, double kcorr) {
	return( appmag - 5 * log10(lumdist) - 25.0 - kcorr );
}
