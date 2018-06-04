#include "cosmology.h"
#include "galaxy.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>


void Galaxy::set(double ra, double dec){

	const double deg_to_rad = M_PI/180.0;

	x_norm[0] = cos(deg_to_rad*dec)*cos(deg_to_rad*ra);
	x_norm[1] = cos(deg_to_rad*dec)*sin(deg_to_rad*ra);
	x_norm[2] = sin(deg_to_rad*dec);

	if (d_comov >= 0.0) {
		x[0] = d_comov * x_norm[0];
		x[1] = d_comov * x_norm[1];
		x[2] = d_comov * x_norm[2];
	}
	else {
		x[0] = 0.0;
		x[1] = 0.0;
		x[2] = 0.0;
	}

}

GalaxyCatalogue::GalaxyCatalogue(int n_gals, int n_jk_regions) : n_jk_regions(n_jk_regions), n_gals(n_gals) {

	gals = new Galaxy[n_gals];

}

GalaxyCatalogue::~GalaxyCatalogue() {

	delete[] gals;

}