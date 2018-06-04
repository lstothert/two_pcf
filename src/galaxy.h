#pragma once
#include "cosmology.h"

// Only used if equatorial coordinate system chosen
class GalaxyEq {

public:

	double ra;
	double dec;
	double z;

};

// Minimum values needed for all parts of calculation
class Galaxy {

public:
	double x[3];
	float weight;
	int jk_region;
	double d_comov;
	double x_norm[3];

#ifdef _USE_INV_WEIGHTS
	long long bitmask[_N_MASK_INTS];
	float input_weight;
#ifdef _USE_DITHER_WEIGHTS
	long long dither_mask[_N_MASK_INTS];
#endif
#endif

	void set(double ra, double dec);

};

class GalaxyCatalogue {

public:
	int n_jk_regions;
	int n_gals;
	Galaxy *gals;

	GalaxyCatalogue(int n_gals, int n_jk_regions);
	~GalaxyCatalogue();

};
