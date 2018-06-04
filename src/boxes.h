#pragma once
#include "galaxy.h"
#include "healpix/healpix_base.h"

class GalaxyBox {

public:

	int numgals;
	int id;
	int loc3d[3];
	Galaxy *gal;

};

class GalaxyBoxes {

public:

	double xmin[3];
	double xmax[3];
	double rmax;
	int n_boxes_3d[3];
	int n_boxes_total;
	int n_boxes_used;
	int n_gals;

	int *box_start_index;
	GalaxyBox *boxes;

	int *id_to_box_index;

	GalaxyBoxes(Galaxy input_galaxies[], int input_n_gals, Galaxy other_cat[], int other_cat_n_gals, double input_rmax);
	~GalaxyBoxes();

};

bool are_neighbours(int *loc_1, int *loc_2);
int adj_id_to_box_index(const int adj_id,const int *loc_3d,const int *id_to_box_index,const int *n_boxes_3d);

class GalaxyPixel {

public:

	int numgals;
	int healpix;
	Galaxy *gal;

};

class GalaxyPixels {

public:

	int n_pixels_used;
	int n_gals;
	double theta_max;

	int *healpix_to_pixel_id;

	GalaxyPixel *pixels;
	Healpix_Base healpix;

	GalaxyPixels(Galaxy input_galaxies[], int input_n_gals, int healpix_order, double input_theta_max, bool require_redshift);
	~GalaxyPixels();

};