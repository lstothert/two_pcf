#include <cstdlib>
#include <iostream>
#include "galaxy.h"
#include "boxes.h"
#include <math.h>
#include <vector>
#include "healpix/healpix_base.h"
#include "healpix/vec3.h"

GalaxyBoxes::GalaxyBoxes(Galaxy input_galaxies[], int input_n_gals, Galaxy other_cat[], int other_cat_n_gals, double input_rmax){

	std::cout << "Placing galaxies in 3d boxes\n";
	rmax = input_rmax;

	// Count the number of object with good redshifts and shift them to the start
	n_gals = 0;
	for(int ii = 0; ii < input_n_gals - ii + n_gals; ++ii){
		if (input_galaxies[n_gals].d_comov < 0.0) {
			Galaxy buffer(input_galaxies[n_gals]);
			input_galaxies[n_gals] = input_galaxies[input_n_gals - 1 - ii + n_gals];
			input_galaxies[input_n_gals - 1 - ii + n_gals] = buffer;
		}
		else {
			++n_gals;
		}
	}

	// Find the min and max of the galaxies
	for (int ii = 0; ii < 3; ++ii) {
		xmin[ii] = 1e50;
		xmax[ii] = -1e50;
	}

	for (int ii = 0; ii < n_gals; ++ii) {
		for (int jj = 0; jj < 3; ++jj) {

			if (input_galaxies[ii].x[jj] < xmin[jj])
				xmin[jj] = input_galaxies[ii].x[jj];
			else if (input_galaxies[ii].x[jj] > xmax[jj])
				xmax[jj] = input_galaxies[ii].x[jj];

		}
	}

	for (int ii = 0; ii < other_cat_n_gals; ++ii) {
		for (int jj = 0; jj < 3; ++jj) {

			if (other_cat[ii].x[jj] < xmin[jj] && other_cat[ii].d_comov > 0.0)
				xmin[jj] = other_cat[ii].x[jj];
			else if (other_cat[ii].x[jj] > xmax[jj] && other_cat[ii].d_comov > 0.0)
				xmax[jj] = other_cat[ii].x[jj];

		}
	}

	for (int ii = 0; ii < 3; ++ii)
		n_boxes_3d[ii] = (int)(((xmax[ii] - xmin[ii]) / rmax) + 1);

	n_boxes_total = n_boxes_3d[0] * n_boxes_3d[1] * n_boxes_3d[2];

	std::cout << "Min    Max\n";
	std::cout << xmin[0] << " " << xmax[0] << "\n";
	std::cout << xmin[1] << " " << xmax[1] << "\n";
	std::cout << xmin[2] << " " << xmax[2] << "\n";

	std::cout << "Num Boxes: \n" << n_boxes_3d[0] << " " << n_boxes_3d[1] << " " << n_boxes_3d[2] << "\n";
	std::cout << "Total:" << n_boxes_total << "\n";

	// Find the number of galaxies in each box
	int box_id[3], box_full_id;
	int *n_gals_in_box = (int *)calloc(n_boxes_total, sizeof(int));

	n_boxes_used = 0;

	for (int ii = 0; ii < n_gals; ++ii) {

		box_id[0] = (int)((input_galaxies[ii].x[0] - xmin[0]) / rmax);
		box_id[1] = (int)((input_galaxies[ii].x[1] - xmin[1]) / rmax);
		box_id[2] = (int)((input_galaxies[ii].x[2] - xmin[2]) / rmax);
		box_full_id = box_id[0] * n_boxes_3d[1] * n_boxes_3d[2] + box_id[1] * n_boxes_3d[2] + box_id[2];

		if (n_gals_in_box[box_full_id] == 0)
			++n_boxes_used;		

		++n_gals_in_box[box_full_id];

	}

	std::cout << "Num boxes used: " << n_boxes_used << "\n";

	// Allocate boxes
	boxes = (GalaxyBox *)malloc(n_boxes_used * sizeof(GalaxyBox));

	id_to_box_index = (int*)malloc(n_boxes_total * sizeof(int));
	for (int ii = 0; ii < n_boxes_total; ++ii)
		id_to_box_index[ii] = -1;
	box_start_index = (int*)malloc(n_boxes_used * sizeof(int));

	int box_count = 0;
	int n_gals_count = 0;
	for (int ii = 0; ii < n_boxes_total; ++ii) {

		if (n_gals_in_box[ii] > 0) {

			boxes[box_count].numgals = n_gals_in_box[ii];
			
			// Set to zero for use later
			n_gals_in_box[ii] = 0;

			boxes[box_count].id = ii;
			boxes[box_count].loc3d[0] = ii / (n_boxes_3d[1] * n_boxes_3d[2]);
			boxes[box_count].loc3d[1] = (ii / n_boxes_3d[2]) % n_boxes_3d[1];
			boxes[box_count].loc3d[2] = ii % n_boxes_3d[2];

			box_start_index[box_count] = n_gals_count;
			boxes[box_count].gal = &input_galaxies[n_gals_count];

			id_to_box_index[ii] = box_count;

			n_gals_count += boxes[box_count].numgals;
			++box_count;
		}
	}
		
	// Shuffle galaxies so boxes point to correct place
	bool *correct_position = (bool*)calloc(n_gals, sizeof(bool));
	int box_index;

	int gal_index = 0;
	while (gal_index < n_gals) {

		if (correct_position[gal_index]) {
			++gal_index;
		}
		else {

			box_id[0] = (int)((input_galaxies[gal_index].x[0] - xmin[0]) / rmax);
			box_id[1] = (int)((input_galaxies[gal_index].x[1] - xmin[1]) / rmax);
			box_id[2] = (int)((input_galaxies[gal_index].x[2] - xmin[2]) / rmax);
			box_full_id = box_id[0] * n_boxes_3d[1] * n_boxes_3d[2] + box_id[1] * n_boxes_3d[2] + box_id[2];

			box_index = id_to_box_index[box_full_id];

			if (box_start_index[box_index] + n_gals_in_box[box_index] != gal_index) {

				Galaxy buffer(input_galaxies[gal_index]);
				input_galaxies[gal_index] = boxes[box_index].gal[n_gals_in_box[box_index]];
				boxes[box_index].gal[n_gals_in_box[box_index]] = buffer;

			}

			correct_position[box_start_index[box_index] + n_gals_in_box[box_index]] = 1;
			++n_gals_in_box[box_index];

		}

	}

	free(correct_position);
	free(n_gals_in_box);
	std::cout << "\n";

}

GalaxyBoxes::~GalaxyBoxes() {

	free(box_start_index);
	free(id_to_box_index);
	free(boxes);

}

bool are_neighbours(int *loc_1, int *loc_2) {

	//std::cout << loc_1[0] << " " << loc_2[0] << "\n";
	//std::cout << loc_1[1] << " " << loc_2[1] << "\n";
	//std::cout << loc_1[2] << " " << loc_2[2] << "\n\n";

	if (abs(loc_1[0] - loc_2[0]) <= 1 && \
		abs(loc_1[1] - loc_2[1]) <= 1 && \
		abs(loc_1[2] - loc_2[2]) <= 1) {
		return(1);
	}
	else {
		return(0);
	}

}

//Convert 0 to 26 (Num of neighbours of one box) to the actual box index of the comparison catalogue
int adj_id_to_box_index(const int adj_id,const int *loc_3d,const int *id_to_box_index,const int *n_boxes_3d) {

	int new_loc_3d[3];
	new_loc_3d[0] = loc_3d[0] + adj_id / 9 - 1;
	if (new_loc_3d[0] < 0 || new_loc_3d[0] >= n_boxes_3d[0])
		return(-1);
	new_loc_3d[1] = loc_3d[1] + (adj_id % 9) / 3 - 1;
	if (new_loc_3d[1] < 0 || new_loc_3d[1] >= n_boxes_3d[1])
		return(-1);
	new_loc_3d[2] = loc_3d[2] + (adj_id % 3) - 1;
	if (new_loc_3d[2] < 0 || new_loc_3d[2] >= n_boxes_3d[2])
		return(-1);

	return id_to_box_index[new_loc_3d[0] * n_boxes_3d[1] * n_boxes_3d[2] + new_loc_3d[1] * n_boxes_3d[2] + new_loc_3d[2]];

}

GalaxyPixels::~GalaxyPixels() {

	free(healpix_to_pixel_id);
	free(pixels);

}

GalaxyPixels::GalaxyPixels(Galaxy input_galaxies[], int input_n_gals, int healpix_order, double input_theta_max, bool require_redshift) {

	std::cout << "Placing galaxies in pixels\n";

	theta_max = input_theta_max;
	
	// If redshift required shift galaxies with no redshift to the back
	if (require_redshift) {
		n_gals = 0;
		for (int ii = 0; ii < input_n_gals - ii + n_gals; ++ii) {
			if (input_galaxies[n_gals].d_comov < 0.0) {
				Galaxy buffer(input_galaxies[n_gals]);
				input_galaxies[n_gals] = input_galaxies[input_n_gals - 1 - ii + n_gals];
				input_galaxies[input_n_gals - 1 - ii + n_gals] = buffer;
			}
			else {
				++n_gals;
			}
		}
	}
	else {
		n_gals = input_n_gals;
	}

	std::cout << "N gals used: " << n_gals << "\n";

	// Construct the healpix object
	healpix = Healpix_Base(healpix_order, RING);

	// Find the number of galaxies in each pixel
	int id;
	int *n_gals_in_pixel = (int *)calloc(healpix.Npix(), sizeof(int));

	n_pixels_used = 0;

	for (int ii = 0; ii < n_gals; ++ii) {

		vec3 vec(input_galaxies[ii].x[0], input_galaxies[ii].x[1], input_galaxies[ii].x[2]);
		id = healpix.vec2pix(vec);

		if (n_gals_in_pixel[id] == 0)
			++n_pixels_used;

		++n_gals_in_pixel[id];

	}

	std::cout << "Total pixels on the sky: " << healpix.Npix() << "\n";
	std::cout << "Number of pixels used: " << n_pixels_used << "\n";

	// Allocate boxes
	pixels = (GalaxyPixel *)malloc(n_pixels_used * sizeof(GalaxyPixel));

	healpix_to_pixel_id = (int*)malloc(healpix.Npix() * sizeof(int));
	for (int ii = 0; ii < healpix.Npix(); ++ii)
		healpix_to_pixel_id[ii] = -1;
	int *pixel_start_index = (int*)malloc(n_pixels_used * sizeof(int));

	int pixel_count = 0;
	int n_gals_count = 0;
	for (int ii = 0; ii < healpix.Npix(); ++ii) {

		if (n_gals_in_pixel[ii] > 0) {

			pixels[pixel_count].numgals = n_gals_in_pixel[ii];

			// Set to zero for use later
			n_gals_in_pixel[ii] = 0;

			pixels[pixel_count].healpix = ii;

			pixel_start_index[pixel_count] = n_gals_count;
			pixels[pixel_count].gal = &input_galaxies[n_gals_count];

			healpix_to_pixel_id[ii] = pixel_count;

			n_gals_count += pixels[pixel_count].numgals;
			++pixel_count;
		}
	}

	// Shuffle galaxies so boxes point to correct place
	bool *correct_position = (bool*)calloc(n_gals, sizeof(bool));

	int gal_index = 0;
	while (gal_index < n_gals) {

		if (correct_position[gal_index]) {
			++gal_index;
		}
		else {

			vec3 vec(input_galaxies[gal_index].x[0], input_galaxies[gal_index].x[1], input_galaxies[gal_index].x[2]);
			id = healpix.vec2pix(vec);
			id = healpix_to_pixel_id[id];

			if (pixel_start_index[id] + n_gals_in_pixel[id] != gal_index) {

				Galaxy buffer(input_galaxies[gal_index]);
				input_galaxies[gal_index] = pixels[id].gal[n_gals_in_pixel[id]];
				pixels[id].gal[n_gals_in_pixel[id]] = buffer;

			}

			correct_position[pixel_start_index[id] + n_gals_in_pixel[id]] = 1;
			++n_gals_in_pixel[id];

		}

	}

	free(correct_position);
	free(n_gals_in_pixel);
	free(pixel_start_index);

}
