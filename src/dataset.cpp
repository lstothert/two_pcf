#include <cstdlib>
#include "dataset.h"


DatasetXY::DatasetXY(int input_len, bool is_null, double *input_x, double *input_y):
	len(input_len)
{

	if (is_null) {

		x = (double *)calloc(len, sizeof(double));
		y = (double *)calloc(len, sizeof(double));

	}
	else {

		x = (double *)malloc(len * sizeof(double));
		y = (double *)malloc(len * sizeof(double));

		for (int ii = 0; ii < len; ++ii) {
			x[ii] = input_x[ii];
			y[ii] = input_y[ii];
		}

	}

}


int DatasetXY::lower_index_search(double input_x) {

	if (input_x < x[0] || input_x > x[len - 1])
		return -1;

	int first = 0;
	int last = len - 1;
	int middle = (first + last) / 2;

	while (first < last - 1) {

		if (input_x < x[middle]) {

			last = middle;

		}
		else {

			first = middle;

		}

		middle = (first + last) / 2;

	}

	return(first);

}


double DatasetXY::lin_interp(double input_x) {

	if (input_x <= x[0]) {
		return x[0];
	}
	else if (input_x >= x[len - 1]) {
		return x[len - 1];
	}
	else
	{

		int lower_index = lower_index_search(input_x);
		return y[lower_index] + ((input_x - x[lower_index]) / (x[lower_index+1] - x[lower_index]))*(y[lower_index + 1] - y[lower_index]);

	}

}

DatasetFixedXY::DatasetFixedXY(int input_len, bool is_null, double input_x_step, double input_x_0, double *input_y):
	len(input_len),
	x_step(input_x_step),
	x_0(input_x_0),
	x_n(x_0 + (len-1)*x_step)
{

	if (is_null) {
		y = (double *)calloc(len, sizeof(double));
	}
	else {
		y = (double *)malloc(len * sizeof(double));

		for (int ii = 0; ii < len; ++ii) 
			y[ii] = input_y[ii];

	}
}


int DatasetFixedXY::lower_index_search(double input_x) {

	if (input_x < x_0 || input_x > x_n) {
		return -1;
	}
	else {
		return (int)((input_x - x_0) / x_step);
	}
}


double DatasetFixedXY::lin_interp(double input_x) {

	if (input_x <= x_0) {
		return x_0;
	}
	else if (input_x >= x_n){
		return x_n;
	}
	else
	{

		int lower_index = lower_index_search(input_x);
		return y[lower_index] + ((input_x - (x_0 + lower_index*x_step)) / x_step)*(y[lower_index + 1] - y[lower_index]);

	}

}