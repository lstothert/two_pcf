#pragma once

class DatasetFixedXY {

public:

	int len;

	double x_step, x_0, x_n;

	double *y;

	DatasetFixedXY(int input_len, bool is_null = 1, double input_x_step = 0.0, double input_x_0 = 0.0, double *input_y = 0);

	int lower_index_search(double input_x);
	double lin_interp(double input_x);

};

class DatasetXY {

public:

	int len;

	double *x;
	double *y;

	DatasetXY(int input_len, bool is_null = 1, double *input_x = 0, double *input_y = 0);

	int lower_index_search(double input_x);
	double lin_interp(double input_x);

};