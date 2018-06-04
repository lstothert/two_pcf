#include "io.h"
#include "cosmology.h"
#include "galaxy.h"
#include "correlators.h"
#include "correlators_invpweight.h"
#include <iostream>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[])
{

	std::cout << "--TWOPCF correlation function calculator--\n\n";

  if(argc == 2){
    std::cout << "Correct number of arguments passed\n";
  } 
  else {
    std::cout << "Incorrect number of arguments given, exiting\n";
    return(1);
  }

  //Read parameter file
  Params params = Params(argv[1]);
	
  //Set number of threads if using openmp
#ifdef _OPENMP
  if (params.n_threads <= 0) {
	  params.n_threads = omp_get_max_threads();
  }
  omp_set_num_threads(params.n_threads);
  std::cout << "\nUsing " << params.n_threads << " threads.\n";
#endif	

  //Initialise cosmology if needed
  Cosmology *cosmo = NULL;
  if (params.coord_system == "equatorial") 
	  cosmo = new Cosmology(params.omega_m, params.h, params.z_min, params.z_max, 0.0001, 100);
	
  //Read in data and randoms
   GalaxyCatalogue *data = read_file(params.data_filename, params.data_file_type, cosmo, params);
   GalaxyCatalogue *randoms = read_file(params.random_filename, params.random_file_type, cosmo, params);

  if (params.coord_system == "equatorial") 
	  delete cosmo;

#ifdef _USE_INV_WEIGHTS

	calc_angular_upweight(data, randoms, params); // This function will only calculate them if needed

  if (params.plot_monopole)
	  auto_correlate_monopole_invpweights(data, randoms, params);

  if (params.plot_sigma_pi)
	  auto_correlate_sigma_pi_invpweights(data, randoms, params);

  if (params.plot_s_mu)
	  auto_correlate_s_mu_invpweights(data, randoms, params);

#else
  if (params.plot_monopole)
	  auto_correlate_monopole(data, randoms, params);

  if (params.plot_sigma_pi)
	  auto_correlate_sigma_pi(data, randoms, params);

  if (params.plot_s_mu)
	  auto_correlate_s_mu(data, randoms, params);

#endif

  return 0;
}