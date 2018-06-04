Two point correlation function calculator. 
Lee Stothert - 04/06/2018

To do:
- Periodic box
- Cross correlation support
- Non-flat cosmology support
- fits file support
- Input ra dec & comoving instead of ra dec z
- Angular correlation function output (Can be output already using some of the functions but isn't easy)
- Random catalogue generation given mask/nz
- Auto assign jackknife regions

1. Overview ------------------------------------------------

This C++ code calculates the two point galaxy correlation function given a set of data points 
and a set of randoms using the LS estimator. 

It includes:
- Simultaneous/togglable calculation of monopole, 2d split by sigma/pi and by s/mu
- Option to give a weight for each galaxy/random.
- Linear/log (any base) binning in all modes
- File reading in ascii (sep by spaces, tabs or CSV) and hdf5.
- Support for both equatorial and cartesian coordinate systems
- Openmp multithreading
- Bianci and Percival inverse pair weighting scheme
- Angular upweighting
- On the fly jackknife calculations
- python scripts for easy output reading and calculation of wprp/multipoles (Note these may now not work for latest version out of the box)

2. Input ----------------------------------------------------

ASCII columns separated by spaces, tabs or commas accepted.

ASCII cartesian:
x y z <optional weight>

ASCII equatorial:
ra dec z <optional weight>

HDF5: Add the dataset names to the parameter file (Only hdf5 if using inverse p-weights scheme)


3. Example parameter file -------------------------------------------

#Correlation function paramter file. Comments start with #.

data_filename = data.hdf5
data_file_type = hdf5      # ascii/hdf5
random_filename = randoms.hdf5
random_file_type = hdf5    # ascii/hdf5

coord_system = equatorial            # equatorial/cartesian

ra_x_dataset_name = ra        # hdf5 dataset names
dec_y_dataset_name = dec      # ra/dec/z for equatorial
z_z_dataset_name = z          # x/y/z for cartesian
weight_dataset_name = empty    # Name for weight dataset if needed
jk_dataset_name = jk_region

use_weights = 0    # Boolean 0/1, assumes column 4 if reading ascii file
n_threads = 0       # Set to zero for automatic thread detection

n_jk_regions = 0

omega_m = 0.25
h = 1.0
z_min = 0.0
z_max = 1.0

plot_monopole = 0     # Boolean 0/1
monopole_filename = none
monopole_output_type = hdf5
monopole_log_base = 1.3 # Set to 1 for linear, any float above 1.1 valid
monopole_min = 0.0
monopole_max = 100.0
monopole_n_bins = 30

plot_sigma_pi = 0        # Boolean 0/1
sigma_pi_filename = none
sigma_pi_output_type = hdf5
sigma_log_base = 1.0    # Set to 1 for linear, any float above 1.1 valid
sigma_min = 0.0
sigma_max = 50.0
sigma_n_bins = 50
pi_log_base = 1.0             # Set to 1 for linear, any float above 1.1 valid
pi_min = 0.0
pi_max = 50.0
pi_n_bins = 50

plot_s_mu = 1        # Boolean 0/1
s_mu_filename = s_mu.hdf5
s_mu_output_type = hdf5
s_log_base = 1.3      # Set to 1 for linear, any float above 1.1 valid
s_min = 0.0
s_max = 100.0
s_n_bins = 40
mu_n_bins = 50


# All below used in Bianchi Percival inv p weights scheme (Turn on in makefile)

angular_dd_filename = angDD.hdf5
calculate_angular_dd = 0
angular_dd_invpweights_filename = angDDinvpweights.hdf5
calculate_angular_dd_invpweights = 0
angular_dr_filename = angDR.hdf5
calculate_angular_dr = 0
angular_dr_invpweights_filename = angDRinvpweights.hdf5
calculate_angular_dr_invpweights = 0
theta_max = 0.175     # Radians
theta_n_bins = 40
theta_log_base = 1.6
healpix_order = 5
n_bitwise_runs = 2048
bitwise_weight_dataset_name = BITWEIGHT

# Below used if extra dither weights used (Turn on in makefile)
dither_weight_dataset_name = DITHERMASK


4. Output --------------------------------------------------------

Monopole, ascii columns:
<Bin number (0 to n-1)> <Bin centre> <bin width> <xi> <DD> <DR> <RR>

Sigma pi, ascii columns:
<Sigma bin number (0 to n_sigma-1)> <Pi bin number (0 to n_pi-1)> <Sigma bin centre> <Pi bin centre> <Sigma bin width> <Pi bin width> <xi> <DD> <DR> <RR>

S mu, ascii columns:
<S bin number (0 to n_sigma-1)> <mu bin number (0 to n_pi-1)> <S bin centre> <mu bin centre> <S bin width> <mu bin width> <xi> <DD> <DR> <RR>

Can also output in hdf5 where dataset names are simliar to column names above


5. Log binning ---------------------------------------------------

If log base set to < 1.05 the code will use linear bins.
Else the code will use a scheme where each subsequent bin is log_base times the size of the previous bin.

The recommended method for choosing the bin parameters is to choose the min and max scale in the dimension,
then the ratio between the largest and smallest bin is given by log_base^(n_bins - 1). Then vary one of the parameters
roughtly keeping the ratio fixed to find the required bin sizes.

(Will add python script to make thi easier)

6. Inverse p-weighting scheme -----------------------------------

This code supports the Bianchi and Percival 2018 inverse p-weight scheme. This must be turned on with the flag in the makefile _USE_INV_WEIGHTS and extra parameters are needed in the parameter file.

hdf5 input must be used

The parameter bitwise_dataset_name will have integers appended to read in the bitwise masks.
 Multiple bitwise masks will be needed if the number of runs is above 64, the number of bits in a signed int64. The number of these integers needed is set with the makefile option -D_N_MASK_INTS.

So -D_N_MASK_INTS=2 and bitwise_dataset_name=mask will read in two signed int64s from hdf5 datasets mask0 and mask1. The number of realisations doesn't need to be a multiple of 64. The number of bitwise runs is set by the parameter n_bitwise_runs. -D_N_MASK_INTS * 64 must however be larger than n_bitwise runs.

The masks are stored in signed 64 bit integers. So you must first split the array of bits into 64 bit chunks, and convert each of these arrays into the corresponding signed 64 bit integers.

The code will first calculate the angular correlation functions used for angular upweighting. People have often called these functions individually by selecting parameter file flags to turn them on and off and changing the data files as needed. The angular correlation function must be calculated with the parent catalogue, then with the targeted catalogue with invpweights. The DDs and DRs for each are calculated and output ready to be read in by the autocorrelator for the monopole etc. No RRs are calculated but passing RRs as the data and turning on the DDs can give you that histogram to calculate the angular correlaiton function. In the future an option for the angular correlation function will be added. In theory everything can be done at once, by giving the not targeted objects a negative redshift the code can tell which galaxies are targeted and not and use the relevant galaxies for each bit. (This hasn't been fully tested, but it should work and is the easiest way to do things). 

