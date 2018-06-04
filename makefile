# Compile time options

OPTIONS = -D_HAVE_HDF5
OPTIONS += -D_ANGULAR_PI_DEFN
#OPTIONS += -D_USE_INV_WEIGHTS # Turn on for inverse p-weights
#OPTIONS += -D_USE_DITHER_WEIGHTS # Turn on for extra dither weights
OPTIONS += -D_N_MASK_INTS=1

######################

COMP = icc
OPT = -Wall -std=c++11 -openmp $(OPTIONS) -O3 -inline-level=2 -ipo1 -ansi-alias -vec-report=2 -align #icc
#OPT = -Wall -std=c++11 -fopenmp #g
#INCLUDES AND LIBRARIES
INCLUDECOM = -I./src -I./src/healpix
LIB = -lm -lhdf5

#.o FILES
HIST = src/hist.o
COSMO = src/cosmology.o
GALAXY = src/galaxy.o
DATA = src/dataset.o
BOX = src/boxes.o
IO = src/io.o
CORRELATORS = src/correlators.o
CORRELATORS_INVPWEIGHT = src/correlators_invpweight.o

HEALPIX_ERROR = src/healpix/error_handling.o
HEALPIX_BASE = src/healpix/healpix_base.o
HEALPIX_POINTING = src/healpix/pointing.o
HEALPIX_GEOM = src/healpix/geom_utils.o
HEALPIX_TABLES = src/healpix/healpix_tables.o
HEALPIX_STRING = src/healpix/string_utils.o

MAIN = src/main.cpp
OFILES = $(HIST) $(COSMO) $(GALAXY) $(DATA) $(BOX) $(IO) $(CORRELATORS) $(CORRELATORS_INVPWEIGHT) \
         $(HEALPIX_ERROR) $(HEALPIX_BASE) $(HEALPIX_POINTING) \
         $(HEALPIX_GEOM) $(HEALPIX_TABLES) $(HEALPIX_STRING) \
         $(MAIN)

#FINAL GOAL
EXE = TWOPCF

#RULES
default : $(EXE) 
#RULE TO MAKE .o's FROM .cpp's
$(HIST) : src/hist.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(COSMO) : src/cosmology.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(GALAXY) : src/galaxy.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(DATA) : src/dataset.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(IO) : src/io.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(BOX) : src/boxes.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(CORRELATORS) : src/correlators.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(CORRELATORS_INVPWEIGHT) : src/correlators_invpweight.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)

$(HEALPIX_ERROR) : src/healpix/error_handling.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(HEALPIX_BASE) : src/healpix/healpix_base.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(HEALPIX_POINTING) : src/healpix/pointing.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(HEALPIX_GEOM) : src/healpix/geom_utils.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(HEALPIX_TABLES) : src/healpix/healpix_tables.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)
$(HEALPIX_STRING) : src/healpix/string_utils.cpp
	$(COMP) $(OPT) -c $< -o $@ $(INCLUDECOM) $(LIB)

#RULES TO MAKE THE FINAL EXECUTABLES
$(EXE) : $(OFILES)
	$(COMP) $(OPT) $(OFILES) -o $(EXE) $(INCLUDECOM) $(LIB)

#CLEANING RULES
clean :
	rm -f ./src/*.o ./src/healpix/*.o

cleaner :
	rm -f ./src/*.o ./src/*~ *~ ./src/healpix/*.o ./src/healpix/*~ $(EXE)
