### for debugging try:
OPT=-O0 -g
### for performance related builds use
# OPT=-O6
AR=ar
CXX=g++
INC=-I. 

# REQUIREMENTS
# 1. lapacke and GNU GSL must be installed (e.g. in ubuntu: "sudo apt-get install liblapacke liblapacke-dev libgsl-dev libgsl23 libgslcblas0" )
# 2. Make sure that the path pointing to the static lapacke library is included in the linker command
LIB=-L. -L /usr/lib/x86_64-linux-gnu/\
		-Wl,--start-group \
		-llapacke			\
		-lgsl				\
		-Wl,--end-group

#############################################################################################################
#############################################################################################################

default:	lib

lib:		util interpolation cubic data data_util  dist_func  obj_func TI RI
	$(AR) r libconinter.a			\
	        cubic_interpolant.o		\
	        concentric_interpolation.o	\
	        util.o				\
	        data.o				\
	        data_util.o			\
	        distance_functions.o		\
	        tangential_interpolation.o	\
	        radial_interpolation.o		\
	        objective_functions.o
	        


clean:
	rm -rf *.o

distclean:
	make clean
	rm -rf *.a

cubic:	cubic_interpolant.h cubic_interpolant.cxx
	$(CXX) $(OPT) -c cubic_interpolant.cxx $(INC)

util:	util.h util.cxx
	$(CXX) $(OPT) -c util.cxx $(INC)

interpolation:	concentric_interpolation.h concentric_interpolation.cxx ../examples/analytical_functions.h ../examples/analytical_functions.cxx
	$(CXX) $(OPT) -c concentric_interpolation.cxx ../examples/analytical_functions.cxx $(INC)

#############################################################################################################
#############################################################################################################

data:
	$(CXX) $(OPT) -c data.cxx $(INC)
	
data_util:
	$(CXX) $(OPT) -c data_util.cxx $(INC)

dist_func:
	$(CXX) $(OPT) -c distance_functions.cxx $(INC)
	
obj_func:
	$(CXX) $(OPT) -c objective_functions.cxx $(INC)

TI:
	$(CXX) $(OPT) -c tangential_interpolation.cxx $(INC)

RI:
	$(CXX) $(OPT) -c radial_interpolation.cxx $(INC)
