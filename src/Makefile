### for debugging try:
# OPT=-O0 -g -std=c++11
### for performance related builds use
OPT=-O6 -std=c++11 -fPIC
AR=ar
CXX=g++
INC=-I.

# REQUIREMENTS
# 1. lapacke and GNU GSL must be installed (e.g. in ubuntu: "sudo apt install liblapacke liblapacke-dev libgsl-dev libgsl23 libgslcblas0" )
# 2. Make sure that the path pointing to the static lapacke library is included in the linker command
LIB=-L. -L /usr/lib/x86_64-linux-gnu/\
		-Wl,--start-group \
		-llapacke			\
		-lgsl				\
		-Wl,--end-group

#############################################################################################################
#############################################################################################################

default:	lib

lib:		CI RI TI data util
	$(AR) r libconinter.a			\
	        concentric_interpolation.o	\
	        radial_interpolation.o		\
	        tangential_interpolation.o	\
	        data.o				\
	        util.o

clean:
	rm -rf *.o

CI:	concentric_interpolation.h concentric_interpolation.cxx
	$(CXX) $(OPT) -c concentric_interpolation.cxx $(INC)

RI:	radial_interpolation.h radial_interpolation.cxx
	$(CXX) $(OPT) -c radial_interpolation.cxx $(INC)

TI:	tangential_interpolation.h tangential_interpolation.cxx
	$(CXX) $(OPT) -c tangential_interpolation.cxx $(INC)

data:
	$(CXX) $(OPT) -c data.cxx $(INC)

util:	util.h util.cxx
	$(CXX) $(OPT) -c util.cxx $(INC)

#############################################################################################################
#############################################################################################################
