# REQUIREMENTS
# lapacke and GNU GSL must be installed (e.g. in ubuntu: "sudo apt-get install liblapacke liblapacke-dev libgsl-dev libgsl23 libgslcblas0" )
INC=-I. -I ../src 
LIB=-L. -L ../src -L/usr/lib/x86_64-linux-gnu/\
		-Wl,--start-group \
		-llapacke			\
		-lgsl				\
		-lconinter				\
		-Wl,--end-group
### for debugging try:
# OPT=-O0 -g
### for performance related builds use
OPT=-O6
CXX=g++
BINDIR=../bin

default:	demo

demo:	interpolation_demo.cxx
	make -C ../src lib
	$(CXX) $(OPT) -c interpolation_demo.cxx $(INC)
	$(CXX) $(OPT) interpolation_demo.o -o $(BINDIR)/interpolation_demo $(LIB)

clean:
	rm -rf *.o

distclean:
	make clean
	make -C ../src distclean
