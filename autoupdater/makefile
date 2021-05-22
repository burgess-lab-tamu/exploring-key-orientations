# need to complile on makoshark for executables that use GSL (interface_pairs)

GSL = ./gsl
SRC = ../src/

interface_sizes: interface_sizes.cpp ${SRC}/triplet_lib.o ${SRC}/gsl_utils.o 
	g++ -std=c++0x -Wno-deprecated -I${GSL} -I${SRC} -L${GSL}/lib -o interface_sizes interface_sizes.cpp ${SRC}/triplet_lib.o ${SRC}/gsl_utils.o -lm -lgsl -lgslcblas




