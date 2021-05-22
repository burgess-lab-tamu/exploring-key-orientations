GSL = /pacific/home/ioerger/gsl-1.15

all: triplet_search triplet_search_db loop_search

loop_search: loop_search.cpp triplet_lib.o gsl_utils.o 
	g++ -std=c++0x -Wno-deprecated -o loop_search loop_search.cpp gsl_utils.o triplet_lib.o -lm -lgsl -lgslcblas

gsl_utils.o: gsl_utils.c gsl_utils.h
	g++ -std=c++0x -Wno-deprecated -I${GSL} -c gsl_utils.c

triplet_lib.o: triplet_lib.cpp triplet_lib.hpp
	g++ -std=c++0x -Wno-deprecated -I${GSL} -c triplet_lib.cpp 

triplet_search: triplet_search.cpp triplet_lib.o gsl_utils.o 
	g++ -std=c++0x -Wno-deprecated -I${GSL} -L${GSL}/lib -o triplet_search triplet_search.cpp gsl_utils.o triplet_lib.o -lm -lgsl -lgslcblas

triplet_search_db: triplet_search_db.cpp triplet_lib.o gsl_utils.o 
	g++ -std=c++0x -Wno-deprecated -I${GSL} -L${GSL}/lib -o triplet_search_db triplet_search_db.cpp gsl_utils.o triplet_lib.o -lm -lgsl -lgslcblas




