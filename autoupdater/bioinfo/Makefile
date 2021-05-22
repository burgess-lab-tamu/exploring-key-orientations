GSL = /home/ioerger/gsl

target: superpose_Calphas

gsl_utils.o: gsl_utils.c gsl_utils.h
	g++ -c -I${GSL} gsl_utils.c -lm

superpose_Calphas: superpose_Calphas.cpp gsl_utils.o
	g++ -std=c++0x -Wno-deprecated -I${GSL} -L${GSL}/lib -o superpose_Calphas superpose_Calphas.cpp gsl_utils.o -lm -lgsl -lgslcblas

temp: superpose_Calphas.cpp gsl_utils.o
	g++ -std=c++0x -Wno-deprecated -I${GSL} -L${GSL}/lib -o temp superpose_Calphas.cpp gsl_utils.o -lm -lgsl -lgslcblas




search: search.cpp
	g++ -o search search.cpp

calcsurf: calcsurf.cpp
	g++ -o calcsurf calcsurf.cpp -lm

test:
	javac test.java

all:
	javac *.java
