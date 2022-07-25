# CC and CFLAGS are varilables
CC = g++
CFLAGS = -Wall -c -std=c++11
# -c option ask g++ to compile the source files, but do not link.
OPTFLAGS = -O3

all	: bin/cswm
	@echo -n ""

# bin/cswm	: main.o Declare.o Init.o Iteration.o Outputfile.o Plot.o
# 			$(CC) $(OPTFLAGS) main.o Declare.o Init.o Outputfile.o Iteration.o Plot.o -I /usr/local/include/matplotlib-cpp-master/ -I /Users/Aaron/miniconda3/lib/python3.8/ -I /Users/Aaron/miniconda3/lib/python3.8/site-packages/numpy/core/include/ -L /usr/local/Frameworks/Python.framework/Versions/3.9/lib -lpython3.9 -o bin/vvm2d
bin/cswm	: main.o Declare.o Init.o Iteration.o Output.o
			$(CC) $(OPTFLAGS) main.o Declare.o Init.o Iteration.o Output.o -o bin/cswm 
main.o 	   	: src/main.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
Declare.o	: src/Declare.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
Init.o		: src/Init.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
Output.o: src/Output.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
# Plot.o		: src/Plot.cpp
# 			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
Iteration.o	: src/Iteration.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@


# clean all the .o and executable files
clean:
		rm -rf *.o bin/*
		rm -rf outputs/h/*  outputs/u/* outputs/v/* outputs/u_lon_lat/* outputs/v_lon_lat/* 
		rm -rf outputs/grid/*
		rm -rf graphs/h/sphere/* graphs/h/curvilinear/* graphs/wind/* graphs/test/h_1/* graphs/test/h_2/* graphs/test/h/*