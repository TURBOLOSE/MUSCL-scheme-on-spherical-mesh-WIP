.PHONY: clean, time

MAIN_DEPEND = test.cpp  sph_gen.h MUSCL_base.hpp MUSCL_geometry.hpp
SOURCE_FILES = test.cpp pmp/surface_mesh.cpp pmp/algorithms/subdivision.cpp pmp/algorithms/differential_geometry.cpp
OBJ_FILES = $(SOURCE_FILES:.cpp=.o)
CC = g++
CFLAGS = -std=c++20 -O3

main: $(OBJ_FILES) main.o
	$(CC) $(CFLAGS) -o test $(OBJ_FILES)

main.o: $(MAIN_DEPEND)
	$(CC) $(CFLAGS) -c test.cpp

pmp/surface_mesh.o: pmp/surface_mesh.h

pmp/algorithms/subdivision.o: pmp/algorithms/subdivision.h

pmp/algorithms/differential_geometry.o: pmp/algorithms/differential_geometry.h


clean: 
	rm -f test
	rm -f *.o
	rm -f pmp/*.o
	rm -f pmp/algorithms/*.o
time: run_test
	time ./run_test
