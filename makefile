.PHONY: clean, time


NAME = test
MAIN_DEPEND = test.cpp sph_gen.h MUSCL_base.hpp MUSCL_geometry.hpp HLLE.hpp HLLE_p.hpp HLLC.hpp HLLCplus.hpp
SOURCE_FILES = test.cpp pmp/surface_mesh.cpp pmp/algorithms/subdivision.cpp pmp/algorithms/differential_geometry.cpp
OBJ_FILES = $(SOURCE_FILES:.cpp=.o)
CC = g++
CFLAGS = -std=c++23 -O3 -g -fopenmp


$(NAME): $(OBJ_FILES) $(NAME).o
	$(CC) $(CFLAGS) -o $(NAME) $(OBJ_FILES)

$(NAME).o: $(MAIN_DEPEND)
	$(CC) $(CFLAGS) -c $(NAME).cpp

pmp/surface_mesh.o: pmp/surface_mesh.h

pmp/algorithms/subdivision.o: pmp/algorithms/subdivision.h

pmp/algorithms/differential_geometry.o: pmp/algorithms/differential_geometry.h


clean: 
	rm -f $(NAME)
	rm -f *.o
	rm -f pmp/*.o
	rm -f pmp/algorithms/*.o
time:
	time ./$(NAME)
