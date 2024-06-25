.PHONY: clean, time


NAME = test
MAIN_DEPEND = test.cpp src/geometry/sph_gen.h src/MUSCL_base/MUSCL_base.hpp \
 src/geometry/MUSCL_geometry.hpp $(wildcard src/Riemann_solvers/*.hpp) $(wildcard src/physics/*.hpp) 
SOURCE_FILES = test.cpp src/pmp/surface_mesh.cpp \
src/pmp/algorithms/subdivision.cpp src/pmp/algorithms/differential_geometry.cpp
OBJ_FILES = $(SOURCE_FILES:.cpp=.o)
CC = g++
CFLAGS = -std=c++23 -O3 -g -fopenmp


$(NAME): $(OBJ_FILES) $(NAME).o
	$(CC) $(CFLAGS) -o $(NAME) $(OBJ_FILES)

$(NAME).o: $(MAIN_DEPEND)
	$(CC) $(CFLAGS) -c $(NAME).cpp

src/pmp/surface_mesh.o: src/pmp/surface_mesh.h

src/pmp/algorithms/subdivision.o: src/pmp/algorithms/subdivision.h

src/pmp/algorithms/differential_geometry.o: src/pmp/algorithms/differential_geometry.h


clean: 
	rm -f $(NAME)
	rm -f *.o
	rm -f src/pmp/*.o
	rm -f src/pmp/algorithms/*.o
time:
	time ./$(NAME)
