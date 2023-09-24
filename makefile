.PHONY: clean, time
all: run_test 
	./run_test
run_test: test.cpp  sph_gen.h
	g++ test.cpp pmp/surface_mesh.cpp pmp/algorithms/subdivision.cpp pmp/algorithms/differential_geometry.cpp  -std=c++20 -O3 -o run_test

clean: 
	rm -f run_test
time: run_test
	time ./run_test
