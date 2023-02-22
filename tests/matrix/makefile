all: openmp_matrix_test.exe

openmp_matrix_test.exe: openmp_matrix_test.cpp 
	mpicxx $< -std=c++17  -fopenmp -lgtest  -o $@

run: openmp_matrix_test.exe
	./openmp_matrix_test.exe

clean:
	@/bin/rm -f *.exe

