all: test.exe test_mpi.exe test_result.exe

test.exe: test.cpp 
	mpicxx $< -std=c++17 -fopenmp -lgtest  -o $@

test_result.exe: test_result.cpp
	mpicxx $< -std=c++17 -fopenmp -lgtest  -o $@

test_mpi.exe: test_mpi.cpp
	mpicxx $< -std=c++17 -fopenmp -lgtest  -o $@

run: test.exe
	rm -rf Information
	rm -rf mx_out
	./test.exe

clean:
	@/bin/rm -f *.mod
	@/bin/rm -rf mx_out
	@/bin/rm -rf Information
	@/bin/rm -f *.optrpt
	@/bin/rm -f *.exe
	@echo "cleanall"
