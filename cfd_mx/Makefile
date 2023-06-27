ALL: EXE

OCL_FLAG = -lOpenCL

STD_FLAG = -std=c++17

OMP_ANAY = -qopt-report= 5 -qopt-report-phase=vec 

# ------------------------
CC = mpicxx 						# !for amgcl
OPT = $(STD_FLAG) $(STD_FLAG)  -fopenmp -O3 -march=native
# -------------------------------

# CFLAGS= -std=c++17 -g -lOpenCL -fopenmp
# OPT = -O3 -g -fopenmp -std=c++17 -march=native -mcmodel=large# for g++(debug)
# OPT = -O3 -fopenmp -std=c++17 -march=native -mcmodel=large# for g++


# !====================INC

Linking =

INC = -I inc/ \
	$(INCBOOST) \
	-I inc/ocl \
	-I inc/matrix \
	-I inc/Dfib \
	-I inc/Debug \
	-I /mnt/c/Users/user/includes/amgcl \
	-I /mnt/c/Users/user/includes

INCBOOST = -I /home/Jimmy/installBoost/boost_1_76_0

# !====================INC


# -ipo -xCORE-AVX2 -align array32byte 
# -mcmodel=large
# -heap-arrays 64 
# -mcmodel=medium 


DIR_OBJ = ./obj
DIR_BIN = ./bin


# Defining the objects (OBJS) variables
OBJS = \
	main.o \
	# RUN_GPU.o \



# 2_0_2_0_GridGenerator

# Linking object files
EXE = mx

exe :  $(OBJS) 
	$(CC) -o $(EXE) $(OBJS) $(OPT) $(Linking)

# echo something
	@echo "    |  Author:  | Chen Jing Ming                  "
	@echo "    |  Version: | 2.0                             "
	@echo "    |  Web:     | http://smetana.me.ntust.edu.tw/ "
	@echo "    |  Editor:  | JM                              "
	@echo "    |  run:     | ./mx                            "

# Defining the flags of objects
%.o: src/%.cpp
	@$(CC) $< $(OPT) $(INC) -c 

%.o: src/%.c
	@$(C) $< $(OPT) $(INC) -c

# Removeing data file 
.PHONY: cleanData
clean :
	@/bin/rm -rf mx_out
	@/bin/rm -rf Information
	@/bin/rm -f *.optrpt
	@echo "clean"


cleanopt :
	@/bin/rm -f *.optrpt
	@echo "clean opt"


# Removing object files
.PHONY: cleanall
cleanall : 
	@/bin/rm -f $(OBJS) $(EXE)  *.mod
	@/bin/rm -rf mx_out
	@/bin/rm -rf Information
	@/bin/rm -f *.optrpt
	@echo "cleanall"

config :
	if [ ! -d obj ] ; then mkdir obj ; fi
	if [ ! -d run ] ; then mkdir bin ; fi

nohup :
	nohup ./mx &>> nohup.out &

.PHONY: perf-test
perf-test :


.PHONY: run
run: mx
	./mx

.PHONY: test
test : 
	valgrind ./mx

# .PHONY:plot 
# plot :
# 	./
