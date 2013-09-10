#FC = mpif90
FC = gfortran 
CC = mpicc
CXX = mpiCC
FFLAGS += -O3 
CFLAGS += -O3 

LAPACKLIB = -llapack

NESTLIBDIR = ./
 
FOPT= -O3 -ffixed-line-length-none -ffixed-form
# Determine where to install stuff (prefix is set in configure)
prefix=/scratch/xuewei/study/Dark/darksusy-5.1.1
# DS_INSTALL is where the library and data files will be installed
DS_INSTALL=${prefix}

LIB=$(DS_INSTALL)/lib
INC=-I./ -I$(DS_INSTALL)/include -I$(DS_INSTALL)/src/templates
cfitsio=.


export FC CC CXX FFLAGS CFLAGS LAPACKLIB FOPT LIB INC cfitsio

 
AR = ar r  
LINKLIB = ld -shared
 
NSOBJECTS = utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^ 

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^ 

default: libnest3.a

all: libnest3.a obj_detect eggboxC eggboxC++ gaussian gauss_shell \
rosenbrock himmelblau ackley eptest
 
libnest3.so: $(NSOBJECTS) 
	$(LINKLIB) -o $(LIBS) $@ $^ 
 
libnest3.a: $(NSOBJECTS) 
	$(AR) $@ $^ 
 
obj_detect:
	make -C example_obj_detect
 
gaussian:
	make -C example_gaussian
 
rosenbrock:
	make -C example_rosenbrock
 
ackley:
	make -C example_ackley
 
himmelblau:
	make -C example_himmelblau
 
gauss_shell:
	make -C example_gauss_shell
	
eggboxC:
	make -C example_eggbox_C

eptest:
	make -C eptestF 
	
eggboxC++:
	make -C example_eggbox_C++

clean: 
	-rm $(NESTLIBDIR)/libnest3.*  *.o *.mod
	
cleanall: clean_exec clean clean_obj_detect clean_gaussian clean_gauss_shell \
clean_example_eggbox_C clean_example_eggbox_C++ clean_rosenbrock clean_himmelblau \
clean_ackley clean_eptest

clean_exec:
	-rm obj_detect gaussian rosenbrock ackley himmelblau gauss_shell eggboxC eggboxC++

clean_obj_detect:
	make -C example_obj_detect clean
	
clean_gaussian:
	make -C example_gaussian clean
	
clean_rosenbrock:
	make -C example_rosenbrock clean
	
clean_ackley:
	make -C example_ackley clean
	
clean_himmelblau:
	make -C example_himmelblau clean
	
clean_gauss_shell:
	make -C example_gauss_shell clean
	
clean_example_eggbox_C:
	make -C example_eggbox_C clean

clean_eptest:
	make -C eptestF clean
	
clean_example_eggbox_C++:
	make -C example_eggbox_C++ clean
