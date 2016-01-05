# makefile
#
#
#  Directives :
#
#    

SHELL=/bin/bash
FC=gfortran
FCUNI=gfortran
LIB=-L/usr/local/lapack -llapack -lrefblas
INC=
#OPT=-O3 -fopenmp
OPT=-Og -fcheck=all 
WFLAGS= 
SRCDIR=/home/joe/Desktop/work/software/SEM_v2.0/src_2.1/

# specify the directory main directory for this experiment
DIR=$(shell pwd | rev | cut -d "/" -f2- | rev )
INIDIR=${DIR}/init/





# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# Test program list #
TestLL : LinkedListClass.o TestLinkedList.o
	${FC} ${OPT} LinkedListClass.o TestLinkedList.o -o $@

TestLinkedList.o : LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}common/Testing/TestLinkedList.f90 -o $@

# Module list #
ModelPrecision.o :
	${FC} ${OPT} -c ${SRCDIR}common/ModelPrecision.f90 -o $@

ModelFlags.o :
	${FC} ${OPT} -c ${SRCDIR}common/ModelFlags.f90 -o $@

LinkedListClass.o :
	${FC} ${OPT} -c ${SRCDIR}common/LinkedListClass.f90 -o $@
	
HashTableClass.o : LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}common/HashTableClass.f90 -o $@

ConstantsDictionary.o : ModelPrecision.o
	${FC} ${OPT} -c ${SRCDIR}common/ConstantsDictionary.f90 -o $@

CommonRoutines.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}common/CommonRoutines.f90 -o $@

RunParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/RunParamsClass.f90 -o $@

Timing.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/Timing.f90 -o $@


# ---------------------------------------------------------------------- #
#                           ~/src/interp/                                #
# ---------------------------------------------------------------------- #

# Binary generation #
LagInterp1DTest : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Legendre.o \
                  Lagrange_1D_Class.o TestLagrangeInterpolation_1D.o
	${FC} ${OPT} ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Legendre.o \
                Lagrange_1D_Class.o TestLagrangeInterpolation_1D.o -o $@
                   
                  
# Test program list #
TestLagrangeInterpolation_1D.o : ModelPrecision.o CommonRoutines.o Legendre.o Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}interp/Testing/TestLagrangeInterpolation_1D.f90 -o $@

# Module list #
Chebyshev.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}interp/Chebyshev.f90 -o $@

Legendre.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}interp/Legendre.f90 -o $@

Lagrange_1D_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}interp/Lagrange_1D_Class.f90 -o $@

Lagrange_2D_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}interp/Lagrange_2D_Class.f90 -o $@
	
Lagrange_3D_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}interp/Lagrange_3D_Class.f90 -o $@



# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod
























