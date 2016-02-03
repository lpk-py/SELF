# makefile
#
#
#  Directives :
#
#    

SHELL=/bin/bash
FC=gfortran
FCUNI=gfortran
LIB=
INC=
#OPT=-O3 -fopenmp
OPT=-O0 -fcheck=all 
WFLAGS= 
SRCDIR=/nh/nest/u/jschoonover/Software/SELF/src/

# specify the directory main directory for this experiment
DIR=$(shell pwd | rev | cut -d "/" -f2- | rev )
INIDIR=${DIR}/init/





# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# Test program list #
TestLL : LinkedListClass.o TestLinkedList.o
	${FC} ${OPT} LinkedListClass.o TestLinkedList.o -o $@

TestKR : CommonRoutines.o KeyRingClass.o TestKeyRing.o 
	${FC} ${OPT} CommonRoutines.o KeyRingClass.o TestKeyRing.o -o $@

TestLinkedList.o : LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}common/Testing/TestLinkedList.f90 -o $@

TestKeyRing.o : KeyRingClass.o
	${FC} ${OPT} -c ${SRCDIR}common/Testing/TestKeyRing.f90 -o $@

# Module list #
ModelPrecision.o :
	${FC} ${OPT} -c ${SRCDIR}common/ModelPrecision.f90 -o $@

ModelFlags.o :
	${FC} ${OPT} -c ${SRCDIR}common/ModelFlags.f90 -o $@

LinkedListClass.o :
	${FC} ${OPT} -c ${SRCDIR}common/LinkedListClass.f90 -o $@

KeyRingClass.o : CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/KeyRingClass.f90 -o $@

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
























