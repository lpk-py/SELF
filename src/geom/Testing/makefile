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
#OPT=-O3 -fOPENmp
OPT=-Og -fcheck=all 
WFLAGS= 
SRCDIR=/home/joe/Desktop/work/software/SEM_v2.0/src_2.1/

# specify the directory main directory for this experiment
DIR=$(shell pwd | rev | cut -d "/" -f2- | rev )
INIDIR=${DIR}/init/





# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

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

# Module list #
Chebyshev.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}interp/Chebyshev.f90 -o $@

Legendre.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}interp/Legendre.f90 -o $@

Lagrange_1D_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}interp/Lagrange_1D_Class.f90 -o $@

Lagrange_2D_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}interp/Lagrange_2D_Class.f90 -o $@
	
Lagrange_3D_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}interp/Lagrange_3D_Class.f90 -o $@


# ---------------------------------------------------------------------- #
#                            ~/src/nodal/                                #
# ---------------------------------------------------------------------- #

# Module list #
NodalStorage_1D_Class.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o Legendre.o \
                          Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}nodal/NodalStorage_1D_Class.f90 -o $@

NodalStorage_2D_Class.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o Legendre.o \
                          Lagrange_1D_Class.o Lagrange_2D_Class.o
	${FC} ${OPT} -c ${SRCDIR}nodal/NodalStorage_2D_Class.f90 -o $@
	
NodalStorage_3D_Class.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o Legendre.o \
                          Lagrange_1D_Class.o Lagrange_3D_Class.o
	${FC} ${OPT} -c ${SRCDIR}nodal/NodalStorage_3D_Class.f90 -o $@


# ---------------------------------------------------------------------- #
#                            ~/src/dgsem/                                #
# ---------------------------------------------------------------------- #

# Module list #
DGSEMSolutionStorageClass_1D.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}dgsem/DGSEMSolutionStorageClass_1D.f90 -o $@

DGSEMSolutionStorageClass_2D.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}dgsem/DGSEMSolutionStorageClass_2D.f90 -o $@
	
DGSEMSolutionStorageClass_3D.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}dgsem/DGSEMSolutionStorageClass_3D.f90 -o $@
	
# ---------------------------------------------------------------------- #
#                            ~/src/geom/                                 #
# ---------------------------------------------------------------------- #

geomtestobj=ModelPrecision.o \
            ConstantsDictionary.o \
            ModelFlags.o \
            CommonRoutines.o \
            Legendre.o \
            Lagrange_1D_Class.o \
            Lagrange_2D_Class.o \
            NodalStorage_2D_Class.o \
            CurveClass_2D.o \
            VectorClass.o \
            GeometryParamsClass.o \
            GeometrySupportRoutines.o \
            MappedGeometryClass_2D.o \
            TestMappedGeometry_2D.o

Geometry2DTest : ${geomtestobj}
	${FC} ${OPT} ${geomtestobj} -o $@
		
TestMappedGeometry_2D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_2D_Class.o \
                          NodalStorage_2D_Class.o CurveClass_2D.o GeometryParamsClass.o MappedGeometryClass_2D.o
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/TestMappedGeometry_2D.f90 -o $@
	
meshtestobj=ModelPrecision.o \
            ConstantsDictionary.o \
            ModelFlags.o \
            CommonRoutines.o \
            LinkedListClass.o \
            HashTableClass.o \
            Chebyshev.o \
            Legendre.o \
            Lagrange_1D_Class.o \
            Lagrange_2D_Class.o \
            NodalStorage_2D_Class.o \
            CurveClass_2D.o \
            VectorClass.o \
            GeometryParamsClass.o \
            GeometrySupportRoutines.o \
            MappedGeometryClass_2D.o \
            EdgeClass.o \
            NodeClass_2D.o \
            QuadElementClass.o \
            QuadMeshClass.o \
            TestQuadMesh.o
            

QuadMeshTest : ${meshtestobj}
	${FC} ${OPT} ${meshtestobj} -o $@
		
TestQuadMesh.o : ModelPrecision.o ConstantsDictionary.o NodalStorage_2D_Class.o GeometryParamsClass.o QuadMeshClass.o
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/TestQuadMesh.f90 -o $@

# Module list #
NodeClass_1D.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}geom/NodeClass_1D.f90 -o $@
	
NodeClass_2D.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o  LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}geom/NodeClass_2D.f90 -o $@
	
NodeClass_3D.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}geom/NodeClass_3D.f90 -o $@

EdgeClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/EdgeClass.f90 -o $@
	
CurveClass_2D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/CurveClass_2D.f90 -o $@	
	
CurveClass_3D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/CurveClass_3D.f90 -o $@	

VectorClass.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/VectorClass.f90 -o $@

GeometryParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}geom/GeometryParamsClass.f90 -o $@

GeometrySupportRoutines.o : ModelPrecision.o ConstantsDictionary.o Lagrange_1D_Class.o \
                            CurveClass_2D.o CurveClass_3D.o
	${FC} ${OPT} -c ${SRCDIR}geom/GeometrySupportRoutines.f90 -o $@
	
MappedGeometryClass_2D.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o CommonRoutines.o \
                           Lagrange_2D_Class.o CurveClass_2D.o VectorClass.o GeometrySupportRoutines.o
	${FC} ${OPT} -c ${SRCDIR}geom/MappedGeometryClass_2D.f90 -o $@

QuadElementClass.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o \
                     Lagrange_2D_Class.o CurveClass_2D.o GeometrySupportRoutines.o MappedGeometryClass_2D.o
	${FC} ${OPT} -c ${SRCDIR}geom/QuadElementClass.f90 -o $@	

QuadMeshClass.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o LinkedListClass.o HashTableClass.o \
                  Chebyshev.o Lagrange_2D_Class.o CurveClass_2D.o  MappedGeometryClass_2D.o \
                  QuadElementClass.o EdgeClass.o NodeClass_2D.o
	${FC} ${OPT} -c ${SRCDIR}geom/QuadMeshClass.f90 -o $@	

# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod
























