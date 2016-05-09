# makefile
#
#
#  Directives :
#
#    

SHELL=/bin/bash
FC=gfortran
FCUNI=gfortran
#LIB=-L/usr/local/lapack -llapack -lrefblas
LIB=
INC=
OPT=-O3
#OPT=-O0 -fcheck=all -Wall
WFLAGS= 
SRCDIR=../src/

# specify the directory main directory for this experiment
DIR=$(shell pwd | rev | cut -d "/" -f2- | rev )
INIDIR=${DIR}/init/



# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# Modules
ModelPrecision.o :
	${FC} ${OPT} -c ${SRCDIR}common/ModelPrecision.f90 -o $@

ModelFlags.o :
	${FC} ${OPT} -c ${SRCDIR}common/ModelFlags.f90 -o $@

LinkedListClass.o :
	${FC} ${OPT} -c ${SRCDIR}common/LinkedListClass.f90 -o $@

KeyRingClass.o : CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/KeyRingClass.f90 -o $@

HashTableClass.o : ConstantsDictionary.o LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}common/HashTableClass.f90 -o $@

ConstantsDictionary.o : ModelPrecision.o
	${FC} ${OPT} -c ${SRCDIR}common/ConstantsDictionary.f90 -o $@

CommonRoutines.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}common/CommonRoutines.f90 -o $@

RunParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/RunParamsClass.f90 -o $@

Timing.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/Timing.f90 -o $@

# Test programs

# ---------------------------------------------------------------------- #
#                       ~/src/iterativesolve/                            #
# ---------------------------------------------------------------------- #

# Modules
IterativeSolvers_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o \
                           IterativeSolversParams_Class.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/IterativeSolvers_Class.f90 -o $@

# Test programs
IterativeSolversParams_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/IterativeSolversParams_Class.f90 -o $@

EllipticFV1D_Class.o : ModelPrecision.o ConstantsDictionary.o IterativeSolvers_Class.o \
                       IterativeSolversParams_Class.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/Testing/CG/EllipticFV1D_Class.f90 -o $@

EllipticFV1D_GMRES_Class.o : ModelPrecision.o ConstantsDictionary.o IterativeSolvers_Class.o \
                       IterativeSolversParams_Class.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/Testing/GMRES/EllipticFV1D_GMRES_Class.f90 -o $@

Driver_EllipticFV1D.o : ModelPrecision.o ConstantsDictionary.o IterativeSolvers_Class.o Timing.o \
                        IterativeSolversParams_Class.o EllipticFV1D_GMRES_Class.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/Testing/GMRES/Driver_EllipticFV1D.f90 -o $@
           
# ---------------------------------------------------------------------- #
#                           ~/src/interp/                                #
# ---------------------------------------------------------------------- #

# Modules
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

# Test Programs

# ---------------------------------------------------------------------- #
#                            ~/src/nodal/                                #
# ---------------------------------------------------------------------- #

# Modules
NodalStorage_1D_Class.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o Legendre.o \
                          Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}nodal/NodalStorage_1D_Class.f90 -o $@

NodalStorage_2D_Class.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o Legendre.o \
                          Lagrange_1D_Class.o Lagrange_2D_Class.o
	${FC} ${OPT} -c ${SRCDIR}nodal/NodalStorage_2D_Class.f90 -o $@

NodalStorage_3D_Class.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o Legendre.o \
                          Lagrange_1D_Class.o Lagrange_3D_Class.o
	${FC} ${OPT} -c ${SRCDIR}nodal/NodalStorage_3D_Class.f90 -o $@

# Test Programs

# ---------------------------------------------------------------------- #
#                          ~/src/filters/                                #
# ---------------------------------------------------------------------- #

# Modules
ModalCutoffFilter1D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                                  Legendre.o Lagrange_1D_Class.o NodalStorage_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}filters/ModalCutoffFilter1D_Class.f90 -o $@

ModalCutoffFilter2D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                          Legendre.o Lagrange_2D_Class.o NodalStorage_2D_Class.o
	${FC} ${OPT} -c ${SRCDIR}filters/ModalCutoffFilter2D_Class.f90 -o $@

RollOffFilter1D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                          Legendre.o Lagrange_1D_Class.o NodalStorage_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}filters/RollOffFilter1D_Class.f90 -o $@

RollOffFilter2D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                          Legendre.o Lagrange_2D_Class.o NodalStorage_2D_Class.o
	${FC} ${OPT} -c ${SRCDIR}filters/RollOffFilter2D_Class.f90 -o $@

RollOffFilter3D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                          Legendre.o Lagrange_3D_Class.o NodalStorage_3D_Class.o
	${FC} ${OPT} -c ${SRCDIR}filters/RollOffFilter3D_Class.f90 -o $@

# TestPrograms
TestFilter_1D.o : ModelPrecision.o ModelFlags.o NodalStorage_1D_Class.o RollOffFilter1D_Class.o 
	${FC} ${OPT} -c ${SRCDIR}filters/Testing/TestFilter_1D.f90 -o $@

TestFilter_2D.o : ModelPrecision.o ModelFlags.o NodalStorage_2D_Class.o RollOffFilter2D_Class.o 
	${FC} ${OPT} -c ${SRCDIR}filters/Testing/TestFilter_2D.f90 -o $@

TestFilter_3D.o : ModelPrecision.o ModelFlags.o NodalStorage_3D_Class.o RollOffFilter3D_Class.o 
	${FC} ${OPT} -c ${SRCDIR}filters/Testing/TestFilter_3D.f90 -o $@


TFOBJS = ModelPrecision.o \
         ModelFlags.o \
         ConstantsDictionary.o \
         CommonRoutines.o \
         Legendre.o \
         Lagrange_1D_Class.o \
         Lagrange_2D_Class.o \
         NodalStorage_2D_Class.o \
         RollOffFilter2D_Class.o \
         TestFilter_2D.o \

tf2d : ${TFOBJS}
	${FC} ${OPT} ${TFOBJS} -o $@

# ---------------------------------------------------------------------- #
#                            ~/src/dgsem/                                #
# ---------------------------------------------------------------------- #

# Modules
DGSEMSolutionStorageClass_1D.o : ModelPrecision.o ConstantsDictionary.o NodalStorage_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}dgsem/DGSEMSolutionStorageClass_1D.f90 -o $@

DGSEMSolutionStorageClass_2D.o : ModelPrecision.o ConstantsDictionary.o NodalStorage_2D_Class.o
	${FC} ${OPT} -c ${SRCDIR}dgsem/DGSEMSolutionStorageClass_2D.f90 -o $@
	
DGSEMSolutionStorageClass_3D.o : ModelPrecision.o ConstantsDictionary.o NodalStorage_3D_Class.o
	${FC} ${OPT} -c ${SRCDIR}dgsem/DGSEMSolutionStorageClass_3D.f90 -o $@

# Test Programs

# ---------------------------------------------------------------------- #
#                            ~/src/cgsem/                                #
# ---------------------------------------------------------------------- #

# Modules
Elliptic1DParams_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}cgsem/Elliptic1DParams_Class.f90 -o $@

Elliptic1D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o Lagrange_1D_Class.o \
                     NodalStorage_1D_Class.o MappedGeometryClass_1D.o SegmentMeshClass.o \
                     IterativeSolvers_Class.o Elliptic1DParams_Class.o 
	${FC} ${OPT} -c ${SRCDIR}cgsem/Elliptic1D_Class.f90 -o $@

CGsemParams_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}cgsem/CGsemParams_Class.f90 -o $@ 

CGsemElliptic_2D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                           NodalStorage_2D_Class.o QuadMeshClass.o MappedGeometryClass_2D.o \
                           CGsemParams_Class.o CGsemPreconditioner_2D_Class.o 
	${FC} ${OPT} -c ${SRCDIR}cgsem/CGsemElliptic_2D_Class.f90 -o $@

CGsemPreconditioner_2D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                                 NodalStorage_2D_Class.o QuadMeshClass.o MappedGeometryClass_2D.o
	${FC} ${OPT} -c ${SRCDIR}cgsem/CGsemPreconditioner_2D_Class.f90 -o $@

# Test Programs
Elliptic1D_Driver.o : ModelPrecision.o IterativeSolvers_Class.o Elliptic1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}cgsem/Elliptic1D_Driver.f90 -o $@

# ---------------------------------------------------------------------- #
#                            ~/src/geom/                                 #
# ---------------------------------------------------------------------- #

# Modules 
NodeClass_1D.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}geom/NodeClass_1D.f90 -o $@
	
NodeClass_2D.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o  LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}geom/NodeClass_2D.f90 -o $@
	
NodeClass_3D.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o LinkedListClass.o
	${FC} ${OPT} -c ${SRCDIR}geom/NodeClass_3D.f90 -o $@

EdgeClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/EdgeClass.f90 -o $@
	
FaceClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/FaceClass.f90 -o $@
	
CurveClass_2D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/CurveClass_2D.f90 -o $@	
	
CurveClass_3D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/CurveClass_3D.f90 -o $@	

SurfaceClass_3D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_1D_Class.o Lagrange_2D_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/SurfaceClass_3D.f90 -o $@	

VectorClass.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/VectorClass.f90 -o $@

GeometryParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}geom/GeometryParamsClass.f90 -o $@

HexMeshParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/HexMeshParamsClass.f90 -o $@

GeometrySupportRoutines.o : ModelPrecision.o ConstantsDictionary.o Lagrange_1D_Class.o \
                            CurveClass_2D.o CurveClass_3D.o
	${FC} ${OPT} -c ${SRCDIR}geom/GeometrySupportRoutines.f90 -o $@

MappedGeometryClass_1D.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o CommonRoutines.o \
                           Lagrange_1D_Class.o 
	${FC} ${OPT} -c ${SRCDIR}geom/MappedGeometryClass_1D.f90 -o $@

MappedGeometryClass_2D.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o CommonRoutines.o \
                           Lagrange_2D_Class.o CurveClass_2D.o VectorClass.o GeometrySupportRoutines.o
	${FC} ${OPT} -c ${SRCDIR}geom/MappedGeometryClass_2D.f90 -o $@

MappedGeometryClass_3D.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o CommonRoutines.o \
                           Lagrange_3D_Class.o SurfaceClass_3D.o VectorClass.o
	${FC} ${OPT} -c ${SRCDIR}geom/MappedGeometryClass_3D.f90 -o $@

QuadElementClass.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o \
                     Lagrange_2D_Class.o CurveClass_2D.o GeometrySupportRoutines.o MappedGeometryClass_2D.o
	${FC} ${OPT} -c ${SRCDIR}geom/QuadElementClass.f90 -o $@	
	
HexElementClass.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o \
                     Lagrange_3D_Class.o SurfaceClass_3D.o MappedGeometryClass_3D.o
	${FC} ${OPT} -c ${SRCDIR}geom/HexElementClass.f90 -o $@

HexMeshClass.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o LinkedListClass.o KeyRingClass.o \
                 Chebyshev.o Lagrange_3D_Class.o MappedGeometryClass_3D.o HexElementClass.o FaceClass.o \
                 NodeClass_3D.o
	${FC} ${OPT} -c ${SRCDIR}geom/HexMeshClass.f90 -o $@

QuadMeshClass.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o LinkedListClass.o HashTableClass.o \
                  Chebyshev.o Lagrange_2D_Class.o CurveClass_2D.o  MappedGeometryClass_2D.o \
                  QuadElementClass.o EdgeClass.o NodeClass_2D.o
	${FC} ${OPT} -c ${SRCDIR}geom/QuadMeshClass.f90 -o $@	

SegmentMeshClass.o : ModelPrecision.o ConstantsDictionary.o ModelFlags.o \
                     Lagrange_1D_Class.o MappedGeometryClass_1D.o 
	${FC} ${OPT} -c ${SRCDIR}geom/SegmentMeshClass.f90 -o $@

TestMappedGeometry_2D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_2D_Class.o \
                          NodalStorage_2D_Class.o CurveClass_2D.o GeometryParamsClass.o MappedGeometryClass_2D.o
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/TestMappedGeometry_2D.f90 -o $@


# Test Programs
TestSegmentMesh.o : ModelPrecision.o ModelFlags.o NodalStorage_1D_Class.o SegmentMeshClass.o 
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/TestSegmentMesh.f90 -o $@

TestMappedGeometry_3D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_3D_Class.o \
                          NodalStorage_3D_Class.o SurfaceClass_3D.o GeometryParamsClass.o MappedGeometryClass_3D.o
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/TestMappedGeometry_3D.f90 -o $@

TestHexMesh.o : ModelFlags.o NodalStorage_3D_Class.o HexMeshParamsClass.o HexMeshClass.o 
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/TestHexMesh.f90 -o $@

DecomposeQuadMesh.o : ModelPrecision.o ConstantsDictionary.o NodalStorage_2D_Class.o QuadMeshClass.o \
                      GeometryParamsClass.o
	${FC} ${OPT} -c ${SRCDIR}geom/DecomposeQuadMesh.f90 -o $@

DQMOBJ = ModelPrecision.o \
         ConstantsDictionary.o \
         ModelFlags.o \
         CommonRoutines.o \
         Legendre.o \
         Chebyshev.o \
         Lagrange_1D_Class.o \
         Lagrange_2D_Class.o \
         NodalStorage_2D_Class.o \
         LinkedListClass.o \
         HashTableClass.o \
         CurveClass_2D.o \
         CurveClass_3D.o \
         VectorClass.o \
         GeometrySupportRoutines.o \
         MappedGeometryClass_2D.o \
         QuadElementClass.o \
         EdgeClass.o \
         NodeClass_2D.o \
         QuadMeshClass.o \
         GeometryParamsClass.o \
         DecomposeQuadMesh.o
         
dcqm : ${DQMOBJ}
	${FC} ${OPT} ${DQMOBJ} -o $@
# ---------------------------------------------------------------------- #
#                     ~/src/highend/shallowwater/                        #
# ---------------------------------------------------------------------- #
         
# Modules         
SWParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}highend/shallowwater/SWParamsClass.f90 -o $@

ShallowWaterClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                      Lagrange_2D_Class.o NodalStorage_2D_Class.o DGSEMSolutionStorageClass_2D.o \
                      EdgeClass.o QuadElementClass.o QuadMeshClass.o SWParamsClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/shallowwater/ShallowWaterClass.f90 -o $@

# Test Programs
ShallowWater_InitialConditions.o : ModelPrecision.o QuadMeshClass.o SWParamsClass.o ShallowWaterClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/shallowwater/Testing/ShallowWater_InitialConditions.f90 -o $@

ShallowWater_Driver.o : ModelPrecision.o QuadMeshClass.o SWParamsClass.o ShallowWaterClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/shallowwater/ShallowWater_Driver.f90 -o $@

# ---------------------------------------------------------------------- #
#                     ~/src/highend/Burgers1D/                           #
# ---------------------------------------------------------------------- #

# Modules         
BurgersParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}highend/burgers1D/BurgersParamsClass.f90 -o $@

BurgersClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                   Lagrange_1D_Class.o NodalStorage_1D_Class.o DGSEMSolutionStorageClass_1D.o \
                   RollOffFilter1D_Class.o BurgersParamsClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/burgers1D/BurgersClass.f90 -o $@

# Test Programs
Burgers_InitialConditions.o : ModelPrecision.o BurgersParamsClass.o BurgersClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/burgers1D/Burgers_InitialConditions.f90 -o $@

Burgers_Driver.o : ModelPrecision.o BurgersParamsClass.o BurgersClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/burgers1D/Burgers_Driver.f90 -o $@

# ---------------------------------------------------------------------- #
#                     ~/src/highend/advection2d /                        #
# ---------------------------------------------------------------------- #

# Modules
AdvectionParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection2d/AdvectionParamsClass.f90 -o $@

AdvectionClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                   Lagrange_2D_Class.o NodalStorage_2D_Class.o DGSEMSolutionStorageClass_2D.o \
                   EdgeClass.o QuadElementClass.o QuadMeshClass.o AdvectionParamsClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection2d/AdvectionClass.f90 -o $@

# Test Programs
Advection_InitialConditions.o : ModelPrecision.o QuadMeshClass.o AdvectionParamsClass.o AdvectionClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection2d/Advection_InitialConditions.f90 -o $@

Advection_Driver.o : ModelPrecision.o QuadMeshClass.o AdvectionParamsClass.o AdvectionClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection2d/Advection_Driver.f90 -o $@

# ---------------------------------------------------------------------- #
#                     ~/src/highend/advection3d /                        #
# ---------------------------------------------------------------------- #

# Modules
Advection3DParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection3d/Advection3DParamsClass.f90 -o $@

Advection3DClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o Timing.o \
                   Lagrange_3D_Class.o NodalStorage_3D_Class.o DGSEMSolutionStorageClass_3D.o \
                   FaceClass.o HexElementClass.o HexMeshClass.o Advection3DParamsClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection3d/Advection3DClass.f90 -o $@

# Test Programs
Advection3D_InitialConditions.o : ModelPrecision.o HexMeshClass.o Advection3DParamsClass.o Advection3DClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection3d/Advection3D_InitialConditions.f90 -o $@

Advection3D_Driver.o : ModelPrecision.o Timing.o HexMeshClass.o Advection3DParamsClass.o Advection3DClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection3d/Advection3D_Driver.f90 -o $@

# ---------------------------------------------------------------------- #
#                 ~/src/highend/hydrostaticprimitives/                   #
# ---------------------------------------------------------------------- #

# Modules
StackedHexMeshClass.o: ModelPrecision.o ConstantsDictionary.o ModelFlags.o Lagrange_2D_Class.o \
                       Lagrange_3D_Class.o QuadMeshClass.o MappedGeometryClass_3D.o \
                       SurfaceClass_3D.o HexElementClass.o FaceClass.o NodeClass_3D.o HexMeshClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/hydrostaticprimitives/StackedHexMeshClass.f90 -o $@


#Advection3DParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
#	${FC} ${OPT} -c ${SRCDIR}highend/advection3d/Advection3DParamsClass.f90 -o $@

#Advection3DClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o Timing.o \
#                   Lagrange_3D_Class.o NodalStorage_3D_Class.o DGSEMSolutionStorageClass_3D.o \
#                   FaceClass.o HexElementClass.o HexMeshClass.o Advection3DParamsClass.o
#	${FC} ${OPT} -c ${SRCDIR}highend/advection3d/Advection3DClass.f90 -o $@

# Test Programs
#Advection3D_InitialConditions.o : ModelPrecision.o HexMeshClass.o Advection3DParamsClass.o Advection3DClass.o
#	${FC} ${OPT} -c ${SRCDIR}highend/advection3d/Advection3D_InitialConditions.f90 -o $@

#Advection3D_Driver.o : ModelPrecision.o Timing.o HexMeshClass.o Advection3DParamsClass.o Advection3DClass.o
#	${FC} ${OPT} -c ${SRCDIR}highend/advection3d/Advection3D_Driver.f90 -o $@


# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod
























