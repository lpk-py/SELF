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
#OPT=-O3 -fopenmp
OPT=-O0 -fcheck=all -Wall
WFLAGS= 
#SRCDIR=/home/joe/Desktop/work/software/SELF/src/
SRCDIR=/nh/nest/u/jschoonover/Software/SELF/src/

# specify the directory main directory for this experiment
DIR=$(shell pwd | rev | cut -d "/" -f2- | rev )
INIDIR=${DIR}/init/





# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# MODULE list #
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

# ---------------------------------------------------------------------- #
#                       ~/src/iterativesolve/                            #
# ---------------------------------------------------------------------- #

IterativeSolvers_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o \
                           IterativeSolversParams_Class.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/IterativeSolvers_Class.f90 -o $@

# Test programs for the iterative solvers
IterativeSolversParams_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/Testing/GMRES/IterativeSolversParams_Class.f90 -o $@

EllipticFV1D_Class.o : ModelPrecision.o ConstantsDictionary.o IterativeSolvers_Class.o \
                       IterativeSolversParams_Class.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/Testing/CG/EllipticFV1D_Class.f90 -o $@

EllipticFV1D_GMRES_Class.o : ModelPrecision.o ConstantsDictionary.o IterativeSolvers_Class.o \
                       IterativeSolversParams_Class.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/Testing/GMRES/EllipticFV1D_GMRES_Class.f90 -o $@

Driver_EllipticFV1D.o : ModelPrecision.o ConstantsDictionary.o IterativeSolvers_Class.o Timing.o \
                        IterativeSolversParams_Class.o EllipticFV1D_GMRES_Class.o
	${FC} ${OPT} -c ${SRCDIR}iterativesolve/Testing/GMRES/Driver_EllipticFV1D.f90 -o $@

CGFV1D_OBJS = ModelPrecision.o \
              ConstantsDictionary.o \
              CommonRoutines.o \
              ModelFlags.o \
              Timing.o \
              IterativeSolversParams_Class.o \
              IterativeSolvers_Class.o \
              EllipticFV1D_GMRES_Class.o \
              Driver_EllipticFV1D.o

solveFV1D : ${CGFV1D_OBJS}
	${FC} ${OPT} ${CGFV1D_OBJS} -o $@
           
# ---------------------------------------------------------------------- #
#                           ~/src/interp/                                #
# ---------------------------------------------------------------------- #

# MODULE list #
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

# MODULE list #
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

# MODULE list #
DGSEMSolutionStorageClass_1D.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}dgsem/DGSEMSolutionStorageClass_1D.f90 -o $@

DGSEMSolutionStorageClass_2D.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}dgsem/DGSEMSolutionStorageClass_2D.f90 -o $@
	
DGSEMSolutionStorageClass_3D.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}dgsem/DGSEMSolutionStorageClass_3D.f90 -o $@
	
# ---------------------------------------------------------------------- #
#                            ~/src/geom/                                 #
# ---------------------------------------------------------------------- #

Geometry2DTest : ${geomtestobj}
	${FC} ${OPT} ${geomtestobj} -o $@
		
TestMappedGeometry_2D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_2D_Class.o \
                          NodalStorage_2D_Class.o CurveClass_2D.o GeometryParamsClass.o MappedGeometryClass_2D.o
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/TestMappedGeometry_2D.f90 -o $@


Geometry3DTest : ${geom3dtestobj}
	${FC} ${OPT} ${geom3dtestobj} -o $@
		
TestMappedGeometry_3D.o : ModelPrecision.o ConstantsDictionary.o Lagrange_3D_Class.o \
                          NodalStorage_3D_Class.o SurfaceClass_3D.o GeometryParamsClass.o MappedGeometryClass_3D.o
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/TestMappedGeometry_3D.f90 -o $@

hexobjs = ModelPrecision.o \
          ModelFlags.o \
          ConstantsDictionary.o \
          CommonRoutines.o \
          LinkedListClass.o \
          KeyRingClass.o \
          Lagrange_1D_Class.o \
          Lagrange_2D_Class.o \
          Lagrange_3D_Class.o \
          Legendre.o \
          Chebyshev.o \
          NodalStorage_3D_Class.o \
          SurfaceClass_3D.o \
          VectorClass.o \
          NodeClass_3D.o \
          FaceClass.o \
          MappedGeometryClass_3D.o \
          HexElementClass.o \
          HexMeshClass.o \
          HexMeshParamsClass.o \
          TestHexMesh.o \

hextest : ${hexobjs}
	${FC} ${OPT} ${hexobjs} -o $@

TestHexMesh.o : ModelFlags.o NodalStorage_3D_Class.o HexMeshParamsClass.o HexMeshClass.o 
	${FC} ${OPT} -c ${SRCDIR}geom/Testing/TestHexMesh.f90 -o $@


# MODULE list #
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

# ---------------------------------------------------------------------- #
#                     ~/src/highend/shallowwater/                        #
# ---------------------------------------------------------------------- #
         
initialize_sw : ${SWOBJS} ShallowWater_InitialConditions.o
	${FC} ${OPT} ${SWOBJS} ShallowWater_InitialConditions.o -o $@

integrate_sw : ${SWOBJS} ShallowWater_Driver.o
	${FC} ${OPT} ${SWOBJS} ShallowWater_Driver.o -o $@
         
SWParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}highend/shallowwater/SWParamsClass.f90 -o $@

ShallowWaterClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                      Lagrange_2D_Class.o NodalStorage_2D_Class.o DGSEMSolutionStorageClass_2D.o \
                      EdgeClass.o QuadElementClass.o QuadMeshClass.o SWParamsClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/shallowwater/ShallowWaterClass.f90 -o $@

ShallowWater_InitialConditions.o : ModelPrecision.o QuadMeshClass.o SWParamsClass.o ShallowWaterClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/shallowwater/Testing/ShallowWater_InitialConditions.f90 -o $@

ShallowWater_Driver.o : ModelPrecision.o QuadMeshClass.o SWParamsClass.o ShallowWaterClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/shallowwater/ShallowWater_Driver.f90 -o $@


# ---------------------------------------------------------------------- #
#                     ~/src/highend/advection2d /                        #
# ---------------------------------------------------------------------- #

ADVOBJS = ModelPrecision.o \
          ModelFlags.o \
          ConstantsDictionary.o \
          CommonRoutines.o \
          LinkedListClass.o \
          HashTableClass.o \
          Lagrange_1D_Class.o \
          Lagrange_2D_Class.o \
          Legendre.o \
          Chebyshev.o \
          NodalStorage_2D_Class.o \
          DGSEMSolutionStorageClass_2D.o \
          GeometrySupportRoutines.o \
          CurveClass_2D.o \
          VectorClass.o \
          NodeClass_2D.o \
          EdgeClass.o \
          MappedGeometryClass_2D.o \
          QuadElementClass.o \
          QuadMeshClass.o \
          AdvectionParamsClass.o \
          AdvectionClass.o


initialize_adv : ${ADVOBJS} Advection_InitialConditions.o
	${FC} ${OPT} ${ADVOBJS} Advection_InitialConditions.o -o $@

integrate_adv : ${ADVOBJS} Advection_Driver.o
	${FC} ${OPT} ${ADVOBJS} Advection_Driver.o -o $@
         
AdvectionParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection2d/AdvectionParamsClass.f90 -o $@

AdvectionClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                   Lagrange_2D_Class.o NodalStorage_2D_Class.o DGSEMSolutionStorageClass_2D.o \
                   EdgeClass.o QuadElementClass.o QuadMeshClass.o AdvectionParamsClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection2d/AdvectionClass.f90 -o $@

Advection_InitialConditions.o : ModelPrecision.o QuadMeshClass.o AdvectionParamsClass.o AdvectionClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection2d/Advection_InitialConditions.f90 -o $@

Advection_Driver.o : ModelPrecision.o QuadMeshClass.o AdvectionParamsClass.o AdvectionClass.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection2d/Advection_Driver.f90 -o $@
# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod
























