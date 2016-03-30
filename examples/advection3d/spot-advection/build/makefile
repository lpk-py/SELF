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
#OPT=-O3
OPT=-O0 -fcheck=all -Wall
WFLAGS= 
SRCDIR=../src/

# specify the directory main directory for this experiment
DIR=$(shell pwd | rev | cut -d "/" -f2- | rev )
INIDIR=${DIR}/

echo "SRDIR set to :" ${SRCDIR}
echo "INIDIR set to : " ${INIDIR}




# ---------------------------------------------------------------------- #
#               ~/examples/advection3d/spot-advection                    #
# ---------------------------------------------------------------------- #

ADV3DOBJS = ModelPrecision.o \
            ModelFlags.o \
            ConstantsDictionary.o \
            CommonRoutines.o \
            Timing.o \
            LinkedListClass.o \
            KeyRingClass.o \
            Lagrange_1D_Class.o \
            Lagrange_2D_Class.o \
            Lagrange_3D_Class.o \
            Legendre.o \
            Chebyshev.o \
            NodalStorage_3D_Class.o \
            DGSEMSolutionStorageClass_3D.o \
            RollOffFilter3D_Class.o \
            SurfaceClass_3D.o \
            VectorClass.o \
            NodeClass_3D.o \
            FaceClass.o \
            MappedGeometryClass_3D.o \
            HexElementClass.o \
            HexMeshClass.o \
            Advection3DParamsClass.o \
            Advection3DClass.o


initialize_adv3d : ${ADV3DOBJS} Advection3D_InitialConditions.o
	${FC} ${OPT} ${ADV3DOBJS} Advection3D_InitialConditions.o -o $@

integrate_adv3d : ${ADV3DOBJS} Advection3D_Driver.o
	${FC} ${OPT} ${ADV3DOBJS} Advection3D_Driver.o -o $@
         
Advection3DParamsClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}highend/advection3d/Advection3DParamsClass.f90 -o $@

Advection3DClass.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o Timing.o \
                   Lagrange_3D_Class.o NodalStorage_3D_Class.o DGSEMSolutionStorageClass_3D.o \
                   FaceClass.o HexElementClass.o HexMeshClass.o Advection3DParamsClass.o
	${FC} ${OPT} -c ${INIDIR}/localmods/Advection3DClass.f90 -o $@

Advection3D_InitialConditions.o : ModelPrecision.o HexMeshClass.o Advection3DParamsClass.o Advection3DClass.o
	${FC} ${OPT} -c ${INIDIR}/localmods/Advection3D_InitialConditions.f90 -o $@

Advection3D_Driver.o : ModelPrecision.o Timing.o HexMeshClass.o Advection3DParamsClass.o Advection3DClass.o
	${FC} ${OPT} -c ${INIDIR}/localmods/Advection3D_Driver.f90 -o $@



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

RollOffFilter1D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                          Legendre.o Lagrange_1D_Class.o NodalStorage_1D_Class.o
	${FC} ${OPT} -c ${SRCDIR}filters/RollOffFilter1D_Class.f90 -o $@

RollOffFilter3D_Class.o : ModelPrecision.o ModelFlags.o ConstantsDictionary.o CommonRoutines.o \
                          Legendre.o Lagrange_3D_Class.o NodalStorage_3D_Class.o
	${FC} ${OPT} -c ${SRCDIR}filters/RollOffFilter3D_Class.f90 -o $@

# TestPrograms
TestFilter_1D.o : ModelPrecision.o ModelFlags.o NodalStorage_1D_Class.o RollOffFilter1D_Class.o 
	${FC} ${OPT} -c ${SRCDIR}filters/Testing/TestFilter_1D.f90 -o $@

TestFilter_3D.o : ModelPrecision.o ModelFlags.o NodalStorage_3D_Class.o RollOffFilter3D_Class.o 
	${FC} ${OPT} -c ${SRCDIR}filters/Testing/TestFilter_3D.f90 -o $@


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



# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod
























