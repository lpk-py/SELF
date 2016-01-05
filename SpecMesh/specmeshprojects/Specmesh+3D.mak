# Modify these for a given installation
#
F90 = gfortran
SpecMeshSourcePath = /home/joe/Desktop/work/SEM/SpecMesh/specmesh-3d
FTOLibPath = /home/joe/Desktop/work/SEM/SpecMesh/ftobjectlibrary/Source
# end
#
FFLAGS = -O2
##########################
# Object Files for build #
##########################

OBJS = \
3DMeshController.o \
BoundaryEdgeCleaning.o \
ChainedSegmentedCurveClass.o \
CommandLineReader.o \
Connections.o \
CurveConversions.o \
CurveInterpolantClass.o \
ElementOperations.o \
EquationEvaluatorClass.o \
FatalErrorException.o \
FileAndStringProcessing.o \
fmin.o \
FRSegmentedCurveClass.o \
FTDataClass.o \
FTDictionaryClass.o \
FTExceptionClass.o \
FTHashTableClass.o \
FTLinkedListClass.o \
FTObjectArrayClass.o \
FTObjectClass.o \
FTStackClass.o \
FTValueClass.o \
FTValueDictionaryClass.o \
Geometry.o \
Geometry3D.o \
Hash.o \
HexMeshObjects.o \
InterfaceElementMethods.o \
InterpolationAndDerivatives.o \
LaplaceMeshSmoother.o \
Mesh3DOutputMethods.o \
MeshBoundaryMethods.o \
MeshCleaner.o \
MeshGeneratorMethods.o \
MeshOperationsModule.o \
MeshOutputMethods.o \
MeshProject.o \
MeshQualityAnalysis.o \
MeshSmoother.o \
Misc.o \
NodesTemplates.o \
ObjectArrayAdditions.o \
ParametricEquationCurveClass.o \
ProgramGlobals.o \
QuadTreeGridClass.o \
QuadTreeGridGenerator.o \
QuadTreeTemplateOperations.o \
ReaderExceptions.o \
SegmentedCurveArrayClass.o \
SimpleSweep.o \
Sizer.o \
SizerControls.o \
SMChainedCurveClass.o \
SMConstants.o \
SMCurveClass.o \
SMLine.o \
SMMeshClass.o \
SMMeshObjects.o \
SMModel.o \
SMSplineCurveClass.o \
SpecMesh+3DMain.o \
SpringMeshSmoother.o \
Templates.o \
TimerClass.o \
TransfiniteMapClass.o \
Utilities.o \

Specmesh+3D : $(OBJS)
	 ${F90}  -o $@ $(OBJS)

#######################################
# Object dependencies and compilation #
#######################################
3DMeshController.o : $(SpecMeshSourcePath)/3DSource/3DMeshController.f90 \
FTValueDictionaryClass.o \
SimpleSweep.o \
HexMeshObjects.o \
SMMeshObjects.o \
MeshProject.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/3DMeshController.f90

BoundaryEdgeCleaning.o : $(SpecMeshSourcePath)/Mesh/BoundaryEdgeCleaning.f90 \
Connections.o \
ProgramGlobals.o \
MeshBoundaryMethods.o \
SMModel.o \
SMMeshClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/BoundaryEdgeCleaning.f90

ChainedSegmentedCurveClass.o : $(SpecMeshSourcePath)/Curves/DiscreteCurves/ChainedSegmentedCurveClass.f90 \
FTExceptionClass.o \
FTValueClass.o \
FTExceptionClass.o \
FRSegmentedCurveClass.o \
FTObjectArrayClass.o \
FTDataClass.o \
SMConstants.o \
ProgramGlobals.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/DiscreteCurves/ChainedSegmentedCurveClass.f90

CommandLineReader.o : $(SpecMeshSourcePath)/Foundation/CommandLineReader.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Foundation/CommandLineReader.f90

Connections.o : $(SpecMeshSourcePath)/Mesh/Connections.f90 \
SMMeshClass.o \
FTLinkedListClass.o \
SMMeshObjects.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/Connections.f90

CurveConversions.o : $(SpecMeshSourcePath)/Curves/DiscreteCurves/CurveConversions.f90 \
SegmentedCurveArrayClass.o \
SMChainedCurveClass.o \
ChainedSegmentedCurveClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/DiscreteCurves/CurveConversions.f90

CurveInterpolantClass.o : $(SpecMeshSourcePath)/3DSource/Geometry/CurveInterpolantClass.f90 \
InterpolationAndDerivatives.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/Geometry/CurveInterpolantClass.f90

ElementOperations.o : $(SpecMeshSourcePath)/MeshObjects/ElementOperations.f90 \
SMMeshClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/MeshObjects/ElementOperations.f90

EquationEvaluatorClass.o : $(SpecMeshSourcePath)/Curves/ContinuousCurves/EquationEvaluatorClass.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/ContinuousCurves/EquationEvaluatorClass.f90

FatalErrorException.o : $(SpecMeshSourcePath)/Foundation/FatalErrorException.f90 \
FTExceptionClass.o \
FTValueClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Foundation/FatalErrorException.f90

FileAndStringProcessing.o : $(SpecMeshSourcePath)/IO/FileAndStringProcessing.f90 \
SMConstants.o \
ProgramGlobals.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/IO/FileAndStringProcessing.f90

fmin.o : $(SpecMeshSourcePath)/Curves/fmin.f90 \
SMCurveClass.o \
ProgramGlobals.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/fmin.f90

FRSegmentedCurveClass.o : $(SpecMeshSourcePath)/Curves/DiscreteCurves/FRSegmentedCurveClass.f90 \
SMConstants.o \
Geometry.o \
SMCurveClass.o \
ProgramGlobals.o \
FTObjectArrayClass.o \
FTLinkedListClass.o \
ObjectArrayAdditions.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/DiscreteCurves/FRSegmentedCurveClass.f90

FTDataClass.o : $(FTOLibPath)/FTObjects/FTDataClass.f90 \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTDataClass.f90

FTDictionaryClass.o : $(FTOLibPath)/FTObjects/FTDictionaryClass.f90 \
FTObjectArrayClass.o \
FTLinkedListClass.o \
FTObjectClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTDictionaryClass.f90

FTExceptionClass.o : $(FTOLibPath)/FTObjects/FTExceptionClass.f90 \
FTDictionaryClass.o \
FTValueDictionaryClass.o \
FTStackClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTExceptionClass.f90

FTHashTableClass.o : $(FTOLibPath)/FTObjects/FTHashTableClass.f90 \
FTLinkedListClass.o \
FTObjectClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTHashTableClass.f90

FTLinkedListClass.o : $(FTOLibPath)/FTObjects/FTLinkedListClass.f90 \
FTObjectArrayClass.o \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTLinkedListClass.f90

FTObjectArrayClass.o : $(FTOLibPath)/FTObjects/FTObjectArrayClass.f90 \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTObjectArrayClass.f90

FTObjectClass.o : $(FTOLibPath)/FTObjects/FTObjectClass.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTObjectClass.f90

FTStackClass.o : $(FTOLibPath)/FTObjects/FTStackClass.f90 \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTStackClass.f90

FTValueClass.o : $(FTOLibPath)/FTObjects/FTValueClass.f90 \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTValueClass.f90

FTValueDictionaryClass.o : $(FTOLibPath)/FTObjects/FTValueDictionaryClass.f90 \
FTValueClass.o \
FTDictionaryClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTValueDictionaryClass.f90

Geometry.o : $(SpecMeshSourcePath)/Foundation/Geometry.f90 \
SMConstants.o \
ProgramGlobals.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Foundation/Geometry.f90

Geometry3D.o : $(SpecMeshSourcePath)/3DSource/Geometry/Geometry3D.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/Geometry/Geometry3D.f90

Hash.o : $(FTOLibPath)/FTObjects/Hash.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/Hash.f90

HexMeshObjects.o : $(SpecMeshSourcePath)/3DSource/HexMeshObjects.f90 \
FTValueDictionaryClass.o \
SMConstants.o \
SMMeshObjects.o \
MeshProject.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/HexMeshObjects.f90

InterfaceElementMethods.o : $(SpecMeshSourcePath)/Mesh/InterfaceElementMethods.f90 \
MeshProject.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/InterfaceElementMethods.f90

InterpolationAndDerivatives.o : $(SpecMeshSourcePath)/3DSource/Geometry/InterpolationAndDerivatives.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/Geometry/InterpolationAndDerivatives.f90

LaplaceMeshSmoother.o : $(SpecMeshSourcePath)/Mesh/LaplaceMeshSmoother.f90 \
Connections.o \
MeshSmoother.o \
SMModel.o \
SMMeshClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/LaplaceMeshSmoother.f90

Mesh3DOutputMethods.o : $(SpecMeshSourcePath)/3DSource/Mesh3DOutputMethods.f90 \
HexMeshObjects.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/Mesh3DOutputMethods.f90

MeshBoundaryMethods.o : $(SpecMeshSourcePath)/Mesh/MeshBoundaryMethods.f90 \
CurveConversions.o \
FatalErrorException.o \
ProgramGlobals.o \
Connections.o \
SMChainedCurveClass.o \
SMModel.o \
SMMeshClass.o \
fmin.o \
Sizer.o \
Geometry.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/MeshBoundaryMethods.f90

MeshCleaner.o : $(SpecMeshSourcePath)/Mesh/MeshCleaner.f90 \
MeshQualityAnalysis.o \
InterfaceElementMethods.o \
FatalErrorException.o \
Connections.o \
ElementOperations.o \
SMModel.o \
SMMeshClass.o \
fmin.o \
Geometry.o \
MeshBoundaryMethods.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/MeshCleaner.f90

MeshGeneratorMethods.o : $(SpecMeshSourcePath)/Mesh/MeshGeneratorMethods.f90 \
QuadTreeGridGenerator.o \
MeshBoundaryMethods.o \
FatalErrorException.o \
fmin.o \
Geometry.o \
MeshProject.o \
BoundaryEdgeCleaning.o \
MeshOperationsModule.o \
ProgramGlobals.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/MeshGeneratorMethods.f90

MeshOperationsModule.o : $(SpecMeshSourcePath)/Mesh/MeshOperationsModule.f90 \
FTObjectClass.o \
QuadTreeGridClass.o \
FTLinkedListClass.o \
SMMeshObjects.o \
SMMeshClass.o \
FTLinkedListClass.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/MeshOperationsModule.f90

MeshOutputMethods.o : $(SpecMeshSourcePath)/IO/MeshOutputMethods.f90 \
MeshOperationsModule.o \
SMModel.o \
SMMeshObjects.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/IO/MeshOutputMethods.f90

MeshProject.o : $(SpecMeshSourcePath)/Project/MeshProject.f90 \
CurveConversions.o \
SMModel.o \
SMConstants.o \
MeshSmoother.o \
LaplaceMeshSmoother.o \
Sizer.o \
QuadTreeGridClass.o \
SMMeshClass.o \
SpringMeshSmoother.o \
FTExceptionClass.o \
FTValueClass.o \
SMChainedCurveClass.o \
ChainedSegmentedCurveClass.o \
SizerControls.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Project/MeshProject.f90

MeshQualityAnalysis.o : $(SpecMeshSourcePath)/Mesh/MeshQualityAnalysis.f90 \
SMMeshObjects.o \
SMMeshClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/MeshQualityAnalysis.f90

MeshSmoother.o : $(SpecMeshSourcePath)/Mesh/MeshSmoother.f90 \
MeshBoundaryMethods.o \
SMModel.o \
SMMeshClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/MeshSmoother.f90

Misc.o : $(SpecMeshSourcePath)/Foundation/Misc.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Foundation/Misc.f90

NodesTemplates.o : $(SpecMeshSourcePath)/QuadTreeGrid/NodesTemplates.f90 \
ProgramGlobals.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/QuadTreeGrid/NodesTemplates.f90

ObjectArrayAdditions.o : $(SpecMeshSourcePath)/Categories/ObjectArrayAdditions.f90 \
FTObjectArrayClass.o \
FTLinkedListClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Categories/ObjectArrayAdditions.f90

ParametricEquationCurveClass.o : $(SpecMeshSourcePath)/Curves/ContinuousCurves/ParametricEquationCurveClass.f90 \
FTExceptionClass.o \
FTValueClass.o \
SMConstants.o \
EquationEvaluatorClass.o \
SMCurveClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/ContinuousCurves/ParametricEquationCurveClass.f90

ProgramGlobals.o : $(SpecMeshSourcePath)/Foundation/ProgramGlobals.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Foundation/ProgramGlobals.f90

QuadTreeGridClass.o : $(SpecMeshSourcePath)/QuadTreeGrid/QuadTreeGridClass.f90 \
Sizer.o \
SMConstants.o \
SMMeshObjects.o \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/QuadTreeGrid/QuadTreeGridClass.f90

QuadTreeGridGenerator.o : $(SpecMeshSourcePath)/QuadTreeGrid/QuadTreeGridGenerator.f90 \
Sizer.o \
QuadTreeGridClass.o \
QuadTreeTemplateOperations.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/QuadTreeGrid/QuadTreeGridGenerator.f90

QuadTreeTemplateOperations.o : $(SpecMeshSourcePath)/QuadTreeGrid/QuadTreeTemplateOperations.f90 \
QuadTreeGridClass.o \
SMConstants.o \
SMMeshObjects.o \
Templates.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/QuadTreeGrid/QuadTreeTemplateOperations.f90

ReaderExceptions.o : $(SpecMeshSourcePath)/3DSource/ReaderExceptions.f90 \
FTExceptionClass.o \
FTValueDictionaryClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/ReaderExceptions.f90

SegmentedCurveArrayClass.o : $(SpecMeshSourcePath)/Curves/DiscreteCurves/SegmentedCurveArrayClass.f90 \
Geometry.o \
ProgramGlobals.o \
SMConstants.o \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/DiscreteCurves/SegmentedCurveArrayClass.f90

SimpleSweep.o : $(SpecMeshSourcePath)/3DSource/SimpleSweep.f90 \
FTExceptionClass.o \
HexMeshObjects.o \
FTExceptionClass.o \
CurveInterpolantClass.o \
TransfiniteMapClass.o \
FTValueDictionaryClass.o \
MeshProject.o \
ReaderExceptions.o \
SMConstants.o \
ProgramGlobals.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/SimpleSweep.f90

Sizer.o : $(SpecMeshSourcePath)/Project/Sizer/Sizer.f90 \
Geometry.o \
FTLinkedListClass.o \
FTLinkedListClass.o \
SizerControls.o \
ChainedSegmentedCurveClass.o \
SMConstants.o \
ProgramGlobals.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Project/Sizer/Sizer.f90

SizerControls.o : $(SpecMeshSourcePath)/Project/Sizer/SizerControls.f90 \
SMConstants.o \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Project/Sizer/SizerControls.f90

SMChainedCurveClass.o : $(SpecMeshSourcePath)/Curves/ContinuousCurves/SMChainedCurveClass.f90 \
ProgramGlobals.o \
SMCurveClass.o \
FTExceptionClass.o \
FTValueClass.o \
FTObjectArrayClass.o \
FTLinkedListClass.o \
FTExceptionClass.o \
FTObjectClass.o \
FTLinkedListClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/ContinuousCurves/SMChainedCurveClass.f90

SMConstants.o : $(SpecMeshSourcePath)/Foundation/SMConstants.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Foundation/SMConstants.f90

SMCurveClass.o : $(SpecMeshSourcePath)/Curves/ContinuousCurves/SMCurveClass.f90 \
ProgramGlobals.o \
SMConstants.o \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/ContinuousCurves/SMCurveClass.f90

SMLine.o : $(SpecMeshSourcePath)/Curves/ContinuousCurves/SMLine.f90 \
SMCurveClass.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/ContinuousCurves/SMLine.f90

SMMeshClass.o : $(SpecMeshSourcePath)/Project/Mesh/SMMeshClass.f90 \
FTHashTableClass.o \
SegmentedCurveArrayClass.o \
FTObjectArrayClass.o \
FTLinkedListClass.o \
SMMeshObjects.o \
FTObjectClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Project/Mesh/SMMeshClass.f90

SMMeshObjects.o : $(SpecMeshSourcePath)/MeshObjects/SMMeshObjects.f90 \
FTObjectArrayClass.o \
ProgramGlobals.o \
SMConstants.o \
FTLinkedListClass.o \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/MeshObjects/SMMeshObjects.f90

SMModel.o : $(SpecMeshSourcePath)/Project/Model/SMModel.f90 \
FTExceptionClass.o \
SMConstants.o \
FTExceptionClass.o \
SMChainedCurveClass.o \
SMLine.o \
FTLinkedListClass.o \
FTValueClass.o \
ParametricEquationCurveClass.o \
SMSplineCurveClass.o \
ProgramGlobals.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Project/Model/SMModel.f90

SMSplineCurveClass.o : $(SpecMeshSourcePath)/Curves/ContinuousCurves/SMSplineCurveClass.f90 \
SMCurveClass.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Curves/ContinuousCurves/SMSplineCurveClass.f90

SpecMesh+3DMain.o : $(SpecMeshSourcePath)/SpecMesh+3DMain.f90 \
MeshQualityAnalysis.o \
MeshCleaner.o \
FTExceptionClass.o \
Mesh3DOutputMethods.o \
CommandLineReader.o \
TimerClass.o \
MeshProject.o \
MeshOutputMethods.o \
3DMeshController.o \
MeshGeneratorMethods.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/SpecMesh+3DMain.f90

SpringMeshSmoother.o : $(SpecMeshSourcePath)/Mesh/SpringMeshSmoother.f90 \
MeshSmoother.o \
MeshBoundaryMethods.o \
SMModel.o \
SMMeshClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Mesh/SpringMeshSmoother.f90

Templates.o : $(SpecMeshSourcePath)/QuadTreeGrid/Templates.f90 \
QuadTreeGridClass.o \
SMConstants.o \
SMMeshObjects.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/QuadTreeGrid/Templates.f90

TimerClass.o : $(SpecMeshSourcePath)/Foundation/TimerClass.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/Foundation/TimerClass.f90

TransfiniteMapClass.o : $(SpecMeshSourcePath)/3DSource/Geometry/TransfiniteMapClass.f90 \
CurveInterpolantClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/Geometry/TransfiniteMapClass.f90

Utilities.o : $(SpecMeshSourcePath)/3DSource/Geometry/Utilities.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(SpecMeshSourcePath)/3DSource/Geometry/Utilities.f90

