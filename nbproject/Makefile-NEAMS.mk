#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=NEAMS
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/source/BetaSimulationResults.o \
	${OBJECTDIR}/source/CompositeMicroCell.o \
	${OBJECTDIR}/source/Compound.o \
	${OBJECTDIR}/source/DelayedNeutronSet.o \
	${OBJECTDIR}/source/Element.o \
	${OBJECTDIR}/source/EnumsAndFunctions.o \
	${OBJECTDIR}/source/HomogenousMesh.o \
	${OBJECTDIR}/source/HomogenousMicroCell.o \
	${OBJECTDIR}/source/HomogenousMonteCarlo.o \
	${OBJECTDIR}/source/HomogenousWorth.o \
	${OBJECTDIR}/source/InfiniteCompositeMonteCarlo.o \
	${OBJECTDIR}/source/InfiniteCompositeReactor.o \
	${OBJECTDIR}/source/InfiniteCompositeWorth.o \
	${OBJECTDIR}/source/InfiniteHomogenousReactor.o \
	${OBJECTDIR}/source/InputFileParser.o \
	${OBJECTDIR}/source/Isotope.o \
	${OBJECTDIR}/source/MaterialDataPacket.o \
	${OBJECTDIR}/source/MaterialLibrary.o \
	${OBJECTDIR}/source/Mesh.o \
	${OBJECTDIR}/source/MicroCell.o \
	${OBJECTDIR}/source/MicroCellBoundaryCondition.o \
	${OBJECTDIR}/source/MicroGeometry.o \
	${OBJECTDIR}/source/MicroSolution.o \
	${OBJECTDIR}/source/Mixture.o \
	${OBJECTDIR}/source/PythonPlot.o \
	${OBJECTDIR}/source/ReactivityInsertion.o \
	${OBJECTDIR}/source/Reactor.o \
	${OBJECTDIR}/source/ReactorKinetics.o \
	${OBJECTDIR}/source/ReactorMonteCarlo.o \
	${OBJECTDIR}/source/SimulationResults.o \
	${OBJECTDIR}/source/SphericalMesh.o \
	${OBJECTDIR}/source/Tally.o \
	${OBJECTDIR}/source/TallyGroup.o \
	${OBJECTDIR}/source/WorthStudy.o \
	${OBJECTDIR}/source/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/promptdopplerfeedback

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/promptdopplerfeedback: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/promptdopplerfeedback ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/source/BetaSimulationResults.o: source/BetaSimulationResults.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/BetaSimulationResults.o source/BetaSimulationResults.cpp

${OBJECTDIR}/source/CompositeMicroCell.o: source/CompositeMicroCell.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/CompositeMicroCell.o source/CompositeMicroCell.cpp

${OBJECTDIR}/source/Compound.o: source/Compound.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/Compound.o source/Compound.cpp

${OBJECTDIR}/source/DelayedNeutronSet.o: source/DelayedNeutronSet.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/DelayedNeutronSet.o source/DelayedNeutronSet.cpp

${OBJECTDIR}/source/Element.o: source/Element.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/Element.o source/Element.cpp

${OBJECTDIR}/source/EnumsAndFunctions.o: source/EnumsAndFunctions.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/EnumsAndFunctions.o source/EnumsAndFunctions.cpp

${OBJECTDIR}/source/HomogenousMesh.o: source/HomogenousMesh.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/HomogenousMesh.o source/HomogenousMesh.cpp

${OBJECTDIR}/source/HomogenousMicroCell.o: source/HomogenousMicroCell.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/HomogenousMicroCell.o source/HomogenousMicroCell.cpp

${OBJECTDIR}/source/HomogenousMonteCarlo.o: source/HomogenousMonteCarlo.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/HomogenousMonteCarlo.o source/HomogenousMonteCarlo.cpp

${OBJECTDIR}/source/HomogenousWorth.o: source/HomogenousWorth.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/HomogenousWorth.o source/HomogenousWorth.cpp

${OBJECTDIR}/source/InfiniteCompositeMonteCarlo.o: source/InfiniteCompositeMonteCarlo.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/InfiniteCompositeMonteCarlo.o source/InfiniteCompositeMonteCarlo.cpp

${OBJECTDIR}/source/InfiniteCompositeReactor.o: source/InfiniteCompositeReactor.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/InfiniteCompositeReactor.o source/InfiniteCompositeReactor.cpp

${OBJECTDIR}/source/InfiniteCompositeWorth.o: source/InfiniteCompositeWorth.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/InfiniteCompositeWorth.o source/InfiniteCompositeWorth.cpp

${OBJECTDIR}/source/InfiniteHomogenousReactor.o: source/InfiniteHomogenousReactor.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/InfiniteHomogenousReactor.o source/InfiniteHomogenousReactor.cpp

${OBJECTDIR}/source/InputFileParser.o: source/InputFileParser.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/InputFileParser.o source/InputFileParser.cpp

${OBJECTDIR}/source/Isotope.o: source/Isotope.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/Isotope.o source/Isotope.cpp

${OBJECTDIR}/source/MaterialDataPacket.o: source/MaterialDataPacket.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/MaterialDataPacket.o source/MaterialDataPacket.cpp

${OBJECTDIR}/source/MaterialLibrary.o: source/MaterialLibrary.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/MaterialLibrary.o source/MaterialLibrary.cpp

${OBJECTDIR}/source/Mesh.o: source/Mesh.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/Mesh.o source/Mesh.cpp

${OBJECTDIR}/source/MicroCell.o: source/MicroCell.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/MicroCell.o source/MicroCell.cpp

${OBJECTDIR}/source/MicroCellBoundaryCondition.o: source/MicroCellBoundaryCondition.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/MicroCellBoundaryCondition.o source/MicroCellBoundaryCondition.cpp

${OBJECTDIR}/source/MicroGeometry.o: source/MicroGeometry.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/MicroGeometry.o source/MicroGeometry.cpp

${OBJECTDIR}/source/MicroSolution.o: source/MicroSolution.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/MicroSolution.o source/MicroSolution.cpp

${OBJECTDIR}/source/Mixture.o: source/Mixture.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/Mixture.o source/Mixture.cpp

${OBJECTDIR}/source/PythonPlot.o: source/PythonPlot.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/PythonPlot.o source/PythonPlot.cpp

${OBJECTDIR}/source/ReactivityInsertion.o: source/ReactivityInsertion.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/ReactivityInsertion.o source/ReactivityInsertion.cpp

${OBJECTDIR}/source/Reactor.o: source/Reactor.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/Reactor.o source/Reactor.cpp

${OBJECTDIR}/source/ReactorKinetics.o: source/ReactorKinetics.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/ReactorKinetics.o source/ReactorKinetics.cpp

${OBJECTDIR}/source/ReactorMonteCarlo.o: source/ReactorMonteCarlo.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/ReactorMonteCarlo.o source/ReactorMonteCarlo.cpp

${OBJECTDIR}/source/SimulationResults.o: source/SimulationResults.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/SimulationResults.o source/SimulationResults.cpp

${OBJECTDIR}/source/SphericalMesh.o: source/SphericalMesh.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/SphericalMesh.o source/SphericalMesh.cpp

${OBJECTDIR}/source/Tally.o: source/Tally.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/Tally.o source/Tally.cpp

${OBJECTDIR}/source/TallyGroup.o: source/TallyGroup.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/TallyGroup.o source/TallyGroup.cpp

${OBJECTDIR}/source/WorthStudy.o: source/WorthStudy.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/WorthStudy.o source/WorthStudy.cpp

${OBJECTDIR}/source/main.o: source/main.cpp
	${MKDIR} -p ${OBJECTDIR}/source
	${RM} "$@.d"
	$(COMPILE.cc) -g -DNEAMS -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/source/main.o source/main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
