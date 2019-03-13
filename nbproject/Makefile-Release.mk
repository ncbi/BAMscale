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
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/BAMcoverage.o \
	${OBJECTDIR}/src/BAMstructs.o \
	${OBJECTDIR}/src/BEDstruct.o \
	${OBJECTDIR}/src/CHROMstruct.o \
	${OBJECTDIR}/src/Inputs.o \
	${OBJECTDIR}/src/Writer.o \
	${OBJECTDIR}/src/binning.o \
	${OBJECTDIR}/src/main.o \
	${OBJECTDIR}/src/multithreads.o \
	${OBJECTDIR}/src/scale.o \
	${OBJECTDIR}/src/segmenter.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk bin/BAMscale

bin/BAMscale: ${OBJECTFILES}
	${MKDIR} -p bin
	${LINK.c} -o bin/BAMscale ${OBJECTFILES} ${LDLIBSOPTIONS} -lBigWig -lhts -lz -lm -lbz2 -llzma -lcurl -ldl -lpthread

${OBJECTDIR}/src/BAMcoverage.o: src/BAMcoverage.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/BAMcoverage.o src/BAMcoverage.c

${OBJECTDIR}/src/BAMstructs.o: src/BAMstructs.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/BAMstructs.o src/BAMstructs.c

${OBJECTDIR}/src/BEDstruct.o: src/BEDstruct.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/BEDstruct.o src/BEDstruct.c

${OBJECTDIR}/src/CHROMstruct.o: src/CHROMstruct.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/CHROMstruct.o src/CHROMstruct.c

${OBJECTDIR}/src/Inputs.o: src/Inputs.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Inputs.o src/Inputs.c

${OBJECTDIR}/src/Writer.o: src/Writer.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Writer.o src/Writer.c

${OBJECTDIR}/src/binning.o: src/binning.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/binning.o src/binning.c

${OBJECTDIR}/src/main.o: src/main.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/main.o src/main.c

${OBJECTDIR}/src/multithreads.o: src/multithreads.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/multithreads.o src/multithreads.c

${OBJECTDIR}/src/scale.o: src/scale.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/scale.o src/scale.c

${OBJECTDIR}/src/segmenter.o: src/segmenter.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/segmenter.o src/segmenter.c

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
