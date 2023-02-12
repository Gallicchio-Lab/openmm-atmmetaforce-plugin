#---------------------------------------------------
# OpenMM ATMMetaForce Plugin
#----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.17)

PROJECT(OpenMMATMMetaForcePlugin)

# We need to know where OpenMM is installed so we can access the headers and libraries.
SET(OPENMM_DIR "/usr/local/openmm" CACHE PATH "Where OpenMM is installed")
INCLUDE_DIRECTORIES("${OPENMM_DIR}/include")
LINK_DIRECTORIES("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

# Specify the C++ version we are building for.
SET (CMAKE_CXX_STANDARD 11)

# Set flags for linking on mac
IF(APPLE)
    SET (CMAKE_INSTALL_NAME_DIR "@rpath")
    SET(EXTRA_COMPILE_FLAGS "-msse2 -stdlib=libc++")
ENDIF(APPLE)

# Select where to install
IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
    IF(WIN32)
        SET(CMAKE_INSTALL_PREFIX "$ENV{ProgramFiles}/OpenMM" CACHE PATH "Where to install the plugin" FORCE)
    ELSE(WIN32)
        SET(CMAKE_INSTALL_PREFIX "/usr/local/openmm" CACHE PATH "Where to install the plugin" FORCE)
    ENDIF(WIN32)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
#SET(ATMMETAFORCE_PLUGIN_SOURCE_SUBDIRS openmmapi serialization)
SET(ATMMETAFORCE_PLUGIN_SOURCE_SUBDIRS openmmapi)

# Set the library name
SET(ATMMETAFORCE_LIBRARY_NAME ATMMetaForcePlugin)
SET(SHARED_ATMMETAFORCE_TARGET ${ATMMETAFORCE_LIBRARY_NAME})

# These are all the places to search for header files which are to be part of the API.
SET(API_INCLUDE_DIRS "openmmapi/include" "openmmapi/include/internal")

# Locate header files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
ENDFOREACH(dir)

# Collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)
FOREACH(subdir ${ATMMETAFORCE_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# Create the library.

ADD_LIBRARY(${SHARED_ATMMETAFORCE_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_ATMMETAFORCE_TARGET}
    PROPERTIES COMPILE_FLAGS "-DATMMETAFORCE_BUILDING_SHARED_LIBRARY ${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_COMPILE_FLAGS}")
TARGET_LINK_LIBRARIES(${SHARED_ATMMETAFORCE_TARGET} OpenMM)
INSTALL(TARGETS ${SHARED_ATMMETAFORCE_TARGET} DESTINATION lib)

# install headers
FILE(GLOB API_ONLY_INCLUDE_FILES "openmmapi/include/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES} DESTINATION include)
FILE(GLOB API_ONLY_INCLUDE_FILES_INTERNAL "openmmapi/include/internal/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES_INTERNAL} DESTINATION include/internal)

# Enable testing

ENABLE_TESTING()

# Build the implementations for different platforms

ADD_SUBDIRECTORY(platforms/reference)
ADD_SUBDIRECTORY(platforms/common)

FIND_PACKAGE(OpenCL QUIET)
IF(OpenCL_FOUND)
    SET(ATMMETAFORCE_BUILD_OPENCL_LIB ON CACHE BOOL "Build implementation for OpenCL")
ELSE(OpenCL_FOUND)
    SET(ATMMETAFORCE_BUILD_OPENCL_LIB OFF CACHE BOOL "Build implementation for OpenCL")
ENDIF(OpenCL_FOUND)
IF(ATMMETAFORCE_BUILD_OPENCL_LIB)
    ADD_SUBDIRECTORY(platforms/opencl)
ENDIF(ATMMETAFORCE_BUILD_OPENCL_LIB)

FIND_PACKAGE(CUDA QUIET)
IF(CUDA_FOUND)
    SET(ATMMETAFORCE_BUILD_CUDA_LIB ON CACHE BOOL "Build implementation for CUDA")
ELSE(CUDA_FOUND)
    SET(ATMMETAFORCE_BUILD_CUDA_LIB OFF CACHE BOOL "Build implementation for CUDA")
ENDIF(CUDA_FOUND)
IF(ATMMETAFORCE_BUILD_CUDA_LIB)
    ADD_SUBDIRECTORY(platforms/cuda)
ENDIF(ATMMETAFORCE_BUILD_CUDA_LIB)

#builds for the HIP platform if the plugin is installed
IF ( EXISTS ${OPENMM_DIR}/lib/plugins/libOpenMMHIP.so )
   SET(ATMMETAFORCE_BUILD_HIP_LIB ON CACHE BOOL "Build implementation for HIP")
   ADD_SUBDIRECTORY(platforms/hip)
ELSE()
   SET(ATMMETAFORCE_BUILD_HIP_LIB OFF CACHE BOOL "Build implementation for HIP")
ENDIF( )

# Build the Python API

FIND_PROGRAM(PYTHON_EXECUTABLE python)
FIND_PROGRAM(SWIG_EXECUTABLE swig)
IF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(ATMMETAFORCE_BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build wrappers for Python")
ELSE(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(ATMMETAFORCE_BUILD_PYTHON_WRAPPERS OFF CACHE BOOL "Build wrappers for Python")
ENDIF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
IF(ATMMETAFORCE_BUILD_PYTHON_WRAPPERS)
    ADD_SUBDIRECTORY(python)
ENDIF(ATMMETAFORCE_BUILD_PYTHON_WRAPPERS)
