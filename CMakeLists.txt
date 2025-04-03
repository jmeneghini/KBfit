cmake_minimum_required(VERSION 3.10)

if(APPLE)
    set(LLVM_DIR "/opt/homebrew/opt/llvm")
    set(CMAKE_C_COMPILER "${LLVM_DIR}/bin/clang")
    set(CMAKE_CXX_COMPILER "${LLVM_DIR}/bin/clang++")
else()
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER g++)
endif()

project(KBfit)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Sets the build type to Release (optimized) if it is not set
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Set maximum orbital angular momentum LMAX 
# for AR (at rest), OA (on-axis total momentum (0,0,n)),
# PD (planar diagonal total momentum (0,n,n)),
# CD (cubic diagonal total momentum (n,n,n))
# LMAX<=6 currently supported; decreasing speeds up compile time.
# Setting any one to a negative value eliminates that momentum ray.
add_compile_definitions(AR_LMAX=6)
add_compile_definitions(OA_LMAX=6)
add_compile_definitions(PD_LMAX=6)
add_compile_definitions(CD_LMAX=6)

# Set maximum total intrinsic spin (times two) SX2MAX 
# for AR (at rest), OA (on-axis total momentum (0,0,n)),
# PD (planar diagonal total momentum (0,n,n)),
# CD (cubic diagonal total momentum (n,n,n))
# SX2MAX<=4 currently supported for AR,OA; SX2MAX<=3 for PD,CD.
# Decreasing these values speeds up compile time.
# Setting less than 0 eliminates all routines for the momentum ray.
add_compile_definitions(AR_SX2MAX=4)
add_compile_definitions(OA_SX2MAX=4)
add_compile_definitions(PD_SX2MAX=3)
add_compile_definitions(CD_SX2MAX=3)

# Find the source files
set(KBFIT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SOURCE_DIR ${KBFIT_DIR}/source)

# Collect box matrix sources and exclude them from the main sources
file(GLOB_RECURSE BOX_MATRIX_SOURCES ${SOURCE_DIR}/box_matrix_*.cc)
file(GLOB_RECURSE SOURCES ${SOURCE_DIR}/*.cc)
list(REMOVE_ITEM SOURCES ${BOX_MATRIX_SOURCES})
file(GLOB_RECURSE HEADERS ${SOURCE_DIR}/*.h)

# Create an object library for the box matrix sources (do not optimize)
add_library(BoxMatrix OBJECT ${BOX_MATRIX_SOURCES})
target_compile_options(BoxMatrix PRIVATE -O0)

# Create an object library for the main sources
add_library(OptimizedSources OBJECT ${SOURCES} ${HEADERS})

# Create the executable
add_executable(KBfit
    $<TARGET_OBJECTS:OptimizedSources>
    $<TARGET_OBJECTS:BoxMatrix>
)

# Include the header directories
file(GLOB SUBDIRS LIST_DIRECTORIES true RELATIVE ${SOURCE_DIR} "${SOURCE_DIR}/*")

foreach(subdir IN LISTS SUBDIRS)
    if(IS_DIRECTORY "${SOURCE_DIR}/${subdir}")
        target_include_directories(OptimizedSources PUBLIC "${SOURCE_DIR}/${subdir}")
    endif()
endforeach()

# set compiler to g++
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)

# External libraries (versions on cluster: LAPACK 3.9.0, Minuit2 6.34.04)
# LAPACK
find_package(LAPACK REQUIRED)
# MINUIT2
# if on linux, usr/local/, if on mac, opt/homebrew/
if(APPLE)
    # set(MINUIT2_INCDIR "/opt/homebrew/include")
    # set(MINUIT2_LIBDIR "/opt/homebrew/lib")
    # set(MINUIT2_LIB "libMinuit2.dylib")

    # OpenMP is a major pain in the ass on mac
    # For C
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${LLVM_DIR}/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${LLVM_DIR}/lib/libomp.dylib")

    # For C++
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LLVM_DIR}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${LLVM_DIR}/lib/libomp.dylib")
# else()
#     set(MINUIT2_INCDIR /usr/local/include)
#     set(MINUIT2_LIBDIR /usr/local/lib)
#     set(MINUIT2_LIB "libMinuit2.a")
endif()

find_package(OpenMP REQUIRED)

# Headers
target_include_directories(KBfit PUBLIC ${MINUIT2_INCDIR})

# Libraries
target_link_directories(KBfit PUBLIC ${MINUIT2_LIBDIR})
target_link_libraries(KBfit PUBLIC Minuit2)

# target_include_directories(KBfit PUBLIC ${MINUIT2_INCDIR})

target_link_libraries(KBfit PUBLIC
    OpenMP::OpenMP_CXX
    LAPACK::LAPACK
    # ${MINUIT2_LIBDIR}/${MINUIT2_LIB}
)

# Put binary in my local bin directory
set(CMAKE_INSTALL_PREFIX /home/jmeneghini/latqcd_miniconda3)

install(TARGETS KBfit 
        RUNTIME DESTINATION bin)
