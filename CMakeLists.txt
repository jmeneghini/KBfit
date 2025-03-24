cmake_minimum_required(VERSION 3.10)
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


# External libraries
# OpenMP
find_package(OpenMP REQUIRED)
# LAPACK
find_package(LAPACK REQUIRED)
# MINUIT2
set(MINUIT2_INCDIR /usr/local/include)
set(MINUIT2_LIBDIR /usr/local/lib)
target_include_directories(KBfit PUBLIC ${MINUIT2_INCDIR})

target_link_libraries(KBfit PUBLIC
    OpenMP::OpenMP_CXX
    LAPACK::LAPACK
    ${MINUIT2_LIBDIR}/libMinuit2.so
)

# Put binary in my local bin directory
set(CMAKE_INSTALL_PREFIX /home/jmeneghini/latqcd_miniconda3)

install(TARGETS KBfit 
        RUNTIME DESTINATION bin)
