cmake_minimum_required(VERSION 3.10)

if(APPLE)
    set(LLVM_DIR "/opt/homebrew/opt/llvm")
    set(CMAKE_C_COMPILER "${LLVM_DIR}/bin/clang")
    set(CMAKE_CXX_COMPILER "${LLVM_DIR}/bin/clang++")
else()
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER g++)
endif()

# Check if lsb_release is available (it is not on the cluster, but is on my desktop)
execute_process(
    COMMAND which lsb_release
    RESULT_VARIABLE LSB_RELEASE_EXISTS
    OUTPUT_QUIET
    ERROR_QUIET
)

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

# ccache for faster compilation
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_PROGRAM})
    message(STATUS "Using ccache for compilation: ${CCACHE_PROGRAM}")
else()
    message(STATUS "ccache not found, not using it.")
endif()

# Set maximum orbital angular momentum LMAX 
# for AR (at rest), OA (on-axis total momentum (0,0,n)),
# PD (planar diagonal total momentum (0,n,n)),
# CD (cubic diagonal total momentum (n,n,n))
# LMAX<=6 currently supported; decreasing speeds up compile time.
# Setting any one to a negative value eliminates that momentum ray.
add_compile_definitions(AR_LMAX=6 OA_LMAX=6 PD_LMAX=6 CD_LMAX=6)
#add_compile_definitions(AR_LMAX=1 OA_LMAX=1 PD_LMAX=1 CD_LMAX=1)

# Set maximum total intrinsic spin (times two) SX2MAX 
# for AR (at rest), OA (on-axis total momentum (0,0,n)),
# PD (planar diagonal total momentum (0,n,n)),
# CD (cubic diagonal total momentum (n,n,n))
# SX2MAX<=4 currently supported for AR,OA; SX2MAX<=3 for PD,CD.
# Decreasing these values speeds up compile time.
# Setting less than 0 eliminates all routines for the momentum ray.
add_compile_definitions(AR_SX2MAX=4 OA_SX2MAX=4 PD_SX2MAX=3 CD_SX2MAX=3)
#add_compile_definitions(AR_SX2MAX=1 OA_SX2MAX=1 PD_SX2MAX=1 CD_SX2MAX=1)


# Find the source files
set(KBFIT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SOURCE_DIR ${KBFIT_DIR}/source)
set(TEST_SOURCE_DIR ${KBFIT_DIR}/source_testing)

# Collect box matrix sources and exclude them from the main sources
file(GLOB_RECURSE BOX_MATRIX_SOURCES ${SOURCE_DIR}/box_matrix_*.cc)
file(GLOB_RECURSE SOURCES ${SOURCE_DIR}/*.cc)
list(REMOVE_ITEM SOURCES ${BOX_MATRIX_SOURCES} ${SOURCE_DIR}/main/KBfit.cc)
file(GLOB_RECURSE HEADERS ${SOURCE_DIR}/*.h)


# External libraries (versions on cluster: LAPACK 3.9.0, Minuit2 6.34.04)
# LAPACK
# set(BLAS_DIR "/usr/lib/x86_64-linux-gnu")
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)
if(HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_LIBRARIES}")
    message(STATUS "HDF5 Include Dirs: ${HDF5_INCLUDE_DIRS}")
    add_compile_definitions(HDF5)
else()
    message(FATAL_ERROR "HDF5 not found - required for KBfit")
endif()

# MINUIT2
# if on linux, usr/local/, if on mac, opt/homebrew/

set(USER_HOME $ENV{HOME})
message(STATUS "User's home folder: ${USER_HOME}")

set(CMAKE_INSTALL_PREFIX ${USER_HOME}/.local)

if(APPLE)
    message(STATUS "On MacOS")
    set(MINUIT2_INCDIR "/opt/homebrew/include/Minuit2")
    set(ROOT_INCDIR "/opt/homebrew/include")
    set(MINUIT2_LIBDIR "/opt/homebrew/lib")
    set(MINUIT2_LIB "libMinuit2.dylib")

    # OpenMP is a major pain in the ass on mac
    # For C
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${LLVM_DIR}/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${LLVM_DIR}/lib/libomp.dylib")

    # For C++
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LLVM_DIR}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${LLVM_DIR}/lib/libomp.dylib")
elseif(LSB_RELEASE_EXISTS EQUAL 0)
    message(STATUS "On Ubuntu Desktop")
    set(MINUIT2_INCDIR "${CMAKE_PREFIX_PATH}/include/Minuit2")
    set(ROOT_INCDIR "")
    set(MINUIT2_LIBDIR "${CMAKE_PREFIX_PATH}/lib")
    set(MINUIT2_LIB "libMinuit2.so")
else()
    message(STATUS "On CMU QCD Cluster")
    set(CMAKE_INSTALL_PREFIX /home/jmeneghini/latqcd_miniconda3)
    set(MINUIT2_INCDIR "/usr/include/root/")
    set(ROOT_INCDIR "")
    set(MINUIT2_LIBDIR "/usr/lib64/root/")
    set(MINUIT2_LIB "libMinuit2.so")
endif()

find_package(OpenMP REQUIRED)


# common headers
set(COMMON_INCLUDE_DIRS " ")
list(APPEND COMMON_INCLUDE_DIRS ${MINUIT2_INCDIR})
if(ROOT_INCDIR AND NOT ROOT_INCDIR STREQUAL MINUIT2_INCDIR)
    list(APPEND COMMON_INCLUDE_DIRS ${ROOT_INCDIR})
endif()
file(GLOB SUBDIRS LIST_DIRECTORIES true RELATIVE ${SOURCE_DIR} "${SOURCE_DIR}/*")
foreach(subdir IN LISTS SUBDIRS)
    if(IS_DIRECTORY "${SOURCE_DIR}/${subdir}")
        list(APPEND COMMON_INCLUDE_DIRS "${SOURCE_DIR}/${subdir}")
    endif()
endforeach()

# Create an object library for the main sources
add_library(OptimizedSources OBJECT ${SOURCES})
target_include_directories(OptimizedSources PUBLIC ${COMMON_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
# if in debug, compile optimized sources as such
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(OptimizedSources PRIVATE -O0 -g)
else()
    target_compile_options(OptimizedSources PRIVATE -O3)
endif()

# Create an static library for the box matrix sources (do not optimize)
add_library(BoxMatrix_lib STATIC ${BOX_MATRIX_SOURCES})
target_compile_options(BoxMatrix_lib PRIVATE -O0)
target_include_directories(BoxMatrix_lib PUBLIC ${COMMON_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

message(STATUS "COMMON_INCLUDE_DIRS: ${COMMON_INCLUDE_DIRS}")

# Create the executable
add_executable(KBfit ${SOURCE_DIR}/main/KBfit.cc
                    $<TARGET_OBJECTS:OptimizedSources>)
target_link_libraries(KBfit PUBLIC
        BoxMatrix_lib
        OpenMP::OpenMP_CXX
        LAPACK::LAPACK
        ${MINUIT2_LIBDIR}/${MINUIT2_LIB}
        ${HDF5_LIBRARIES}
)

# Put binary in my conda bin directory if on AlmaLinux, else usr/local/bin

install(TARGETS KBfit 
        RUNTIME DESTINATION bin)

# testing setup with doctest
set(BUILD_TESTING True)

message(STATUS "BUILD_TESTING: ${BUILD_TESTING}")

if(BUILD_TESTING)
    enable_testing()
    set(DOCTEST_INCLUDE_DIR "${TEST_SOURCE_DIR}/include" CACHE PATH "Path to Doctest include directory")
    message(STATUS "Test source directory: ${TEST_SOURCE_DIR}")
    message(STATUS "Doctest include directory: ${DOCTEST_INCLUDE_DIR}")

    # collect test files
    file(GLOB_RECURSE TEST_SOURCES ${TEST_SOURCE_DIR}/test_*.cc)

    if(TEST_SOURCES STREQUAL "")
        message(WARNING "No test files found in ${TEST_SOURCE_DIR} matching test_*.cc. Test executable will not be built correctly.")
    endif()

    message(STATUS "Test files found: ${TEST_SOURCES}")

    # for cnpy
    find_package(ZLIB REQUIRED)

    add_executable(KBfit_tests
            "${TEST_SOURCE_DIR}/run_tests.cc" # Main test runner
            ${TEST_SOURCES}                     # All test files
           "${TEST_SOURCE_DIR}/cnpy.cc"
    )

    target_include_directories(KBfit_tests PRIVATE
            ${DOCTEST_INCLUDE_DIR}
            ${COMMON_INCLUDE_DIRS}
            ${SOURCE_DIR}
            ${ZLIB_INCLUDE_DIRS}
            ${TEST_SOURCE_DIR}
            ${HDF5_INCLUDE_DIRS}
    )

    target_link_libraries(KBfit_tests PRIVATE
            BoxMatrix_lib                     # Static library being tested
            $<TARGET_OBJECTS:OptimizedSources> # Object files from the other part of your code
            OpenMP::OpenMP_CXX                # Link OpenMP if tests use it or tested code needs its symbols at link time
            LAPACK::LAPACK                    # Link LAPACK for dpotrf_, dsyev_, etc.
            ${ZLIB_LIBRARIES}
            "${MINUIT2_LIBDIR}/${MINUIT2_LIB}" # Link Minuit2 for ROOT::Minuit2:: symbols
            ${HDF5_LIBRARIES}                 # Link HDF5 libraries for tests
    )

    set_target_properties(KBfit_tests PROPERTIES CXX_STANDARD ${CMAKE_CXX_STANDARD})
    set_target_properties(KBfit_tests PROPERTIES CXX_STANDARD_REQUIRED ${CMAKE_CXX_STANDARD_REQUIRED})

    # Define project root path and pass it to the test executable
    target_compile_definitions(KBfit_tests PRIVATE 
        PROJECT_ROOT_DIR="${KBFIT_DIR}"
    )

    # Add the test to CTest
    # run tests with 'ctest' command in the build directory
    add_test(NAME ProjectTests COMMAND KBfit_tests)

    message(STATUS "Testing enabled (Doctest). Test executable: KBfit_tests")
else()
    message(STATUS "Testing disabled.")
endif()