cmake_minimum_required(VERSION 3.2)

# Adding customized cmake module
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(angpow)

    # Uses system libraries or downloads and build if necessary
    include(BuildFFTW)
    include_directories(${FFTW_INCLUDES})
    link_directories(${FFTW_INCLUDE_DIRS})

    # Compilation flags
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "^(Apple)?Clang$")
      # using Clang, disabling OpenMP support
      set(CMAKE_CXX_FLAGS "-fno-common  -O3 -fPIC -Wa,-q -Wall -Wpedantic -std=c++11 -march=native -ffast-math -m64")
    else()
      set(CMAKE_CXX_FLAGS "-Wall -Wpointer-arith -O3 -fPIC -std=c++11 -fopenmp  -fno-common  -Wall -march=native -ffast-math -m64")
    endif()


    # Defines list of Angpow library src files
    set(ANGPOW_SRC src/walltimer.cc src/walltime_c.cc src/angpow_bessel.cc
	                 src/angpow_chebyshevInt.cc src/angpow_radint.cc src/angpow_kinteg.cc
                   src/angpow_pk2cl.cc src/angpow_tools.cc src/angpow_clbase.cc
                   src/angpow_powspec.cc src/angpow_parameters.cc src/angpow_cosmo.cc
                   src/angpow_ccl.cc src/angpow_ctheta.cc)


    # Adds path to include folder
    include_directories(inc)

    # Builds the main Angpow library
    add_library(libangpow ${ANGPOW_SRC})
    target_link_libraries(libangpow ${FFTW_LIBRARIES} m)
    if(NOT FFTW_FOUND)
      add_dependencies(libangpow FFTW)
    endif()

    add_executable(angpow angpow.cc)
    target_link_libraries(angpow libangpow)
    #
    # add_executable(limber limber.cc)
    # target_link_libraries(limber libangpow)
