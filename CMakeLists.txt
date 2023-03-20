cmake_minimum_required(VERSION 3.2)

# Adding customized cmake module
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(angpow)

    # Uses system libraries or downloads and build if necessary
    include(FindFFTW)
    include(FindGSL)
    include_directories(${FFTW_INCLUDES} ${GSL_INCLUDE_DIRS})
    link_directories(${GSL_LIBRARY_DIRS})

    # Compilation flags
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "^(Apple)?Clang$")
      # using Clang, disabling OpenMP support
      set(CMAKE_CXX_FLAGS "-fno-common  -O3 -fPIC -Wall -Wpedantic -std=c++11 -march=native -ffast-math -m64")
    else()
      set(CMAKE_CXX_FLAGS "-Wall -Wpointer-arith -O3 -fPIC -std=c++11 -fopenmp  -fno-common  -Wall -march=native -ffast-math -m64")
    endif()

    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

    # Defines list of Angpow library src files
    set(ANGPOW_SRC src/walltimer.cc src/walltime_c.cc src/angpow_bessel.cc
	                 src/angpow_chebyshevInt.cc src/angpow_radint.cc src/angpow_kinteg.cc
                   src/angpow_pk2cl.cc src/angpow_tools.cc src/angpow_clbase.cc
                   src/angpow_powspec.cc src/angpow_parameters.cc src/angpow_cosmo.cc
                   src/angpow_ctheta.cc src/angpow_ccl.cc)


    # Adds path to include folder
    include_directories(inc)

    # Builds the main Angpow library
    add_library(libangpow ${ANGPOW_SRC})
    set_target_properties(libangpow PROPERTIES PREFIX "")
    target_link_libraries(libangpow ${FFTW_LIBRARIES} ${GSL_LIBRARIES} m)

    # Compiles angpow executable
    add_executable(angpow angpow.cc)
    target_link_libraries(angpow libangpow)

    # Creates a pkg-config file pour faire pro
    SET(PROJECT_DESCRIPTION "AngPow: Angular Power spectrum")
    SET(PKG_CONFIG_LIBDIR "${CMAKE_INSTALL_PREFIX}/lib")
    SET(PKG_CONFIG_INCLUDEDIR "${CMAKE_INSTALL_PREFIX}/include/")
    SET(PKG_CONFIG_REQUIRES "gsl >= 1.6")
    SET(PKG_CONFIG_REQUIRES_PRIVATE "fftw3")
    SET(PKG_CONFIG_LIBS "-L${PKG_CONFIG_LIBDIR} -langpow -lm")
    SET(PKG_CONFIG_CFLAGS "-I${PKG_CONFIG_INCLUDEDIR}")
    CONFIGURE_FILE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/pkg-config.cmakein"
                    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc" )

    # Installs library and include files
    install(TARGETS libangpow DESTINATION lib)
    install(TARGETS angpow DESTINATION bin)
    install(DIRECTORY inc/Angpow DESTINATION include
            FILES_MATCHING PATTERN "*.h")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/angpow.pc DESTINATION lib/pkgconfig)
