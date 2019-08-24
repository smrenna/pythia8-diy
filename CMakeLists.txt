project                     (pythia8-diy)
cmake_minimum_required      (VERSION 2.8)

option                      (diy_thread     "Enable diy threading"                 OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -ggdb")# -fPermissive")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_STANDARD 11)

# OSX flags
if                          (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions           (-DMAC_OSX)
  set                       (CMAKE_MACOSX_RPATH on)

    # --- following RPATH settings are for Sierra w/ Clang, hopefully they don't hurt other versions
    # ref: https://cmake.org/Wiki/CMake_RPATH_handling
    # use, i.e. don't skip, the full RPATH for the build tree
    set                     (CMAKE_SKIP_BUILD_RPATH            false)
    # when building, don't use the install RPATH already (but later on when installing)
    set                     (CMAKE_BUILD_WITH_INSTALL_RPATH    false)
    # set RPATH to install path
    set                     (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set                     (CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
    # the RPATH to be used when installing, but only if it's not a system directory
    list                    (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
                                                               "${CMAKE_INSTALL_PREFIX}/lib"
                                                               isSystemDir)
    if                      ("${isSystemDir}" STREQUAL         "-1")
      set                   (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    endif                   ()
endif                       (${CMAKE_SYSTEM_NAME} MATCHES      "Darwin")

 # MPI
find_package                (MPI REQUIRED)
    set                     (libraries ${libraries} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

find_package                (HDF5)
    set                     (libraries ${libraries} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})

# Threads
find_package                (Threads)
find_package                (OpenMP)
if                          (OPENMP_FOUND)
  set                       (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set                       (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else                        ()
  message                   ("OpenMP not found")
  add_definitions           (-DMFA_NO_OPENMP)
  # set                       (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif                       ()
if                          (NOT diy_thread)
    message                 ("Diy threading is disabled; setting diy threads will have no effect")
    add_definitions         (-DDIY_NO_THREADS)
endif                       (NOT diy_thread)

add_definitions         (-DDIY_NO_THREADS)
# DIY
find_path                   (DIY_INCLUDE_DIRS diy/types.h)

set (MY_PUBLIC_HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/CalcConfig.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/GenericBlock.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/opts.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Reduce.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Serialisation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Tools.hpp
    )

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${DIY_INCLUDE_DIRS}
			    ${YODA_DIR}/include
                            ${HEPMC_DIR}/include
                            ${RIVET_DIR}/include
                            ${HIGHFIVE_DIR}/include 
                            ${PYTHIA8_DIR}/include
                            ${FASTJET_DIR}/include
                             SYSTEM ${MPI_INCLUDE_PATH}
                            )

# Figure out correct library suffix
if                          (APPLE)
  set                       (CMAKE_SHARED_LIBRARY_SUFFIX ".dylib")
  message                   ("Shared library suffix set to dylib")
else                        ()
  set                       (CMAKE_SHARED_LIBRARY_SUFFIX ".so")
  message                   ("Shared library suffix set to so")
endif                       ()

# Libraries
set                         (libraries
                             ${libraries}
			     ${YODA_DIR}/lib/libYODA${CMAKE_SHARED_LIBRARY_SUFFIX}
			     ${RIVET_DIR}/lib/libRivet${CMAKE_SHARED_LIBRARY_SUFFIX}
                             ${HEPMC_DIR}/lib/libHepMC${CMAKE_SHARED_LIBRARY_SUFFIX}
                             ${PYTHIA8_DIR}/lib/libpythia8${CMAKE_SHARED_LIBRARY_SUFFIX}
                             ${FASTJET_DIR}/lib/libfastjet${CMAKE_SHARED_LIBRARY_SUFFIX}
                             )

# Test if user has supplied -DHIGHFIVE_DIR                        
if (NOT HIGHFIVE_DIR)
    message                 (FATAL_ERROR " Need HighFive in order to build program --- check -DHIGHFIVE_DIR")
endif                       ()

# Test if libs are found
find_library(YODA_LIB
    NAMES libYODA${CMAKE_SHARED_LIBRARY_SUFFIX}
    HINTS "${YODA_DIR}/lib"
    )
if(NOT YODA_LIB)
    message(FATAL_ERROR "YODA library not found --- check -DYODA_DIR")
endif()

find_library(RIVET_LIB
    NAMES libRivet${CMAKE_SHARED_LIBRARY_SUFFIX}
    HINTS "${RIVET_DIR}/lib"
    )
if(NOT RIVET_LIB)
    message(FATAL_ERROR "RIVET library not found --- check -DRIVET_DIR")
endif()

find_library(HEPMC_LIB
    NAMES libHepMC${CMAKE_SHARED_LIBRARY_SUFFIX}
    HINTS "${HEPMC_DIR}/lib"
    )
if(NOT HEPMC_LIB)
    message(FATAL_ERROR "HEPMC library not found --- check -DHEPMC_DIR")
endif()

find_library(PYTHIA8_LIB
    NAMES libpythia8${CMAKE_SHARED_LIBRARY_SUFFIX}
    HINTS "${PYTHIA8_DIR}/lib"
    )
if(NOT PYTHIA8_LIB)
    message(FATAL_ERROR "PYTHIA8 library not found --- check -DPYTHIA8_DIR")
endif()
find_library(FASTJET_LIB
    NAMES libfastjet${CMAKE_SHARED_LIBRARY_SUFFIX}
    HINTS "${FASTJET_DIR}/lib"
    )
if(NOT FASTJET_LIB)
    message(FATAL_ERROR "FASTJET library not found --- check -DFASTJET_DIR")
endif()

# If all found, add build target
add_subdirectory        (src)


## Test if user has supplied -DHIGHFIVE_DIR                        
if                          (LHEH5_DIR AND HIGHFIVE_DIR)
    message("Configuring to build lhd hdf5 support")
    set                         (libraries
        ${libraries}
        ${LHEH5_DIR}/lib/liblheh5${CMAKE_SHARED_LIBRARY_SUFFIX}
        )
    
    find_library(LHEH5_LIB
        NAMES liblheh5${CMAKE_SHARED_LIBRARY_SUFFIX}
        HINTS "${LHEH5_DIR}/lib"
        )
    if(NOT LHEH5_LIB)
        message(FATAL_ERROR "LHEH5 library not found --- check -DLHEH5_DIR")
    endif()

    find_file(HAVE_H5_H
        NAMES H5File.hpp
        HINTS "${HIGHFIVE_DIR}/include/highfive" 
        )
    if(NOT HAVE_H5_H)
        message(FATAL_ERROR "Specified HighFive dir not found --- check -DHIGHFIVE_DIR")
    endif()

    add_subdirectory        (withlheh5)
else                        ()
    message                 (INFO " Need LHEH5 and HighFive in order to build lhehdf program")
endif                       ()

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "8")
set(CPACK_PACKAGE_VERSION_PATCH "2")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES
    "Docker;benchmark;test;.git;~;build;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
