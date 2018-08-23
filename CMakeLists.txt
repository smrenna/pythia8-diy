project                     (gen)
cmake_minimum_required      (VERSION 2.8)

option                      (diy_thread     "Enable diy threading"                 ON)

# C++11
set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

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


# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${DIY_INCLUDE_DIRS}
			    ${YODA_DIR}/include
                            ${HEPMC_DIR}/include
                            ${RIVET_DIR}/include
                            ${PYTHIA8_DIR}/include
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
                             ${CMAKE_THREAD_LIBS_INIT})


#Generate the shared library from the sources
add_library(diygen SHARED config.cpp CalcConfig.cpp)

add_executable              (pythia8-diy pythia8-diy.cpp)
target_link_libraries       (pythia8-diy ${libraries} diygen)


# Test if user has supplied -DHIGHFIVE_DIR                        
if                          (LHEH5_DIR AND HIGHFIVE_DIR)
    add_subdirectory        (withlheh5)
else                        ()
    message                 (ERROR " Need LHEH5 and HighFive in order to build program")
endif                       ()

# Build the reader example
#add_subdirectory        (examples)
