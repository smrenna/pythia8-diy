include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/../
                            ${DIY_INCLUDE_DIRS}
			    ${YODA_DIR}/include
                            ${HEPMC_DIR}/include
                            ${RIVET_DIR}/include
                            ${PYTHIA8_DIR}/include
                            ${LHEH5_DIR}/include
                            ${HIGHFIVE_DIR}/include
                            SYSTEM ${MPI_INCLUDE_PATH}
                            )


add_executable              (withlheh5 withlheh5.cpp)
target_link_libraries       (withlheh5 ${libraries} diygen ${LHEH5_DIR}/lib/liblheh5${CMAKE_SHARED_LIBRARY_SUFFIX})
