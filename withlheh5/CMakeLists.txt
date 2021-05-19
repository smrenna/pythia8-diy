include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/../include
                            ${DIY_INCLUDE_DIRS}
                            ${YODA_DIR}/include
                            ${HEPMC_DIR}/include
                            ${RIVET_DIR}/include
                            ${PYTHIA8_DIR}/include
#                            ${LHEH5_DIR}/include
                            ${HIGHFIVE_DIR}/include
                            SYSTEM ${MPI_INCLUDE_PATH}
                            )


#target_link_libraries       (withlheh5 ${libraries} diygen ${LHEH5_DIR}/lib/liblheh5${CMAKE_SHARED_LIBRARY_SUFFIX})

add_executable              (withlheh5 withlheh5.cpp)
target_link_libraries (withlheh5 LINK_PUBLIC ${libraries} diygen)
install(TARGETS withlheh5
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
    )
