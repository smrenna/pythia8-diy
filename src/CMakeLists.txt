add_library(diygen SHARED CalcConfig.cpp)
install(TARGETS diygen
    DESTINATION ${CMAKE_INSTALL_PREFIX/lib}
    LIBRARY DESTINATION lib  # <-- Add this line
   COMPONENT library
    )
install(FILES ${MY_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    COMPONENT devel)

find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(YODF5 REQUIRED)
include_directories(${YODF5_INCLUDE_DIR}, ${EIGEN3_INCLUDE_DIR})

add_executable              (pythia8-diy pythia8-diy.cpp)

target_link_libraries (pythia8-diy LINK_PUBLIC ${libraries} diygen YODF5)
install(TARGETS pythia8-diy
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
    )

# add_executable              (pythia8-diy-wcalib pythia8-diy-wcalib.cpp)

# target_link_libraries (pythia8-diy-wcalib LINK_PUBLIC ${libraries} diygen)
# install(TARGETS pythia8-diy-wcalib
#     DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
#     RUNTIME DESTINATION bin
#     )