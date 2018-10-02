
add_library(diygen SHARED config.cpp CalcConfig.cpp)
install(TARGETS diygen
    DESTINATION ${CMAKE_INSTALL_PREFIX/lib}
    LIBRARY DESTINATION lib  # <-- Add this line
   COMPONENT library
    )
install(FILES ${MY_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    COMPONENT devel)


add_executable              (pythia8-diy pythia8-diy.cpp)

target_link_libraries (pythia8-diy LINK_PUBLIC ${libraries} diygen )
install(TARGETS pythia8-diy
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
    )

add_executable              (mg5_py8_diy mg5_py8_diy.cpp)
target_link_libraries (mg5_py8_diy LINK_PUBLIC ${libraries} diygen)
install(TARGETS mg5_py8_diy
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
	)
