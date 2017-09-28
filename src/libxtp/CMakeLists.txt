add_custom_target(gitversion-xtp COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

file(GLOB_RECURSE VOTCA_SOURCES *.cc *.cpp)
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc )
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

find_package(CUDA)
if (CUDA_FOUND)
  message("-- CUDA found, compiling xtp with gpu acceleration.")
  enable_language(CUDA)
  include_directories(${CUDA_INCLUDE_DIRS})
  link_directories(${CUDA_LIBRARY_DIRS})
  include_directories("/home/pranav/dev/cuda-8.0/samples/common/inc/")
  add_definitions(-DWITH_CUDA)
  add_subdirectory(gpu)
  set(VOTCA_SOURCES ${VOTCA_SOURCES} $<TARGET_OBJECTS:votca_xtp_gpu>)
else()
  message("-- CUDA not found, compiling xtp without gpu acceleration.")
endif()

add_library(votca_xtp  ${VOTCA_SOURCES})
set_target_properties(votca_xtp PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(votca_xtp gitversion-xtp)
target_link_libraries(votca_xtp ${VOTCA_CTP_LIBRARIES} ${GSL_LIBRARIES} ${VOTCA_CSG_LIBRARIES} ${VOTCA_TOOLS_LIBRARIES} ${Boost_LIBRARIES} ${LIBXC_LIBRARIES})
install(TARGETS votca_xtp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

configure_file(libvotca_xtp.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_xtp.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_xtp.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
