project(hmin322-tp-idc-3D)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_PYTHON            "Use OpenGL"         OFF)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
file(GLOB SRCFILES tp/*.cpp)
add_executable(tp_bin ${SRCFILES})
target_link_libraries(tp_bin_bin igl::core igl::opengl_glfw)