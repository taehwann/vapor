cmake_minimum_required(VERSION 3.20)
project(vapor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(dependency/glfw)

set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dependency/imgui)

add_library(imgui_impl
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp)
target_include_directories(imgui_impl PUBLIC ${IMGUI_DIR} ${IMGUI_DIR}/backends)
target_link_libraries(imgui_impl PUBLIC glfw)

find_package(OpenGL REQUIRED)
target_link_libraries(imgui_impl PUBLIC OpenGL::GL)

add_executable(vapor main.cpp gl_loader.h)
target_link_libraries(vapor PRIVATE imgui_impl)
