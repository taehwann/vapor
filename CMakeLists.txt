cmake_minimum_required(VERSION 3.20)
project(vapor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(VAPOR_BUILD_APP "Build the interactive OpenGL application" ON)
include(CTest)

set(VAPOR_CONFIG_SOURCES
    src/config/SimulationConfig.hpp
    src/config/SimulationConfig.cpp)
configure_file(vapor.cfg vapor.cfg COPYONLY)

# Source groups, not libraries: all production code is compiled into vapor.
# Tests reuse these lists directly, without a vapor_core target.
set(VAPOR_NUMERICAL_SOURCES
    src/domain/IDomain.hpp
    src/domain/MacGridDomain.hpp
    src/domain/MacGridDomain.cpp
    src/fluid/MacGridState.hpp
    src/fluid/MacGridState.cpp
    src/fluid/MacGridOperators.hpp
    src/fluid/MacGridOperators.cpp
    src/fluid/MacGridSampler.hpp
    src/fluid/MacGridSampler.cpp
    src/fluid/VelocityFieldView.hpp
    src/fluid/advection/AdvectionContext.hpp
    src/fluid/advection/AdvectionContext.cpp
    src/fluid/advection/CpuAdvectionKernels.hpp
    src/fluid/advection/CpuAdvectionKernels.cpp
    src/fluid/advection/IAdvection.hpp
    src/fluid/advection/MacGridAdvection.hpp
    src/fluid/advection/SemiLagrangian.hpp
    src/fluid/advection/SemiLagrangian.cpp
    src/fluid/advection/MacCormack.hpp
    src/fluid/advection/MacCormack.cpp
    src/fluid/projection/IPressureSolver.hpp
    src/fluid/projection/CpuPressureSolver.hpp
    src/fluid/projection/CpuPressureSolver.cpp
    src/math/Vec3.hpp
    src/numerics/ILinearSolver.hpp
    src/numerics/LinearSolveOptions.hpp
    src/numerics/StencilLinearSystem.hpp
    src/numerics/StencilLinearSystem.cpp
    src/numerics/SOR.hpp
    src/numerics/SOR.cpp)

set(VAPOR_SIMULATION_SOURCES
    src/fluid/IFluidSolver.hpp
    src/fluid/MacGridFluidSolver.hpp
    src/fluid/MacGridFluidSolver.cpp)

# Presentation reads simulation state; simulation sources never depend on it.
# These adapters require no graphics context and are also tested headlessly.
set(VAPOR_PRESENTATION_SOURCES
    src/presentation/IPresentationSource.hpp
    src/presentation/MacGridPresentationAdapter.hpp
    src/presentation/MacGridPresentationAdapter.cpp
    src/renderer/RenderData.hpp)

set(VAPOR_2D_SOURCES
    src/domain/MacGridDomain2D.hpp src/domain/MacGridDomain2D.cpp
    src/fluid/MacGridState2D.hpp src/fluid/MacGridState2D.cpp
    src/fluid/MacGridOperators2D.hpp src/fluid/MacGridOperators2D.cpp
    src/fluid/advection/SemiLagrangian2D.hpp src/fluid/advection/SemiLagrangian2D.cpp
    src/fluid/projection/CpuPressureSolver2D.hpp src/fluid/projection/CpuPressureSolver2D.cpp
    src/fluid/MacGridFluidSolver2D.hpp src/fluid/MacGridFluidSolver2D.cpp)
set(VAPOR_COMPOSITION_SOURCES
    src/composition/SimulationInstance.hpp
    src/composition/SimulationInstance.cpp
    src/composition/SimulationRegistry.hpp
    src/composition/SimulationRegistry.cpp
    src/composition/MacSimulationFactory.hpp
    src/composition/MacSimulationFactory.cpp
    src/composition/MacSimulationRegistry.hpp
    src/composition/MacSimulationRegistry.cpp)
set(VAPOR_2D_PRESENTATION_SOURCES
    src/presentation/MacGridPresentationAdapter2D.hpp
    src/presentation/MacGridPresentationAdapter2D.cpp
    src/renderer/ImageRenderData.hpp)

if(VAPOR_BUILD_APP)
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
    find_package(OpenGL REQUIRED)
    target_link_libraries(imgui_impl PUBLIC glfw OpenGL::GL)

    # GLSL remains editable on disk and is embedded for runtime compilation.
    macro(vapor_read_shader symbol path)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/${path}" ${symbol})
        list(APPEND VAPOR_SHADER_FILES "${path}")
    endmacro()
    vapor_read_shader(divSrc src/gpu/shaders/divergence.comp)
    vapor_read_shader(rbgsSrc src/gpu/shaders/sor.comp)
    vapor_read_shader(correctSrc src/gpu/shaders/correct_velocity.comp)
    vapor_read_shader(advectSlSrc src/gpu/shaders/semi_lagrangian.comp)
    vapor_read_shader(advectMcSrc src/gpu/shaders/maccormack.comp)
    vapor_read_shader(copySrc src/gpu/shaders/copy.comp)
    vapor_read_shader(vert src/renderer/shaders/volume.vert)
    vapor_read_shader(frag src/renderer/shaders/volume.frag)
    vapor_read_shader(imageVert src/renderer/shaders/image.vert)
    vapor_read_shader(imageFrag src/renderer/shaders/image.frag)
    set_source_files_properties(${VAPOR_SHADER_FILES} PROPERTIES HEADER_FILE_ONLY TRUE)
    configure_file(src/graphics/ShaderSources.hpp.in generated/ShaderSources.hpp @ONLY)

    set(VAPOR_GPU_SOURCES
        src/gpu/GpuMacBackend.hpp
        src/gpu/GpuMacBackend.cpp
        src/gpu/GpuPressureSolver.hpp
        src/gpu/GpuPressureSolver.cpp
        src/gpu/GpuAdvection.hpp
        src/gpu/GpuAdvection.cpp
        src/graphics/GlLoader.hpp
        src/graphics/ShaderSources.hpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/generated/ShaderSources.hpp
        ${VAPOR_SHADER_FILES})
    set(VAPOR_GPU_REGISTRATION_SOURCES
        src/composition/GpuSimulationRegistration.hpp
        src/composition/GpuMacModuleRegistration.hpp
        src/composition/GpuMacModuleRegistration.cpp)
    set(VAPOR_RENDERER_SOURCES
        src/renderer/Renderer2D.hpp
        src/renderer/Renderer2D.cpp
        src/renderer/Renderer.hpp
        src/renderer/Renderer.cpp
        src/renderer/OrbitCamera.hpp
        src/renderer/OrbitCamera.cpp)

    add_executable(vapor main.cpp
        vapor.cfg
        ${VAPOR_CONFIG_SOURCES}
        src/app/SimulationApp.hpp
        src/app/SimulationApp.cpp
        src/app/SimulationApp3D.cpp
        ${VAPOR_COMPOSITION_SOURCES}
        ${VAPOR_GPU_REGISTRATION_SOURCES}
        ${VAPOR_2D_SOURCES}
        ${VAPOR_2D_PRESENTATION_SOURCES}
        ${VAPOR_NUMERICAL_SOURCES}
        ${VAPOR_SIMULATION_SOURCES}
        ${VAPOR_PRESENTATION_SOURCES}
        ${VAPOR_GPU_SOURCES}
        ${VAPOR_RENDERER_SOURCES})
    target_include_directories(vapor PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
    target_link_libraries(vapor PRIVATE imgui_impl)
    add_custom_command(TARGET vapor POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/vapor.cfg" "$<TARGET_FILE_DIR:vapor>/vapor.cfg")
    set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT vapor)
endif()

if(BUILD_TESTING)
    add_executable(vapor_numerics_tests tests/NumericsTests.cpp ${VAPOR_NUMERICAL_SOURCES})
    foreach(case IN ITEMS
            grid_layout grid_validation closed_boundaries
            sor_manufactured sor_options sor_validation
            projection_rest projection_gradient projection_divergence
            projection_solenoidal projection_injection projection_legacy)
        add_test(NAME ${case} COMMAND vapor_numerics_tests ${case})
    endforeach()

    add_executable(vapor_advection_tests tests/AdvectionTests.cpp ${VAPOR_NUMERICAL_SOURCES})
    foreach(case IN ITEMS sampling advection_stationary advection_transport
                          advection_limiter advection_contracts)
        add_test(NAME ${case} COMMAND vapor_advection_tests ${case})
    endforeach()

    # The full coordinator can now be tested without a graphics context or libs.
    add_executable(vapor_simulation_tests tests/SimulationTests.cpp
        ${VAPOR_NUMERICAL_SOURCES} ${VAPOR_SIMULATION_SOURCES})
    add_test(NAME simulation_legacy_steps COMMAND vapor_simulation_tests)

    add_executable(vapor_interface_tests tests/InterfaceTests.cpp
        ${VAPOR_NUMERICAL_SOURCES} ${VAPOR_SIMULATION_SOURCES}
        ${VAPOR_PRESENTATION_SOURCES})
    add_test(NAME solver_interfaces COMMAND vapor_interface_tests)

    add_executable(vapor_config_tests tests/ConfigTests.cpp ${VAPOR_CONFIG_SOURCES})
    add_test(NAME simulation_config COMMAND vapor_config_tests)

    add_executable(vapor_composition_tests tests/CompositionTests.cpp
        ${VAPOR_CONFIG_SOURCES} ${VAPOR_NUMERICAL_SOURCES}
        ${VAPOR_SIMULATION_SOURCES} ${VAPOR_2D_SOURCES}
        ${VAPOR_COMPOSITION_SOURCES} ${VAPOR_PRESENTATION_SOURCES}
        ${VAPOR_2D_PRESENTATION_SOURCES})
    add_test(NAME simulation_composition COMMAND vapor_composition_tests)

    add_executable(vapor_2d_tests tests/Mac2DTests.cpp
        ${VAPOR_2D_SOURCES} ${VAPOR_2D_PRESENTATION_SOURCES})
    foreach(case IN ITEMS mac2d_grid mac2d_advection mac2d_projection_gradient
            mac2d_projection_solenoidal mac2d_projection_options mac2d_simulation mac2d_interfaces)
        add_test(NAME ${case} COMMAND vapor_2d_tests ${case})
    endforeach()

    if(VAPOR_BUILD_APP)
        add_executable(vapor_graphics_tests tests/GraphicsSmokeTests.cpp
            ${VAPOR_NUMERICAL_SOURCES} ${VAPOR_SIMULATION_SOURCES}
            ${VAPOR_PRESENTATION_SOURCES}
            ${VAPOR_GPU_SOURCES} ${VAPOR_RENDERER_SOURCES})
        target_include_directories(vapor_graphics_tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
        target_link_libraries(vapor_graphics_tests PRIVATE glfw OpenGL::GL)
        add_test(NAME graphics_smoke COMMAND vapor_graphics_tests)
        set_tests_properties(graphics_smoke PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 30)

        add_executable(vapor_graphics_2d_tests tests/Graphics2DTests.cpp
            ${VAPOR_2D_SOURCES} ${VAPOR_2D_PRESENTATION_SOURCES}
            ${VAPOR_RENDERER_SOURCES} src/graphics/GlLoader.hpp
            ${CMAKE_CURRENT_BINARY_DIR}/generated/ShaderSources.hpp ${VAPOR_SHADER_FILES})
        target_include_directories(vapor_graphics_2d_tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
        target_link_libraries(vapor_graphics_2d_tests PRIVATE glfw OpenGL::GL)
        add_test(NAME graphics_2d COMMAND vapor_graphics_2d_tests)
        set_tests_properties(graphics_2d PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 30)
    endif()
endif()

foreach(vapor_target IN ITEMS vapor vapor_numerics_tests vapor_advection_tests
        vapor_simulation_tests vapor_interface_tests vapor_config_tests vapor_composition_tests vapor_graphics_tests
        vapor_2d_tests vapor_graphics_2d_tests)
    if(TARGET ${vapor_target})
        target_include_directories(${vapor_target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
        get_target_property(vapor_target_sources ${vapor_target} SOURCES)
        list(REMOVE_ITEM vapor_target_sources "${CMAKE_CURRENT_BINARY_DIR}/generated/ShaderSources.hpp")
        source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${vapor_target_sources})
    endif()
endforeach()
if(VAPOR_BUILD_APP)
    source_group("Generated" FILES "${CMAKE_CURRENT_BINARY_DIR}/generated/ShaderSources.hpp")
endif()
