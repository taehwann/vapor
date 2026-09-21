cmake_minimum_required(VERSION 3.20)
project(vapor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(VAPOR_BUILD_APP "Build the interactive OpenGL application" ON)
include(CTest)


# Source groups, not libraries: all production code is compiled into vapor.
# Tests reuse these lists directly, without a vapor_core target.
set(VAPOR_NUMERICAL_SOURCES
    src/common/ProjectionOptions.hpp
    src/common/DomainInfo.hpp
    src/solvers/mac3d/MacGridDomain3D.hpp
    src/solvers/mac3d/MacGridDomain3D.cpp
    src/solvers/mac3d/MacGridState3D.hpp
    src/solvers/mac3d/MacGridState3D.cpp
    src/solvers/mac3d/MacGridOperators3D.hpp
    src/solvers/mac3d/MacGridOperators3D.cpp
    src/solvers/mac3d/MacGridSampler3D.hpp
    src/solvers/mac3d/MacGridSampler3D.cpp
    src/solvers/mac3d/VelocityFieldView.hpp
    src/solvers/mac3d/AdvectionContext.hpp
    src/solvers/mac3d/AdvectionContext.cpp
    src/solvers/mac3d/CpuAdvectionKernels.hpp
    src/solvers/mac3d/CpuAdvectionKernels.cpp
    src/solvers/mac3d/Advection3D.hpp
    src/solvers/mac3d/CpuAdvection.hpp
    src/solvers/mac3d/CpuAdvection.cpp
    src/solvers/mac3d/Projection3D.hpp
    src/solvers/mac3d/CpuPressureSolver3D.hpp
    src/solvers/mac3d/CpuPressureSolver3D.cpp
    src/common/Vec3.hpp
    src/common/LinearSolveOptions.hpp
    src/solvers/mac3d/StencilLinearSystem.hpp
    src/solvers/mac3d/StencilLinearSystem.cpp
    src/solvers/mac3d/SOR.hpp
    src/solvers/mac3d/SOR.cpp)

set(VAPOR_SIMULATION_SOURCES
    src/common/FluidDiagnostics.hpp
    src/solvers/mac3d/MacGridFluidSolver3D.hpp
    src/solvers/mac3d/MacGridFluidSolver3D.cpp)

# Presentation reads simulation state; simulation sources never depend on it.
# These adapters require no graphics context and are also tested headlessly.
set(VAPOR_PRESENTATION_SOURCES
    src/renderer/MacRenderData.hpp
    src/renderer/MacRenderData.cpp
    src/renderer/RenderData.hpp)

set(VAPOR_2D_SOURCES
    src/solvers/mac2d/ReflectionMacFluidSolver2D.hpp
    src/solvers/mac2d/MacGridDomain2D.hpp src/solvers/mac2d/MacGridDomain2D.cpp
    src/solvers/mac2d/MacGridState2D.hpp src/solvers/mac2d/MacGridState2D.cpp
    src/solvers/mac2d/MacGridOperators2D.hpp src/solvers/mac2d/MacGridOperators2D.cpp
    src/solvers/mac2d/SemiLagrangian2D.hpp src/solvers/mac2d/SemiLagrangian2D.cpp
    src/solvers/mac2d/CpuPressureSolver2D.hpp src/solvers/mac2d/CpuPressureSolver2D.cpp
    src/solvers/mac2d/MacGridFluidSolver2D.hpp src/solvers/mac2d/MacGridFluidSolver2D.cpp)
set(VAPOR_SIMPLICIAL_2D_SOURCES
    src/solvers/simplicial2d/SimplicialMesh2D.hpp src/solvers/simplicial2d/SimplicialMesh2D.cpp
    src/solvers/simplicial2d/SimplicialFluidState2D.hpp src/solvers/simplicial2d/SimplicialFluidState2D.cpp
    src/solvers/simplicial2d/SimplicialOperators2D.hpp src/solvers/simplicial2d/SimplicialOperators2D.cpp
    src/solvers/simplicial2d/SimplicialFluidSolver2D.hpp src/solvers/simplicial2d/SimplicialFluidSolver2D.cpp)
set(VAPOR_SIMPLICIAL_3D_SOURCES
    src/solvers/simplicial3d/SimplicialMesh3D.hpp src/solvers/simplicial3d/SimplicialMesh3D.cpp
    src/solvers/simplicial3d/SimplicialDual3D.hpp src/solvers/simplicial3d/SimplicialDual3D.cpp
    src/solvers/simplicial3d/SimplicialBoundary3D.hpp src/solvers/simplicial3d/SimplicialBoundary3D.cpp
    src/solvers/simplicial3d/SimplicialFluidSolver3D.hpp src/solvers/simplicial3d/SimplicialFluidSolver3D.cpp)
configure_file(examples/simplicial3d/box.tet examples/simplicial3d/box.tet COPYONLY)
configure_file(examples/simplicial3d-bunny/bunny.tet examples/simplicial3d-bunny/bunny.tet COPYONLY)
configure_file(examples/simplicial3d-bunny/bunny-fluid.tet examples/simplicial3d-bunny/bunny-fluid.tet COPYONLY)
set(VAPOR_2D_PRESENTATION_SOURCES
    src/renderer/MacRenderData2D.hpp
    src/renderer/MacRenderData2D.cpp
    src/renderer/ImageRenderData.hpp)
set(VAPOR_SIMPLICIAL_2D_PRESENTATION_SOURCES
    src/solvers/simplicial2d/SimplicialRasterizer2D.hpp
    src/solvers/simplicial2d/SimplicialRasterizer2D.cpp)

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
    vapor_read_shader(divSrc src/solvers/mac3d/gpu/shaders/divergence.comp)
    vapor_read_shader(rbgsSrc src/solvers/mac3d/gpu/shaders/sor.comp)
    vapor_read_shader(correctSrc src/solvers/mac3d/gpu/shaders/correct_velocity.comp)
    vapor_read_shader(advectSlSrc src/solvers/mac3d/gpu/shaders/semi_lagrangian.comp)
    vapor_read_shader(advectMcSrc src/solvers/mac3d/gpu/shaders/maccormack.comp)
    vapor_read_shader(copySrc src/solvers/mac3d/gpu/shaders/copy.comp)
    vapor_read_shader(residentOps src/solvers/mac3d/gpu/shaders/resident_ops.comp)
    vapor_read_shader(vert src/renderer/shaders/volume.vert)
    vapor_read_shader(frag src/renderer/shaders/volume.frag)
    vapor_read_shader(imageVert src/renderer/shaders/image.vert)
    vapor_read_shader(imageFrag src/renderer/shaders/image.frag)
    vapor_read_shader(simplicialAdvection src/solvers/simplicial3d/gpu/shaders/advection.comp)
    vapor_read_shader(simplicialRecovery src/solvers/simplicial3d/gpu/shaders/recovery.comp)
    set_source_files_properties(${VAPOR_SHADER_FILES} PROPERTIES HEADER_FILE_ONLY TRUE)
    configure_file(src/graphics/ShaderSources.hpp.in generated/ShaderSources.hpp @ONLY)

    set(VAPOR_GPU_SOURCES
        src/solvers/mac3d/gpu/GpuMacSimulation3D.hpp
        src/solvers/mac3d/gpu/GpuMacSimulation3D.cpp
        src/solvers/mac3d/gpu/GpuMacBackend.hpp
        src/solvers/mac3d/gpu/GpuMacBackend.cpp
        src/solvers/mac3d/gpu/GpuPressureSolver.hpp
        src/solvers/mac3d/gpu/GpuPressureSolver.cpp
        src/solvers/mac3d/gpu/GpuAdvection.hpp
        src/solvers/mac3d/gpu/GpuAdvection.cpp
        src/graphics/GlLoader.hpp
        src/graphics/ShaderSources.hpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/generated/ShaderSources.hpp
        ${VAPOR_SHADER_FILES})
    set(VAPOR_APP_SUPPORT_SOURCES
        src/app/WindowContext.hpp src/app/WindowContext.cpp)
    set(VAPOR_RENDERER_SOURCES
        src/renderer/TetrahedralWireData.hpp
        src/renderer/TetrahedralWireRenderer.hpp
        src/renderer/TetrahedralWireRenderer.cpp
        src/renderer/Renderer2D.hpp
        src/renderer/Renderer2D.cpp
        src/renderer/VolumeRenderer.hpp
        src/renderer/VolumeRenderer.cpp
        src/renderer/OrbitCamera.hpp
        src/renderer/OrbitCamera.cpp)

    set(VAPOR_APPLICATION_SOURCES

        src/app/SimulationApp.hpp
        src/app/SimulationApp.cpp
        src/app/SimulationApp2D.cpp
        src/app/SimulationApp3D.cpp
        src/app/SimulationAppSimplicial3D.cpp
        src/app/SimulationAppMesh3D.cpp
        src/renderer/SimplicialRenderData3D.hpp
        ${VAPOR_SIMPLICIAL_3D_SOURCES}

        ${VAPOR_APP_SUPPORT_SOURCES}
        ${VAPOR_2D_SOURCES}
        ${VAPOR_2D_PRESENTATION_SOURCES}
        ${VAPOR_SIMPLICIAL_2D_SOURCES}
        ${VAPOR_SIMPLICIAL_2D_PRESENTATION_SOURCES}
        ${VAPOR_NUMERICAL_SOURCES}
        ${VAPOR_SIMULATION_SOURCES}
        ${VAPOR_PRESENTATION_SOURCES}
        ${VAPOR_GPU_SOURCES}
        ${VAPOR_RENDERER_SOURCES})
    add_executable(vapor main.cpp ${VAPOR_APPLICATION_SOURCES})
    target_include_directories(vapor PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
    target_link_libraries(vapor PRIVATE imgui_impl)
    add_custom_command(TARGET vapor POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:vapor>/examples/simplicial3d"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/simplicial3d/box.tet" "$<TARGET_FILE_DIR:vapor>/examples/simplicial3d/box.tet")
    set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT vapor)
    add_custom_command(TARGET vapor POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/simplicial3d-bunny"
            "$<TARGET_FILE_DIR:vapor>/examples/simplicial3d-bunny")
endif()

if(BUILD_TESTING)
    add_executable(vapor_teapot_2d_tests tests/SimplicialTeapot2DTests.cpp
        ${VAPOR_SIMPLICIAL_2D_SOURCES} ${VAPOR_SIMPLICIAL_2D_PRESENTATION_SOURCES})
    target_include_directories(vapor_teapot_2d_tests PRIVATE src)
    foreach(case IN ITEMS geometry walls recovery stability)
        add_test(NAME teapot2d_${case} COMMAND vapor_teapot_2d_tests ${case})
    endforeach()
    set_tests_properties(teapot2d_stability PROPERTIES TIMEOUT 180)
    if(VAPOR_BUILD_APP)
        add_executable(vapor_teapot_2d_graphics_tests tests/SimplicialTeapot2DGraphicsTests.cpp
            ${VAPOR_SIMPLICIAL_2D_SOURCES} ${VAPOR_SIMPLICIAL_2D_PRESENTATION_SOURCES}
            src/renderer/TetrahedralWireRenderer.cpp src/renderer/Renderer2D.cpp)
        target_include_directories(vapor_teapot_2d_graphics_tests PRIVATE src ${CMAKE_CURRENT_BINARY_DIR}/generated)
        target_link_libraries(vapor_teapot_2d_graphics_tests PRIVATE glfw OpenGL::GL)
        add_test(NAME graphics_teapot2d COMMAND vapor_teapot_2d_graphics_tests)
        set_tests_properties(graphics_teapot2d PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 90)
        add_executable(vapor_startup_tests tests/StartupTests.cpp ${VAPOR_APPLICATION_SOURCES})
        target_include_directories(vapor_startup_tests PRIVATE src ${CMAKE_CURRENT_BINARY_DIR}/generated)
        target_link_libraries(vapor_startup_tests PRIVATE imgui_impl)
        add_test(NAME application_solver_selection COMMAND vapor_startup_tests)
        set_tests_properties(application_solver_selection PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 120)
    endif()

    add_executable(vapor_tetrahedral_wire_tests tests/TetrahedralWireTests.cpp)
    target_include_directories(vapor_tetrahedral_wire_tests PRIVATE src)
    add_test(NAME tetrahedral_wire_geometry COMMAND vapor_tetrahedral_wire_tests)
    add_executable(vapor_simplicial_geometry_tests tests/SimplicialGeometry3DTests.cpp
        src/solvers/simplicial3d/SimplicialMesh3D.cpp src/solvers/simplicial3d/SimplicialDual3D.cpp)
    target_include_directories(vapor_simplicial_geometry_tests PRIVATE src)
    add_test(NAME simplicial3d_geometry COMMAND vapor_simplicial_geometry_tests)
    add_executable(vapor_simplicial_3d_tests tests/Simplicial3DTests.cpp ${VAPOR_SIMPLICIAL_3D_SOURCES})
    target_include_directories(vapor_simplicial_3d_tests PRIVATE src)
    add_executable(vapor_bunny_3d_tests tests/SimplicialBunny3DTests.cpp ${VAPOR_SIMPLICIAL_3D_SOURCES})
    target_include_directories(vapor_bunny_3d_tests PRIVATE src)
    add_test(NAME simplicial3d_bunny COMMAND vapor_bunny_3d_tests)
    set_tests_properties(simplicial3d_bunny PROPERTIES TIMEOUT 600)
    add_executable(vapor_simplicial_3d_performance tests/Simplicial3DPerformance.cpp ${VAPOR_SIMPLICIAL_3D_SOURCES})
    target_include_directories(vapor_simplicial_3d_performance PRIVATE src)
    foreach(case IN ITEMS topology metric recovery circulation forces boundaries simulation contracts stability transport)
        add_test(NAME simplicial3d_${case} COMMAND vapor_simplicial_3d_tests ${case})
    endforeach()
    set_tests_properties(simplicial3d_stability PROPERTIES TIMEOUT 120)
    add_executable(vapor_historical_2d_tests tests/Historical2DTests.cpp
        tests/fixtures/Ba20987Smoke2D.hpp src/solvers/mac2d/ReflectionMacFluidSolver2D.hpp)
    add_test(NAME mac2d_historical_parity COMMAND vapor_historical_2d_tests)
    add_executable(vapor_numerics_tests tests/NumericsTests.cpp ${VAPOR_NUMERICAL_SOURCES})
    foreach(case IN ITEMS
            grid_layout grid_validation closed_boundaries
            sor_manufactured sor_options sor_validation
            projection_rest projection_gradient projection_divergence
            projection_solenoidal projection_system projection_legacy)
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
    add_test(NAME simulation_seeded_replay COMMAND vapor_simulation_tests seeded_replay)
    add_test(NAME simulation_plume COMMAND vapor_simulation_tests plume)

    add_executable(vapor_backend_tests tests/BackendTests.cpp
        ${VAPOR_NUMERICAL_SOURCES} ${VAPOR_SIMULATION_SOURCES}
        ${VAPOR_PRESENTATION_SOURCES})
    add_test(NAME solver_backends COMMAND vapor_backend_tests)


    add_executable(vapor_construction_tests tests/ConstructionTests.cpp
        ${VAPOR_SIMPLICIAL_3D_SOURCES}
 ${VAPOR_NUMERICAL_SOURCES}
        ${VAPOR_SIMULATION_SOURCES} ${VAPOR_2D_SOURCES}
        ${VAPOR_SIMPLICIAL_2D_SOURCES}
 ${VAPOR_PRESENTATION_SOURCES}
        ${VAPOR_2D_PRESENTATION_SOURCES} ${VAPOR_SIMPLICIAL_2D_PRESENTATION_SOURCES})
    add_test(NAME simulation_construction COMMAND vapor_construction_tests)

    add_executable(vapor_2d_tests tests/Mac2DTests.cpp
        ${VAPOR_2D_SOURCES} ${VAPOR_2D_PRESENTATION_SOURCES})
    foreach(case IN ITEMS mac2d_grid mac2d_advection mac2d_projection_gradient
            mac2d_projection_solenoidal mac2d_projection_options mac2d_simulation mac2d_reflection mac2d_presentation)
        add_test(NAME ${case} COMMAND vapor_2d_tests ${case})
    endforeach()

    add_executable(vapor_simplicial_2d_tests tests/Simplicial2DTests.cpp
        ${VAPOR_SIMPLICIAL_2D_SOURCES} ${VAPOR_SIMPLICIAL_2D_PRESENTATION_SOURCES})
    foreach(case IN ITEMS simplicial2d_mesh simplicial2d_velocity simplicial2d_recovery
                          simplicial2d_circulation simplicial2d_simulation simplicial2d_cfl
                          simplicial2d_stability simplicial2d_presentation simplicial2d_walls simplicial2d_refinement)
        add_test(NAME ${case} COMMAND vapor_simplicial_2d_tests ${case})
    endforeach()

    if(VAPOR_BUILD_APP)
        add_executable(vapor_tetrahedral_wire_graphics_tests tests/TetrahedralWireGraphicsTests.cpp
            src/renderer/TetrahedralWireRenderer.cpp src/renderer/OrbitCamera.cpp)
        target_include_directories(vapor_tetrahedral_wire_graphics_tests PRIVATE src)
        target_link_libraries(vapor_tetrahedral_wire_graphics_tests PRIVATE glfw OpenGL::GL)
        add_test(NAME graphics_tetrahedral_wire COMMAND vapor_tetrahedral_wire_graphics_tests)
        set_tests_properties(graphics_tetrahedral_wire PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 30)
        add_executable(vapor_simplicial_3d_graphics_tests tests/Simplicial3DGraphicsTests.cpp
            ${VAPOR_SIMPLICIAL_3D_SOURCES} ${VAPOR_RENDERER_SOURCES})
        target_include_directories(vapor_simplicial_3d_graphics_tests PRIVATE src ${CMAKE_CURRENT_BINARY_DIR}/generated)
        target_link_libraries(vapor_simplicial_3d_graphics_tests PRIVATE glfw OpenGL::GL)
        add_test(NAME graphics_simplicial3d COMMAND vapor_simplicial_3d_graphics_tests)
        set_tests_properties(graphics_simplicial3d PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 90)
        add_executable(vapor_bunny_3d_graphics_tests tests/SimplicialBunny3DGraphicsTests.cpp
            ${VAPOR_SIMPLICIAL_3D_SOURCES} ${VAPOR_RENDERER_SOURCES})
        target_include_directories(vapor_bunny_3d_graphics_tests PRIVATE src ${CMAKE_CURRENT_BINARY_DIR}/generated)
        target_link_libraries(vapor_bunny_3d_graphics_tests PRIVATE glfw OpenGL::GL)
        add_test(NAME graphics_bunny3d COMMAND vapor_bunny_3d_graphics_tests)
        set_tests_properties(graphics_bunny3d PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600)
        add_executable(vapor_gpu_performance tests/GpuPerformance.cpp
            ${VAPOR_NUMERICAL_SOURCES} ${VAPOR_SIMULATION_SOURCES}
            ${VAPOR_PRESENTATION_SOURCES} ${VAPOR_GPU_SOURCES} ${VAPOR_RENDERER_SOURCES})
        target_include_directories(vapor_gpu_performance PRIVATE src ${CMAKE_CURRENT_BINARY_DIR}/generated)
        target_link_libraries(vapor_gpu_performance PRIVATE glfw OpenGL::GL)
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
        vapor_simulation_tests vapor_backend_tests vapor_construction_tests vapor_graphics_tests
        vapor_2d_tests vapor_historical_2d_tests vapor_simplicial_2d_tests vapor_graphics_2d_tests)
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

find_package(Threads REQUIRED)
if(BUILD_TESTING AND VAPOR_BUILD_APP)
    add_executable(vapor_simplicial_gpu_tests tests/SimplicialGpu3DTests.cpp
        ${VAPOR_SIMPLICIAL_3D_SOURCES} src/solvers/simplicial3d/gpu/GpuSimplicial3D.cpp)
    target_include_directories(vapor_simplicial_gpu_tests PRIVATE src ${CMAKE_CURRENT_BINARY_DIR}/generated)
    target_link_libraries(vapor_simplicial_gpu_tests PRIVATE glfw OpenGL::GL Threads::Threads)
    add_test(NAME simplicial3d_gpu COMMAND vapor_simplicial_gpu_tests)
    set_tests_properties(simplicial3d_gpu PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} SKIP_RETURN_CODE 77 TIMEOUT 600)
endif()
foreach(vapor_threaded_target IN ITEMS vapor vapor_startup_tests vapor_construction_tests
        vapor_simplicial_3d_tests vapor_bunny_3d_tests vapor_simplicial_3d_performance
        vapor_simplicial_3d_graphics_tests vapor_bunny_3d_graphics_tests)
    if(TARGET ${vapor_threaded_target})
        target_link_libraries(${vapor_threaded_target} PRIVATE Threads::Threads)
    endif()
endforeach()
