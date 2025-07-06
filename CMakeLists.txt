cmake_minimum_required(VERSION 3.20)
project(sprosim VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create library target with source files
add_library(sprosim
    src/CoffeeBed.cpp
    src/CoffeeParticle.cpp
    src/CoffeeParticle3D.cpp
    src/WaterFlow.cpp
    src/PhysicsSolver.cpp
    src/models/transport/LinearExtractionModel.cpp
)

# Set include directories for the library
target_include_directories(sprosim
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Enable testing
include(CTest)

# Install rules for the library
install(TARGETS sprosim
    EXPORT sprosim-targets
    INCLUDES DESTINATION include
)

# Install headers maintaining directory structure
install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

# Add tests subdirectory
add_subdirectory(tests)

# Add visualization subdirectory
add_subdirectory(visualization)

# Enable VTK XML export (no VTK library dependency)
add_compile_definitions(SPROSIM_HAS_VTK)
message(STATUS "VTK XML export enabled - install ParaView to visualize results")

# Add demo executable
add_executable(sprosim_demo
    demo/main.cpp
)



# Link demo with the library
target_link_libraries(sprosim_demo PRIVATE
    sprosim
)



# Link with visualization library
target_link_libraries(sprosim_demo PRIVATE sprosim_visualization)


# Set include directories for demo
target_include_directories(sprosim_demo PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)



# Code formatting with clang-format
find_program(CLANG_FORMAT_EXECUTABLE clang-format)

if(CLANG_FORMAT_EXECUTABLE)
    # Collect all source files
    file(GLOB_RECURSE ALL_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/demo/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/visualization/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/visualization/*.h
    )

    # Format target - formats code in-place
    add_custom_target(format
        COMMAND ${CLANG_FORMAT_EXECUTABLE} -i ${ALL_SOURCE_FILES}
        COMMENT "Formatting source code with clang-format"
        VERBATIM
    )

    # Format-check target - checks formatting without modifying files
    add_custom_target(format-check
        COMMAND ${CLANG_FORMAT_EXECUTABLE} --dry-run --Werror ${ALL_SOURCE_FILES}
        COMMENT "Checking source code formatting with clang-format"
        VERBATIM
    )

    message(STATUS "clang-format found: ${CLANG_FORMAT_EXECUTABLE}")
    message(STATUS "Added targets: format, format-check")
else()
    message(WARNING "clang-format not found. Install clang-format to enable formatting targets.")
endif()
