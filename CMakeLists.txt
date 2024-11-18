cmake_minimum_required(VERSION 3.20)
project(sprosim VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create library target with source files
add_library(sprosim
    src/CoffeeBed.cpp
    src/CoffeeParticle.cpp
    src/WaterFlow.cpp
    src/PhysicsSolver.cpp
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