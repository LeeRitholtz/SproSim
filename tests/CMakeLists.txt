# File: tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.20)

# Find Catch2
find_package(Catch2 3 REQUIRED)

# Add test executable
add_executable(sprosim_tests
    physics/PhysicsSolverTests.cpp
    physics/CoffeeParticle3DTests.cpp
    permeability/test_KozenyCarmanPermeabilityModel.cpp
    permeability/test_ConstantPermeabilityModel.cpp
    flow/test_DarcyFlowModel.cpp
    transport/test_LinearExtractionModel.cpp
)

# Include directories
target_include_directories(sprosim_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link libraries
target_link_libraries(sprosim_tests PRIVATE
    sprosim
    Catch2::Catch2WithMain
)

# Add the test
add_test(NAME sprosim_tests COMMAND sprosim_tests)
