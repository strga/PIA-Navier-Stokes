cmake_minimum_required(VERSION 3.26)
project(Lid_Driven_Cavity_Navier_Stokes)

set(CMAKE_CXX_STANDARD 17)

# Find OpenMP package
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    # Include OpenMP flags for both compiling and linking
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib ${OpenMP_CXX_FLAGS}")
endif()

# Add your executable
add_executable(Lid_Driven_Cavity_Navier_Stokes main.cpp Simulation.cpp write_to_file.cpp)

# Link OpenMP to your executable
if(OpenMP_CXX_FOUND)
    target_link_libraries(Lid_Driven_Cavity_Navier_Stokes PUBLIC OpenMP::OpenMP_CXX)
    # Explicitly link the OpenMP library
    target_link_libraries(Lid_Driven_Cavity_Navier_Stokes PUBLIC "-lomp")
endif()
