cmake_minimum_required(VERSION 3.26)
project(soot_indentation)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -march=native -flto=auto ${CMAKE_CXX_FLAGS}")

find_package(TBB REQUIRED)

include_directories(deps/eigen)
include_directories(deps/libtimestep/include)
include_directories(deps/libgran/include)

set(UTILITY_SOURCES energy/compute_energy.cpp writer/writer.cpp)

add_executable(01_particles_colliding verification/01_particles_colliding.cpp ${UTILITY_SOURCES})
add_executable(02_particles_colliding verification/02_particles_colliding.cpp ${UTILITY_SOURCES})
add_executable(03_particles_colliding verification/03_particles_colliding.cpp ${UTILITY_SOURCES})
add_executable(04_particles_colliding verification/04_particles_colliding.cpp ${UTILITY_SOURCES})
add_executable(05_hamaker verification/05_hamaker.cpp ${UTILITY_SOURCES})

target_link_libraries(01_particles_colliding TBB::tbb)
target_link_libraries(02_particles_colliding TBB::tbb)
target_link_libraries(03_particles_colliding TBB::tbb)
target_link_libraries(04_particles_colliding TBB::tbb)
target_link_libraries(05_hamaker TBB::tbb)
