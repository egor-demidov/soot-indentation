cmake_minimum_required(VERSION 3.22)
project(soot_indentation)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -march=native -flto=auto ${CMAKE_CXX_FLAGS}")

find_package(TBB REQUIRED)

include_directories(SYSTEM deps/eigen)                  # Eigen, a linear algebra library
include_directories(SYSTEM deps/libtimestep/include)    # module that implements second-order integrators
include_directories(SYSTEM deps/libgran/include)        # module that implements DEM solvers and contact models
include_directories(SYSTEM deps/cgal)                   # Computational Geometry Algorithms Library

set(UTILITY_SOURCES energy/compute_energy.cpp writer/writer.cpp reader/mackowski_reader.cpp)

### VERIFICATION EXECUTABLES ###

add_executable(01_particles_colliding verification/01_particles_colliding.cpp ${UTILITY_SOURCES})
add_executable(02_particles_colliding verification/02_particles_colliding.cpp ${UTILITY_SOURCES})
add_executable(03_particles_colliding verification/03_particles_colliding.cpp ${UTILITY_SOURCES})
add_executable(04_particles_colliding verification/04_particles_colliding.cpp ${UTILITY_SOURCES})
add_executable(05_hamaker verification/05_hamaker.cpp ${UTILITY_SOURCES})
add_executable(06_hamaker verification/06_hamaker.cpp ${UTILITY_SOURCES})
add_executable(07_sintering verification/07_sintering.cpp ${UTILITY_SOURCES})
add_executable(08_sintering verification/08_sintering.cpp ${UTILITY_SOURCES})
add_executable(09_sintering verification/09_sintering.cpp ${UTILITY_SOURCES})
add_executable(10_sintering verification/10_sintering.cpp ${UTILITY_SOURCES})

target_link_libraries(01_particles_colliding TBB::tbb)
target_link_libraries(02_particles_colliding TBB::tbb)
target_link_libraries(03_particles_colliding TBB::tbb)
target_link_libraries(04_particles_colliding TBB::tbb)
target_link_libraries(05_hamaker TBB::tbb)
target_link_libraries(06_hamaker TBB::tbb)
target_link_libraries(07_sintering TBB::tbb)
target_link_libraries(08_sintering TBB::tbb)
target_link_libraries(09_sintering TBB::tbb)
target_link_libraries(10_sintering TBB::tbb)

### PRODUCTION EXECUTABLES ###

add_executable(01_afm_aggregate production/01_afm_aggregate.cpp ${UTILITY_SOURCES})

target_link_libraries(01_afm_aggregate TBB::tbb)
