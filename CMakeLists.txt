cmake_minimum_required(VERSION 3.22)
project(soot_indentation)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -march=native -flto=auto ${CMAKE_CXX_FLAGS}")

find_package(TBB REQUIRED)

set(CGAL_DIR "deps/cgal")
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
find_package(CGAL REQUIRED)

include_directories(SYSTEM deps/eigen)                  # Eigen, a linear algebra library
include_directories(SYSTEM deps/libtimestep/include)    # module that implements second-order integrators
include_directories(SYSTEM deps/libgran/include)        # module that implements DEM solvers and contact models
include_directories(SYSTEM ${CGAL_INCLUDE_DIRS})        # Computational Geometry Algorithms Library
# CGAL dependencies: libmpfr-dev libgmp-dev

set(UTILITY_SOURCES energy/compute_energy.cpp writer/writer.cpp reader/mackowski_reader.cpp afm/anchor_plane.cpp)

set(SHARED_OBJECTS TBB::tbb gmp mpfr)

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

target_link_libraries(01_particles_colliding ${SHARED_OBJECTS})
target_link_libraries(02_particles_colliding ${SHARED_OBJECTS})
target_link_libraries(03_particles_colliding ${SHARED_OBJECTS})
target_link_libraries(04_particles_colliding ${SHARED_OBJECTS})
target_link_libraries(05_hamaker ${SHARED_OBJECTS})
target_link_libraries(06_hamaker ${SHARED_OBJECTS})
target_link_libraries(07_sintering ${SHARED_OBJECTS})
target_link_libraries(08_sintering ${SHARED_OBJECTS})
target_link_libraries(09_sintering ${SHARED_OBJECTS})
target_link_libraries(10_sintering ${SHARED_OBJECTS})

### PRODUCTION EXECUTABLES ###

add_executable(01_afm_aggregate production/01_afm_aggregate.cpp ${UTILITY_SOURCES})
add_executable(02_spiral production/02_spiral.cpp ${UTILITY_SOURCES})
add_executable(03_hooks production/03_hooks.cpp ${UTILITY_SOURCES})

target_link_libraries(01_afm_aggregate ${SHARED_OBJECTS})
target_link_libraries(02_spiral ${SHARED_OBJECTS})
target_link_libraries(03_hooks ${SHARED_OBJECTS})

#add_executable(oscillator_test test/oscillator.cpp)
