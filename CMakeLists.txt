cmake_minimum_required(VERSION 3.3)
project(tsunami Fortran)


list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

# set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# if build type not specified, default to release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "release")
endif()

# handle integer size
if(INT)
  message(STATUS "Configuring build for ${INT}-bit integers")
  add_definitions(-DINT${INT})
else()
  message(STATUS "Configuring build for 32-bit integers")
  add_definitions(-DINT32)
endif()

# handle real size
if(REAL)
  message(STATUS "Configuring build for ${REAL}-bit reals")
  add_definitions(-DREAL${REAL})
else()
  message(STATUS "Configuring build for 32-bit reals")
  add_definitions(-DREAL32)
endif()

# compiler flags for gfortran
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  list(APPEND CMAKE_Fortran_FLAGS -cpp)
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -C -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

# compiler flags for ifort
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  list(APPEND CMAKE_Fortran_FLAGS -coarray -assume realloc_lhs -heap-arrays)
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -C -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

# compiler flags for CRAY ftn
if(CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
  list(APPEND CMAKE_Fortran_FLAGS -h noomp)
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

find_package(Coarray REQUIRED)

add_library(tsunami src/lib/mod_boundary.f90 src/lib/mod_diagnostics.f90 src/lib/mod_diff.f90 src/lib/mod_field.f90 src/lib/mod_io.f90 src/lib/mod_kinds.f90 src/lib/mod_parallel.f90)
target_compile_options(tsunami PUBLIC ${Coarray_COMPILE_OPTIONS})
target_link_libraries(tsunami PUBLIC ${Coarray_LIBRARIES})

add_executable(app src/app/tsunami.f90)
add_executable(app_dt src/app/tsunami_dt.f90)
target_link_libraries(app tsunami)
target_link_libraries(app_dt tsunami)

enable_testing()
foreach(exec field sync_edges sync_on_assign tiling)
  add_executable(test_${exec} src/tests/test_${exec}.f90)
  target_link_libraries(test_${exec} tsunami)
  add_test(NAME ${exec} 
    COMMAND ${Coarray_EXECUTABLE} ${Coarray_NUMPROC_FLAG} ${Coarray_MAX_NUMPROCS} bin/test_${exec})
endforeach()
