#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DopplerBroadParallel)

# Parallel processing using TOP-C
option(TOPC_USE "Use TOP-C parallel architecture" OFF)
option(TOPC_USE_SEQ "Use TOP-C in sequential mode" OFF)
if(TOPC_USE)
  if(TOPC_USE_SEQ)
    set(PARMODE "--seq")
  else()
    set(PARMODE "--mpi")
  endif()

  # Get TOP-C libraries
  execute_process(COMMAND "topc-config" ${PARMODE} "--libs"
          OUTPUT_VARIABLE TOPC_LIBS)

  # Remove trailing whitespace
  string(STRIP ${TOPC_LIBS} TOPC_LIBS)

  # Get TOP-C compiler flags
  execute_process(COMMAND "topc-config" ${PARMODE} "--cflags"
          OUTPUT_VARIABLE TOPC_CFLAGS)
endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
           ${PROJECT_SOURCE_DIR}/src/*.icc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
           ${PROJECT_SOURCE_DIR}/include/*.h)


#----------------------------------------------------------------------------
# Add the executable, and link it to the DopplerBroadParallel libraries
#
add_executable(dopplerbroadpara DopplerBroadParrallel.cc ${sources} ${headers})
target_link_libraries(dopplerbroadpara ${TOPC_LIBS})


#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(DopplerBroadParallel DEPENDS dopplerbroadpara)


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS dopplerbroadpara DESTINATION bin)
