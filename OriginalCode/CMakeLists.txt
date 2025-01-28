cmake_minimum_required(VERSION 3.3)
set(CMAKE_CXX_STANDARD 17)

project(vrp)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif ()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -DNDEBUG")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")


### This Project #####
add_subdirectory(src)

set(SOURCES
        src/main_vrp.cpp
        src/model_data.cpp
        src/probdata_vrp.cpp
        src/pricer_vrp.cpp
        src/label.cpp
        src/labellist.cpp
        src/labeling_algorithm_vrp.cpp
        src/vardata.cpp
        src/initial.cpp
        src/tools_vrp.cpp
        src/branchingrule_arcflow.cpp
        src/ConshdlrArcflow.cpp
        src/branchingrule_dayvar.cpp
        src/ConshdlrDayVar.cpp
        src/branchingrule_vehicle.cpp
        src/ConshdlrVehicle.cpp
        src/ConshdlrNVehicle.cpp
        src/local_search_pricing.cpp
        src/tourVRP.cpp
        src/eventhdlr_nodeInit.cpp
        src/printer.cpp
        src/var_tools.cpp
        src/tsptw.cpp
        src/ConshdlrKPC.cpp
        src/ConshdlrSRC.cpp
        src/heurDayVarRounding.cpp
        src/prop_varfixing.cpp
        src/branchingrule_vehiclearc.cpp
        src/ConshdlrVehicleArc.cpp
)

add_executable(vrp ${SOURCES})

target_include_directories(vrp
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        )

# SCIP
message(STATUS "Looking for SCIP.")
if (NOT SCIP_DIR)
    message(STATUS "No SCIP specified, will look for environment variable.")
    set(SCIP_DIR "$ENV{SCIP_DIR2}" CACHE PATH "Path to SCIP directory")
else ()
    message(STATUS "Will use SCIP_DIR ${SCIP_DIR}.")
endif ()
find_package(SCIP PATHS ${SCIP_DIR} NO_DEFAULT_PATH REQUIRED)
message(STATUS "Using SCIP library from here ${SCIP_DIR}.")


## This project links

target_link_libraries(vrp
        ${SCIP_LIBRARIES}
)




