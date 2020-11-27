cmake_minimum_required (VERSION 3.14.0 FATAL_ERROR)

project(tbmr LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Set to your opencv dir (containing TBMR in xfeatures2d)
# set(OpenCV_DIR .../opencv/build/)
# add_compile_definitions(-DUSE_OCV_TBMR)

find_package(OpenCV CONFIG REQUIRED)
find_package(Catch2 CONFIG REQUIRED)
find_package(glog CONFIG REQUIRED)
find_package(gflags CONFIG REQUIRED)

# Nasty windows version set
if(WIN32)
  macro(get_WIN32_WINNT version)
    if (WIN32 AND CMAKE_SYSTEM_VERSION)
        set(ver ${CMAKE_SYSTEM_VERSION})
        string(REPLACE "." "" ver ${ver})
        string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})

        set(${version} "0x${ver}")
    endif()
  endmacro()
get_WIN32_WINNT(ver)
add_definitions(-D_WIN32_WINNT=${ver})
endif()

set(TBMR_SRC 
"${CMAKE_CURRENT_SOURCE_DIR}/src/tbmr.hpp"
"${CMAKE_CURRENT_SOURCE_DIR}/src/tbmr.cpp"
"${CMAKE_CURRENT_SOURCE_DIR}/src/component_tree.hpp"
"${CMAKE_CURRENT_SOURCE_DIR}/src/component_tree.cpp"
"${CMAKE_CURRENT_SOURCE_DIR}/src/StopWatch.hpp"
"${CMAKE_CURRENT_SOURCE_DIR}/src/StopWatch.cpp"
"${CMAKE_CURRENT_SOURCE_DIR}/src/Logging.hpp"
)

add_executable(tbmr "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp" ${TBMR_SRC})
target_link_libraries(tbmr ${OpenCV_LIBS} glog::glog)

add_executable(unit_tbmr "${CMAKE_CURRENT_SOURCE_DIR}/src/tests.cpp" ${TBMR_SRC})
target_link_libraries(unit_tbmr ${OpenCV_LIBS} glog::glog Catch2::Catch2)

add_executable(affine_display "${CMAKE_CURRENT_SOURCE_DIR}/src/displayAffineFeatures.cpp")
target_link_libraries(affine_display ${OpenCV_LIBS} glog::glog)