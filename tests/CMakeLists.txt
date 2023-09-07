cmake_minimum_required(VERSION 3.11.3 )

project(signal_decomposition_tester)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

############################################################################################
#Fetch Googletest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
##############################################################################################

enable_testing()

#add_executable(exe_name main.cpp additional.cpp)
add_executable(
  signal_decomposition_test
  signal_decomposition_test.cc
  ../src/wave_gen.cpp
)

#target_include_directories(exe_name PUBLIC "Path/To/Headers")
target_include_directories(signal_decomposition_test PUBLIC "../include/")

#target_link_libraries(exe_name Libraries)
target_link_libraries(
  signal_decomposition_test
  GTest::gtest_main
  Eigen3::Eigen
)

include(GoogleTest)
gtest_discover_tests(signal_decomposition_test)
