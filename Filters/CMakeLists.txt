cmake_minimum_required(VERSION 3.5)
project(filters)

#execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOTINCS)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOTINCS} -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -m64 -I/home/jsvirzi/root/include -std=c++11")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -m64 -I/home/jsvirzi/root/include -std=c++11")

# manually throw the ROOT cflags in
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -m64")
include_directories(include)

set(SOURCE_FILES src/filters.cpp)
add_library(filters SHARED ${SOURCE_FILES})

add_executable(lowpass examples/lowpass.c)
target_link_libraries(lowpass filters)
