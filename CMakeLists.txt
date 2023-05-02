cmake_minimum_required(VERSION 3.10)

project(DRX_CPP VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# find_package(glfw REQUIRED)
# find_package(Boost REQUIRED COMPONENTS system)

# add_subdirectory(include/Glad)

set(SOURCE_DIR "src")
set(INCLUDE_DIR "include")

include_directories(include)

add_executable(main "${SOURCE_DIR}/main.cpp")

target_include_directories(main PRIVATE ${INCLUDE_DIR})
# target_link_libraries(main glfw glm::glm)
