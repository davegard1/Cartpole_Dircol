cmake_minimum_required (VERSION 3.10)

# set the project name
project(cartpole VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)


# Need to add cmake rules/module to find adolc
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(PSOPT CONFIG REQUIRED ) 
find_package(Eigen3 REQUIRED )
find_package(adolc REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)






# Create executable
add_executable( cartpole cartpole.cxx )

add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt adolc Eigen3::Eigen)


