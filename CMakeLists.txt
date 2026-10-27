cmake_minimum_required(VERSION 3.22)

project(Richardson_Gaudin_ED2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Wybieramy OpenBLAS z interfejsem LP64.
set(BLA_VENDOR OpenBLAS)
set(BLA_SIZEOF_INTEGER 4)

find_package(LAPACK REQUIRED)

add_executable(Richardson_Gaudin_ED
        Richardson-Gaudin-ED.cpp
        pairingedlib.cpp
        test.cpp
        input.dat
)

# Podłączenie LAPACK/OpenBLAS do programu.
target_link_libraries(Richardson_Gaudin_ED PRIVATE
        LAPACK::LAPACK
)