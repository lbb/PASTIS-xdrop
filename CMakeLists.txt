cmake_minimum_required(VERSION 3.11)
project(pastis LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

include_directories("plog/include")
include_directories(
	${PROJECT_SOURCE_DIR}/seqan/include 
	${PROJECT_SOURCE_DIR}/CombBLAS/install/include
	${PROJECT_SOURCE_DIR}/CombBLAS/Applications
    ${PROJECT_SOURCE_DIR}/ipuma-lib/src
    ${PROJECT_SOURCE_DIR}/ipuma-lib
)


link_directories(${PROJECT_SOURCE_DIR}/CombBLAS/install/lib)
    
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DSEQAN_ENABLE_DEBUG=0")
add_compile_definitions(NDEBUG)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

FILE(GLOB PASTIS_CXX_SRCS src/*.cpp src/kmer/*.cpp src/align/*.cpp)
FILE(GLOB PASTIS_CXX_INCS inc/*.hpp inc/kmer/*.hpp inc/align/*.hpp)
add_executable(pastis
	${PASTIS_CXX_SRCS}
	${PASTIS_CXX_INCS}
	)

add_subdirectory(ipuma-lib)

# target_link_libraries(pastis PRIVATE -fsanitize=address)
target_link_libraries(pastis PUBLIC CombBLAS GraphGenlib Usortlib IPUMA_LIBRARY)
if (MPI_FOUND)
	message(STATUS "MPI found.")
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(pastis PUBLIC ${MPI_LIBRARIES})
endif()

include(CheckLanguage)
check_language(CUDA)
if (CUDA_FOUND AND CMAKE_CUDA_COMPILER_LOADED)
    add_compile_definitions(ADEPT_ENABLED)
    find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    enable_language(CUDA)
    include_directories( ${PROJECT_SOURCE_DIR}/ADEPT)
    link_directories(${PROJECT_SOURCE_DIR}/ADEPT/build/adept)
 	message(STATUS "Cuda found...")
     message(STATUS "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
     message(STATUS "${CUDART_LIBRARY}")
     target_include_directories(pastis PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
     target_link_libraries(pastis PUBLIC ${CUDART_LIBRARY})

     # ADEPT
     target_link_libraries(pastis PUBLIC adept_lib_shared)
 else()
     message("NOT BUILDING ADEPT, CUDA MISSING!")
 endif()

if(OPENMP_FOUND OR OpenMP_CXX_FOUND)
    # Set THREADED if OpenMP is found
    target_compile_definitions(pastis PUBLIC THREADED)
    if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(pastis PUBLIC OpenMP::OpenMP_CXX)
    else()
        target_compile_options(pastis PUBLIC "${OpenMP_CXX_FLAGS}")
        target_link_libraries(pastis PUBLIC "${OpenMP_CXX_FLAGS}")
    endif()
endif()
