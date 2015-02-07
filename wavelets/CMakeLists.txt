
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")
set(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}")
#set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}")

project(WAVELETS C CXX)
set(WAVELETS_MAJOR_VERSION 1)
set(WAVELETS_MINOR_VERSION 4)
set(WAVELETS_PATCH_VERSION 2)
set(WAVELETS_VERSION "${WAVELETS_MAJOR_VERSION}.${WAVELETS_MINOR_VERSION}.${WAVELETS_PATCH_VERSION}")

find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)
find_package(Eigen3 REQUIRED)

#Compiler options
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message("-- Compiling for ${CMAKE_SYSTEM_NAME} with Clang...")
    set(CXX_FLAGS "-std=c++11 -m64")
    set(CXX_WERR "-Wall -Wextra -Wmissing-format-attribute -Wmissing-noreturn -Wredundant-decls -Wsequence-point -Wswitch-default -Wdeprecated -Wunreachable-code  -Wno-sign-conversion -Wno-old-style-cast -Wcovered-switch-default -Wmissing-variable-declarations -Wfloat-equal -Wunknown-warning-option")
    set(CXX_WNOERR "-Wno-weak-vtables -Wno-c++98-compat-pedantic -Wno-unused-parameter -Wno-deprecated-register -Wno-conversion -Wno-shadow -Wno-padded -Wno-global-constructors -Wno-exit-time-destructors -Wno-source-uses-openmp -Wno-effc++")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("-- Compiling for ${CMAKE_SYSTEM_NAME} with GCC...")
    set(CXX_FLAGS "-std=c++11 -m64")
    set(CXX_WERR "-Wall -Wextra -Wmissing-format-attribute -Wmissing-noreturn -Wredundant-decls -Wsequence-point -Wdeprecated -Wunreachable-code -Wno-old-style-cast -Wfloat-equal -Wsuggest-attribute=const -Wsuggest-attribute=pure -Wno-sign-conversion -Wdouble-promotion")
    set(CXX_WNOERR "-Wno-unused-parameter -Wno-conversion -Wno-shadow -Wno-padded -Wno-effc++")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message("-- Compiling for ${CMAKE_SYSTEM_NAME} with Intel C++...")
    message(WARNING "No custom compiler flags set for Intel C++...")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message("-- Compiling for ${CMAKE_SYSTEM_NAME} with Visual Studio C++...")
    message(WARNING "No custom compiler flags set for Visual Studio C++...")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_FLAGS} ${CXX_WERR} ${CXX_WNOERR} ${MPI_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} ")

#Compiler defines
set(COMPILER_DEFINES "-D_REENTRANT")
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug") 
    set(COMPILER_DEFINES "${COMPILER_DEFINES} -D_DEBUG")
    set(BUILD_TYPE "DEBUG")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release") 
    set(COMPILER_DEFINES "${COMPILER_DEFINES} -D_RELEASE")
    set(BUILD_TYPE "RELEASE")
else()
    set(COMPILER_DEFINES "${COMPILER_DEFINES} -D_NORMAL")
    set(BUILD_TYPE "NORMAL")
endif()
add_definitions("${COMPILER_DEFINES}")

# Better Stronger Faster Harder
# Work  it harder
# Make  it better
# Do    it faster
# Makes us stronger
set(BUILD_SHARED_LIBS OFF)

#Set executable sources
list(APPEND executable_sources 
    "${CMAKE_SOURCE_DIR}/src/main.cpp"
    "${CMAKE_SOURCE_DIR}/src/bench.cpp"
    )

#Find source files
file(GLOB_RECURSE header_files src/*.hpp)
file(GLOB_RECURSE source_files src/*.cpp)
file(GLOB_RECURSE template_files src/*.tpp)
file(GLOB_RECURSE configure_files src/*.hpp.in)
file(GLOB_RECURSE test_sources tests/*.cpp)

#Remove executable sources so that mains do not conflict
set(source_files_no_main ${source_files})
foreach(main_source ${executable_sources})
    LIST(REMOVE_ITEM source_files_no_main "${main_source}")
endforeach()

#Configure headers
foreach(configIn ${configure_files})
    string(REPLACE ".hpp.in" ".hpp.out" configOut "${configIn}")
    configure_file("${configIn}" "${configOut}" @ONLY)
    list(APPEND header_files "${configOut}")
endforeach()

#Local includes
list(APPEND LOCAL_INCLUDE_DIRS "src/") #to avoid empty list
foreach(header ${header_files})
    get_filename_component(dir ${header} DIRECTORY)
    list(APPEND LOCAL_INCLUDE_DIRS ${dir})
endforeach()
list(REMOVE_DUPLICATES LOCAL_INCLUDE_DIRS)

#Header directories
include_directories(${LOCAL_INCLUDE_DIRS})
include_directories(SYSTEM 
    ${Boost_INCLUDEDIR}   #Boost
    ${EIGEN3_INCLUDE_DIR} #Eigen
    )

#Library directories
link_directories(
    ${BOOST_LIBRARYDIR} #Boost
    )

# Create static library (all code excluding mains and Cuda precompiled code)
set(CORE_LIB "Core")
add_library(
    "${CORE_LIB}" 
    STATIC 
    ${source_files_no_main}
    ${template_files}
    ${header_files}
)

#Executables
foreach(main_source ${executable_sources})
    get_filename_component(exec_name ${main_source} NAME_WE)
    add_executable(${exec_name} ${header_files} ${main_source} ${template_files})
    target_link_libraries(
        ${exec_name} 
        ${CORE_LIB}          #Core
        ${Boost_LIBRARIES}   #Boost
    )
endforeach()

#Tests
foreach(main_source ${test_sources})
    get_filename_component(exec_name ${main_source} NAME_WE)
    add_executable(${exec_name} ${header_files} ${main_source} ${template_files})
    target_link_libraries(
        ${exec_name} 
        ${CORE_LIB}          #Core
        ${Boost_LIBRARIES}   #Boost
        )
endforeach()

#enable_testing()

