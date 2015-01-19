
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")
set(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}")
#set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}")

project(Boids C CXX)
set(BOIDS_MAJOR_VERSION 0)
set(BOIDS_MINOR_VERSION 1)
set(BOIDS_PATCH_VERSION 0)
set(BOIDS_VERSION ${BOIDS_MAJOR_VERSION}.${BOIDS_MINOR_VERSION}.${BOIDS_PATCH_VERSION})
get_filename_component(SHADER_FOLDER "${CMAKE_SOURCE_DIR}/shaders" ABSOLUTE)

# TODO # TODO # TODO # TODO # TODO # TODO ###
set(MIN_COMPUTE_CAPABILITY "400")          #T#
set(COMPUTE_CAPABILITY "12")              #O#
set(ENABLE_GUI "1")                       ###
set(CONSOLE_LOG_LEVEL "INFO")             #D#
set(FILE_LOG_LEVEL "DEBUG")               #O#
set(COUT_LOG_LEVEL "5")                   ###
# TODO # TODO # TODO # TODO # TODO # TODO ###

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(GLUT REQUIRED)
find_package(Qt4 REQUIRED QtCore QtGui QtXml QtOpenGL)
find_package(QGLViewer REQUIRED)
find_package(Log4cpp REQUIRED)
if(!{OPENGL_GLU_FOUND})
    message(FATAL_ERROR "GLU not found !")
else()
    message("-- Found OpenGL libraries including GLU !")
endif()

find_package(MPI REQUIRED)
find_package(OpenMP)
find_package(CUDA 5.5)
if(${CUDA_FOUND})
    include(CUDA_compute_capability)
else()
    set(CUDA_TOOLKIT_ROOT_DIR "")
    set(CUDA_NVCC_EXECUTABLE "")
    set(CUDA_INCLUDE_DIRS "")
    set(CUDA_CUDART_LIBRARY "")
    set(CUDA_LIBRARIES "")
    set(CUDA_LIBRARY_DIRS "")
endif()

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
    set(COMPILER_DEFINES "${COMPILER_DEFINES} -D_DEBUG -D_CONSOLE_LOG_LEVEL=DEBUG -D_FILE_LOG_LEVEL=DEBUG")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release") 
    set(COMPILER_DEFINES "${COMPILER_DEFINES} -D_CONSOLE_LOG_LEVEL=WARN -D_FILE_LOG_LEVEL=WARN")
else()
    set(COMPILER_DEFINES "${COMPILER_DEFINES} -D_CONSOLE_LOG_LEVEL=INFO -D_FILE_LOG_LEVEL=INFO")
endif()
add_definitions("${COMPILER_DEFINES}")

# Better Stronger Faster Harder
# Work  it harder
# Make  it better
# Do    it faster
# Makes us stronger
set(BUILD_SHARED_LIBS OFF)

if(${CUDA_FOUND})
    set(CUDA_SEPARABLE_COMPILATION ON)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    if (${CUDA_COMPUTE_CAPABILITY} STREQUAL "21")
        set(CUDA_COMPUTE_CAPABILITY "20")
    endif()
    set(COMPUTE_CAPABILITY ${CUDA_COMPUTE_CAPABILITY})
    #set(CUDA_NVCC_FLAGS "-arch=compute_${CUDA_COMPUTE_CAPABILITY};-code=sm_${CUDA_COMPUTE_CAPABILITY};-std=c++11")
    set(CUDA_NVCC_FLAGS "-arch=compute_${CUDA_COMPUTE_CAPABILITY};-code=sm_${CUDA_COMPUTE_CAPABILITY}")
    set(CUDA_NVCC_FLAGS_DEBUG "-G;-g;-Xptxas=-v")
    set(CUDA_NVCC_FLAGS_RELEASE "-O3")
endif()

#Set executable sources
list(APPEND executable_sources "${CMAKE_SOURCE_DIR}/src/main.cpp")
list(APPEND executable_sources "${CMAKE_SOURCE_DIR}/src/serial.cpp")
list(APPEND executable_sources "${CMAKE_SOURCE_DIR}/src/tree.cpp")
if(${ENABLE_GUI})
    list(APPEND executable_sources "${CMAKE_SOURCE_DIR}/src/display.cpp")
endif()

#Find source files
file(GLOB_RECURSE header_files src/*.hpp)
file(GLOB_RECURSE source_files src/*.cpp)
file(GLOB_RECURSE template_files src/*.tpp)
file(GLOB_RECURSE kernel_files src/*.cu)
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
    ${OPENGL_INCLUDE_DIR}
    ${GLEW_INCLUDE_DIRS} 
    ${GLUT_INCLUDE_DIR}
    ${QGLVIEWER_INCLUDE_DIR}
    ${CUDA_INCLUDE_DIRS} 
    ${LOG4CPP_INCLUDE_DIRS} 
    ${MPI_INCLUDE_PATH}  
    ${QT_INCLUDE_DIR}
    ${QT_QTCORE_INCLUDE_DIR}
    ${QT_QTGUI_INCLUDE_DIR}
    ${QT_QTXML_INCLUDE_DIR}
    ${QT_QTOPENGL_INCLUDE_DIR}
    )

#Library directories
link_directories(
    ${CUDA_LIBRARY_DIRS}
    ${LOG4CPP_LIBRARY_DIRS}
    )

#Compile QT as SYSTEM
set_property(TARGET Qt4::QtCore APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES 
    ${QT_INCLUDE_DIR}
    ${QT_QTCORE_INCLUDE_DIR}
    ${QT_MKSPECS_DIR}/default
    )
set_property(TARGET Qt4::QtGui APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES 
    ${QT_QTGUI_INCLUDE_DIR}
    )
set_property(TARGET Qt4::QtXml APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES 
    ${QT_QTXML_INCLUDE_DIR}
    )
set_property(TARGET Qt4::QtOpenGL APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES 
    ${QT_QTOPENGL_INCLUDE_DIR}
    )

#Precompiled CUDA kernels (static library)
if(${CUDA_FOUND} AND NOT "${kernel_files}" STREQUAL "")
    set(CUDA_KERNELS "CudaKernels")
    cuda_include_directories(${include_directories})
    if(${CMAKE_VERSION} VERSION_GREATER 3.0.0)
        cuda_add_library(${CUDA_KERNELS} ${header_files} ${kernel_files})
    else()
        message(WARNING "Trying to work around a CMake 2.8 bug with CUDA flags")
        cuda_add_library(${CUDA_KERNELS} ${header_files} ${kernel_files} OPTIONS -arch sm_20)
    endif()
endif()

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
        ${CUDA_KERNELS}      #Cuda Kernels
        Qt4::QtGui           #QtGui
        Qt4::QtCore          #QtCore
        Qt4::QtXml           #QtXml
        Qt4::QtOpenGL        #QtOpenGL
        ${QGLVIEWER_LIBRARY} #QGlViewer
        ${LOG4CPP_LIBRARIES} #Log4cpp
        ${MPI_CXX_LIBRARIES} #MPI
        ${CUDA_LIBRARIES}    #Cuda
        ${OPENGL_LIBRARIES}  #OpenGl
        ${GLEW_LIBRARIES}    #Glew
        ${GLUT_glut_LIBRARY} #Glut
        )
endforeach()

#Tests
foreach(main_source ${test_sources})
    get_filename_component(exec_name ${main_source} NAME_WE)
    add_executable(${exec_name} ${header_files} ${main_source} ${template_files})
    target_link_libraries(
        ${exec_name} 
        ${CORE_LIB}          #Core
        ${CUDA_KERNELS}      #Cuda Kernels
        Qt4::QtGui           #QtGui
        Qt4::QtCore          #QtCore
        Qt4::QtXml           #QtXml
        Qt4::QtOpenGL        #QtOpenGL
        ${QGLVIEWER_LIBRARY} #QGlViewer
        ${LOG4CPP_LIBRARIES} #Log4cpp
        ${MPI_CXX_LIBRARIES} #MPI
        ${CUDA_LIBRARIES}    #Cuda
        ${OPENGL_LIBRARIES}  #OpenGl
        ${GLEW_LIBRARIES}    #Glew
        ${GLUT_glut_LIBRARY} #Glut
        )
endforeach()

#enable_testing()

add_custom_target(run
    COMMAND mpirun -np 1 ${EXECUTABLE_OUTPUT_PATH}/main && ${EXECUTABLE_OUTPUT_PATH}/script/merge_files.py && ${EXECUTABLE_OUTPUT_PATH}/display -f data/boids_all.xyz
    DEPENDS main
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
)
