cmake_minimum_required(VERSION 3.8)
project(QSF-examples)
# Prepare the library
add_subdirectory(QSF)
add_subdirectory(extern/cxxopts)
# Compile example projects
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(colors)
define_colors()




macro(Builder mod flags)
    get_filename_component(prj ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    string(REPLACE " " "_" prj ${prj})
    cmake_minimum_required(VERSION 3.8)
    project("${prj}${mod}" LANGUAGES CXX)

    set(CMAKE_VERBOSE_MAKEFILE ON)
    set(CMAKE_RULE_MESSAGES OFF)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    set(CMAKE_BUILD_TYPE Release)
    # im + re version
    add_executable(${prj}${mod} "${prj}.cpp")
    target_compile_options(${prj}${mod} PUBLIC -g -O3 -ffast-math $<$<CXX_COMPILER_ID:GNU>:-march=native>) 
    target_link_libraries(${prj}${mod} PUBLIC cxxopts::cxxopts)
    target_link_libraries(${prj}${mod} PUBLIC QSF)
    set_target_properties(${prj}${mod} PROPERTIES OUTPUT_NAME "qsf-${prj}${mod}")
    target_compile_definitions(${prj}${mod} PRIVATE ${flags} PROJNAME=${prj} DEBUG=223)
    message(STATUS "${Yellow}Binaries for project [${Green}${prj}${mod}${Yellow}] will be placed in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}${ColourReset}")
endmacro()

add_subdirectory(examples/template-ini)
add_subdirectory(examples/nitrogen-3e)
add_subdirectory(examples/neon-3e)
add_subdirectory(examples/argon-2e)
add_subdirectory(examples/thiede)




# set(CMAKE_FIND_FRAMEWORK NEVER)
# set(CMAKE_FIND_APPBUNDLE NEVER)
# set(CMAKE_FIND_USE_CMAKE_SYSTEM_PATH FALSE)
# set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH FALSE)
# set(CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}")
# project(template LANGUAGES CXX)

# set(default_build_type "Release")
# if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
#   message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
#   set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
#       STRING "Choose the type of build." FORCE)
#   # Set the possible values of build type for cmake-gui
#   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
#     "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
# endif()


# detect operating system and host processor
# message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
# message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")

# option(ENABLE_MPI "Build with MPI" ON)
# option(BUILD_SUBMODES "Build targets with limited MODES (IM/RE)" ON)
# option(BUILD_SUBDIMS "Build targets with limited DIMS (dimensions)" ON)
# You can change the Prefix, Output Name and Suffix using the set_target_properties() function and the PREFIX / OUTPUT_NAME / SUFFIX property in the following way:
# set_target_properties(new_thing PROPERTIES PREFIX "")
# set_target_properties(new_thing PROPERTIES OUTPUT_NAME "better_name")
# set_target_properties(new_thing PROPERTIES SUFFIX ".so.1")
# set_target_properties(new_thing PROPERTIES PREFIX "" SUFFIX ".so.${NEW_THING_VER}" )
# set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")
# add_executable(${PROJECT_NAME} ${PROJECT_NAME}/${PROJECT_NAME}.cpp)
# find_library (FFTW_LIBRARIES NAMES fftw3)
# target_link_libraries(${PROJECT_NAME} PUBLIC ${FFTW_LIBRARIES}) 



# set(CMAKE_OSX_SYSROOT /opt/homebrew/)
# set(CMAKE_SHARED_LINKER_FLAGS -lfftw3_mpi,-lfftw3)
# set(CMAKE_SHARED_LINKER_FLAGS -Wl,-undefined,dynamic_lookup)
# find_package(MPI REQUIRED)
# if (MPI_FOUND)
#     if (MPI_C)
#         message(STATUS "MPI ${MPI_VERSION} found with MPI_C bindings")
# target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_C)
#     endif()
#     if (MPI_CXX)
#         message(STATUS "MPI ${MPI_VERSION} found with MPI_CXX bindings")
# target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
#     endif()
#     # message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
# else()
# message(STATUS "MPI not found :(")
# endif()





# target_compile_definitions(
    # ${PROJECT_NAME} PRIVATE
    # $<$<BOOL:${YOUR_DEFINE}>:YOUR_DEFINE>)

    # find_package(Git)
# execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                # WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                # OUTPUT_VARIABLE PACKAGE_GIT_VERSION
                # ERROR_QUIET
                # OUTPUT_STRIP_TRAILING_WHITESPACE)
# include(CheckCXXCompilerFlag)
# add_compile_options(-march=native)
# target_compile_definitions(${PROJECT_NAME}
#   PRIVATE
    # DIM=3
    # DEBUG=1
#     USE_SIMD=1
#     INTERNAL_MACRO
#   INTERFACE
#     ONLY_CONSUMERS_CAN_SEE_THIS_DEFINE=42
# )
# if(ENABLE_MPI)
# set(CMAKE_LIBRARY_PATH /opt/homebrew/lib)
    # target_link_libraries(MyTarget PUBLIC MPI::MPI_CXX)
# endif()
# target_link_libraries(QSF INTERFACE ${FFTW_LIBRARIES})
# target_link_libraries(${PROJECT_NAME} PUBLIC QSF)
# target_compile_options(${PROJECT_NAME} PRIVATE -lfftw -lm)
# include(GNUInstallDirs)
# set(CMAKE_SKIP_BUILD_RPATH FALSE)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
# if("${isSystemDir}" STREQUAL "-1")
#     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# endif("${isSystemDir}" STREQUAL "-1")

# find_package(PkgConfig REQUIRED)
# pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
# target_include_directories(${PROJECT_NAME} PkgConfig::FFTW)
# target_compile_options(${PROJECT_NAME} PUBLIC -lm -lmpi -lfftw3_mpi -lfftw3)
# target_link_libraries(${PROJECT_NAME} PUBLIC PkgConfig::FFTW)


# target_include_directories(${PROJECT_NAME} PUBLIC ${FFTW3_INCLUDE_DIRS})
# target_link_libraries(${PROJECT_NAME} PUBLIC ${FFTW3_LIBRARIES})
# target_link_libraries(${PROJECT_NAME} PUBLIC ${FFTW3_LIBRARY_DIRS})
# target_compile_options(${PROJECT_NAME} PUBLIC -lm -lmpi -lfftw3_mpi -lfftw3)
# target_link_libraries(${PROJECT_NAME} PUBLIC fftw3)
# target_link_libraries(${PROJECT_NAME} PUBLIC FFTW3)
# install(TARGETS ${PROJECT_NAME} DESTINATION  ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME})



# ADD LIBRARIES HERE. EXAMPLES:
# # Header Only:
# find_package  (GLM REQUIRED)
# import_library(GLM_INCLUDE_DIRS)
# # Identical Debug and Release:
# find_package  (OpenGL REQUIRED)
# import_library(OPENGL_INCLUDE_DIR OPENGL_LIBRARIES)
# # Separate Debug and Release:
# find_package  (Boost REQUIRED iostreams)
# import_library(Boost_INCLUDE_DIRS Boost_IOSTREAMS_LIBRARY_DEBUG Boost_IOSTREAMS_LIBRARY_RELEASE)


# set_property(GLOBAL PROPERTY USE_FOLDERS ON) # Turn on this feature
# set_property(TARGET my_lib PROPERTY FOLDER "libs")
# set_property(TARGET another_lib PROPERTY FOLDER "libs/3rdParty")
# set_property(TARGET my_app PROPERTY FOLDER "utils")

# ${PROJECT_NAME}

# add_compile_options(
#        -Wall -Werror -Wno-error=maybe-uninitialized
#        $<$<CONFIG:RELEASE>:-Ofast>
#        $<$<CONFIG:DEBUG>-O0>
#        $<$<CONFIG:DEBUG>-ggdb3>
# )
# add_compile_definitions(
#         $<$<CONFIG:RELEASE>:NDEBUG>
#         $<$<CONFIG:RELEASE>:BOOST_DISABLE_ASSERTS>
# )

# ::::::::::::::::::::::::::::Include guards ::::::::::::::::::::::::::::::::::

# if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
# message(FATAL_ERROR
#     "Run with: cmake . -B project/build; cmake --build project/build --target install"
# )
# endif()



# INTERFACE versus PUBLIC versus PRIVATE
# If PRIVATE is specified for a certain option/property, then that option/property will only impact the current target. If PUBLIC is specified, then the option/property impacts both the current target and any others that link to it. If INTERFACE is specified, then the option/property does not impact the current target but will propagate to other targets that link to it.
# This makes INTERFACE very useful for header-only libraries that aren't compiled (on their own), but rely on other libraries they need to link to. Using target_link_library on our header-only library with INTERFACE specified means that the executable that eventually gets compiled with our library will have the appropriate link flags.

#Great source: https://github.com/bast/cmake-example/blob/master/CMakeLists.txt
#Use directory name: https://stackoverflow.com/questions/12344368/automatically-use-the-directory-as-the-project-name-in-cmake
#simple tut: https://www.siliceum.com/en/blog/post/cmake_01_cmake-basics
#in depth: https://pabloariasal.github.io/2018/02/19/its-time-to-do-cmake-right/
#nice examples: https://github.com/ttroy50/cmake-examples
#nifty: https://github.com/pr0g/cmake-examples

#Installation strategies (someday maybe)
#source: https://github.com/pabloariasal/modern-cmake-sample#using-cmakes-user-package-registry
# Normal Installation
# sudo cmake --build build -- install
# This will install the example library under /usr/local/ on UNIX systems.
# Alternatively, you can specify a custom installation directory by setting -DCMAKE_INSTALL_PREFIX in the cmake configure step:
# cmake -Bbuild -DCMAKE_INSTALL_PREFIX=<custom_install_dir>
# sudo cmake --build build -- install
# To uninstall the library, you can run:
# cd build
# xargs rm < install_manifest.txt

#Installation to CMake's User Package Registry
# Instead of actually installing the library, you can just build it and register the build in CMake's User Package Registry
# cd libjsonutils
# cmake -Bbuild -DCMAKE_EXPORT_PACKAGE_REGISTRY
#This will register the library's build in CMake's User Package Registry (on UNIX systems it defaults to ~/.cmake).
#This is convenient, as packages depending on the library (e.g. via find_package) will be able to find it through the registry, even when the library hasn't been installed.