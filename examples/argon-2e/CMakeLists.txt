Builder("-im" "MODES_ENABLED=MODE::IM")
Builder("-re" "MODES_ENABLED=MODE::RE")
Builder("" "")

# Builder("-im" "MODES_ENABLED=1")
# Builder("-re" "MODES_ENABLED=2")

# cmake_minimum_required(VERSION 3.8)
# project(nitrogen-3e LANGUAGES CXX)

# set(CMAKE_VERBOSE_MAKEFILE ON)
# set(CMAKE_RULE_MESSAGES OFF)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
# set(CMAKE_BUILD_TYPE Release)
# add_compile_definitions(PROJNAME=${PROJECT_NAME} DEBUG=223)

# # im + re version
# add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
# target_compile_options(${PROJECT_NAME} PUBLIC -march=native -O3 -ffast-math) 
# target_link_libraries(${PROJECT_NAME} PUBLIC cxxopts::cxxopts)
# target_link_libraries(${PROJECT_NAME} PUBLIC QSF)
# set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "qsf-${PROJECT_NAME}")

# # im version and re version seperately
# add_executable(${PROJECT_NAME}-im ${PROJECT_NAME}.cpp)
# add_executable(${PROJECT_NAME}-re ${PROJECT_NAME}.cpp)
# target_compile_definitions(${PROJECT_NAME}-im PRIVATE MODES_ENABLED=MODE::IM)
# target_compile_definitions(${PROJECT_NAME}-re PRIVATE MODES_ENABLED=MODE::RE)
# target_compile_options(${PROJECT_NAME}-im PUBLIC -march=native -O3 -ffast-math) 
# target_compile_options(${PROJECT_NAME}-re PUBLIC -march=native -O3 -ffast-math) 
# target_link_libraries(${PROJECT_NAME}-im PUBLIC cxxopts::cxxopts)
# target_link_libraries(${PROJECT_NAME}-re PUBLIC cxxopts::cxxopts)
# target_link_libraries(${PROJECT_NAME}-im PUBLIC QSF)
# target_link_libraries(${PROJECT_NAME}-re PUBLIC QSF)
# set_target_properties(${PROJECT_NAME}-im PROPERTIES OUTPUT_NAME "qsf-${PROJECT_NAME}-im")
# set_target_properties(${PROJECT_NAME}-re PROPERTIES OUTPUT_NAME "qsf-${PROJECT_NAME}-re")

# message(STATUS "${Yellow}Binaries for project [${Green}${PROJECT_NAME}${Yellow}] will be placed in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}${ColourReset}")


# install(TARGETS ${PROJECT_NAME}-re DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
# install(TARGETS ${PROJECT_NAME}-im ${PROJECT_NAME}-re ${PROJECT_NAME} 
# DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

# target_compile_options(${PROJECT_NAME}
# PUBLIC -march=native -g -O0 -W -Wall -Wextra -pedantic -Wno-unknown-pragmas -Wno-unused-parameter -Wno-gnu-zero-variadic-macro-arguments -Wno-pedantic) 
# PUBLIC -march=native -O3 -ffast-math -Wno-unknown-pragmas -Wno-unused-parameter -Wno-gnu-zero-variadic-macro-arguments -Wno-pedantic) 

# PUBLIC -O3 -W -Wall -Wextra -pedantic -Wno-unknown-pragmas -Wno-unused-parameter -Wno-gnu-zero-variadic-macro-arguments -Wno-pedantic) 

# target_link_libraries(${PROJECT_NAME} PUBLIC QSF)
# set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "qsf-${PROJECT_NAME}")


# install(TARGETS ${PROJECT_NAME} DESTINATION  ${CMAKE_CURRENT_SOURCE_DIR})

