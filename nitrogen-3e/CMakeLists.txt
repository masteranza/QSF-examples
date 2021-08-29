cmake_minimum_required(VERSION 3.8)
project(nitrogen-3e LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_RULE_MESSAGES OFF)

# imaginary-time version
add_executable(${PROJECT_NAME}-im ${PROJECT_NAME}.cpp)
# real-time version
add_executable(${PROJECT_NAME}-re ${PROJECT_NAME}.cpp)
# im + re version
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)



target_compile_definitions(${PROJECT_NAME}-im PRIVATE MODES_ENABLED=MODE::IM)
target_compile_definitions(${PROJECT_NAME}-re PRIVATE MODES_ENABLED=MODE::RE)
add_compile_definitions(PROJNAME=${PROJECT_NAME} DEBUG=223)

# option(BATCH "Build for batch execution" OFF)
# if (${BATCH})
#     add_compile_definitions(BATCH=1)
#     message(WARNING "BATCH set to enabled")
# endif()

# option(TEST "Build for test execution" OFF)
# if (${TEST})
#     add_compile_definitions(TEST=1)
#     message(WARNING "TEST set to enabled")
# endif()

# add_compile_options(
# PUBLIC -march=native -O3 -ffast-math -Wno-unknown-pragmas -Wno-unused-parameter -Wno-gnu-zero-variadic-macro-arguments -Wno-pedantic) 

# target_compile_options(${PROJECT_NAME}-re PUBLIC -march=native -O0 -g) 
# target_compile_options(${PROJECT_NAME}-im PUBLIC -march=native -O0 -g) 
# target_compile_options(${PROJECT_NAME} PUBLIC -march=native -O0 -g) 

target_compile_options(${PROJECT_NAME}-im PUBLIC -march=native -O3 -ffast-math) 
target_compile_options(${PROJECT_NAME}-re PUBLIC -march=native -O3 -ffast-math) 
target_compile_options(${PROJECT_NAME} PUBLIC -march=native -O3 -ffast-math) 

target_link_libraries(${PROJECT_NAME} PUBLIC cxxopts::cxxopts)
target_link_libraries(${PROJECT_NAME}-im PUBLIC cxxopts::cxxopts)
target_link_libraries(${PROJECT_NAME}-re PUBLIC cxxopts::cxxopts)

target_link_libraries(${PROJECT_NAME} PUBLIC QSF)
target_link_libraries(${PROJECT_NAME}-im PUBLIC QSF)
target_link_libraries(${PROJECT_NAME}-re PUBLIC QSF)
set_target_properties(${PROJECT_NAME}-im PROPERTIES OUTPUT_NAME "qsf-${PROJECT_NAME}-im")
set_target_properties(${PROJECT_NAME}-re PROPERTIES OUTPUT_NAME "qsf-${PROJECT_NAME}-re")
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "qsf-${PROJECT_NAME}")

# install(TARGETS ${PROJECT_NAME}-re DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ${PROJECT_NAME}-im ${PROJECT_NAME}-re ${PROJECT_NAME} 
DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

# target_compile_options(${PROJECT_NAME}
# PUBLIC -march=native -g -O0 -W -Wall -Wextra -pedantic -Wno-unknown-pragmas -Wno-unused-parameter -Wno-gnu-zero-variadic-macro-arguments -Wno-pedantic) 
# PUBLIC -march=native -O3 -ffast-math -Wno-unknown-pragmas -Wno-unused-parameter -Wno-gnu-zero-variadic-macro-arguments -Wno-pedantic) 

# PUBLIC -O3 -W -Wall -Wextra -pedantic -Wno-unknown-pragmas -Wno-unused-parameter -Wno-gnu-zero-variadic-macro-arguments -Wno-pedantic) 

# target_link_libraries(${PROJECT_NAME} PUBLIC QSF)
# set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "qsf-${PROJECT_NAME}")


# install(TARGETS ${PROJECT_NAME} DESTINATION  ${CMAKE_CURRENT_SOURCE_DIR})

