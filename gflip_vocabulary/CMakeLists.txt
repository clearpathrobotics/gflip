cmake_minimum_required(VERSION 2.8.3)
project(gflip_vocabulary)

add_compile_options(-std=c++14 -Wall -Wextra -Werror)

find_package(catkin REQUIRED)
find_package(Boost REQUIRED COMPONENTS random serialization)
find_package(FLIRT REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  DEPENDS
    Boost
    FLIRT
)

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${FLIRT_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
    src/vocabulary.cpp
)
add_dependencies(${PROJECT_NAME} ${catkin_EXPORTED_TARGETS})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})

install(
  TARGETS
    ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h*"
)
