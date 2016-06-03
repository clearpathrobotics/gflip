cmake_minimum_required(VERSION 2.8.3)
project(gflip_engine)

find_package(catkin REQUIRED COMPONENTS gflip_vocabulary cmake_modules)
find_package(FLIRT REQUIRED)
FIND_PACKAGE(Boost REQUIRED COMPONENTS program_options serialization)
find_package(Eigen REQUIRED)

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${FLIRT_INCLUDE_DIRS}
)
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS
    gflip_vocabulary
  DEPENDS
    Boost
    FLIRT
)

add_library(${PROJECT_NAME}
    src/gflip_engine.cpp
    include/gflip_engine/gflip_engine.hpp
)
add_dependencies(${PROJECT_NAME} ${catkin_EXPORTED_TARGETS})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})

# Build Examples
ADD_EXECUTABLE(featureExtractor examples/FeatureExtractor.cpp)
TARGET_LINK_LIBRARIES(featureExtractor ${PROJECT_NAME} ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(learnVocabularyKMeans examples/LearnVocabularyKMeans.cpp)
TARGET_LINK_LIBRARIES(learnVocabularyKMeans ${PROJECT_NAME} ${catkin_LIBRARIES} ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(generateBoW examples/GenerateBoW.cpp)
TARGET_LINK_LIBRARIES(generateBoW ${catkin_LIBRARIES} ${PROJECT_NAME} ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(nnLoopClosingTest examples/NNLoopClosingTest.cpp)
TARGET_LINK_LIBRARIES(nnLoopClosingTest ${PROJECT_NAME} ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(generateNN examples/GenerateNN.cpp)
TARGET_LINK_LIBRARIES(generateNN ${catkin_LIBRARIES} ${PROJECT_NAME} ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(GFPLoopClosingTest examples/GFPLoopClosingTest.cpp)
add_dependencies(${PROJECT_NAME} ${catkin_EXPORTED_TARGETS})
TARGET_LINK_LIBRARIES(GFPLoopClosingTest ${catkin_LIBRARIES} ${PROJECT_NAME} ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(gflip_cl examples/gflip_cl.cpp)
TARGET_LINK_LIBRARIES(gflip_cl ${catkin_LIBRARIES} ${PROJECT_NAME})

ADD_EXECUTABLE(gflip_cl_onequery examples/gflip_cl_onequery.cpp)
TARGET_LINK_LIBRARIES(gflip_cl_onequery ${catkin_LIBRARIES} ${PROJECT_NAME})



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
