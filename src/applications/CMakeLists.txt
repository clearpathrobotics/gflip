ADD_EXECUTABLE(featureExtractor FeatureExtractor.cpp)
TARGET_LINK_LIBRARIES(featureExtractor ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})
ADD_DEPENDENCIES(featureExtractor flirt)

ADD_EXECUTABLE(learnVocabularyKMeans LearnVocabularyKMeans.cpp)
TARGET_LINK_LIBRARIES(learnVocabularyKMeans ${FLIRT_LIBRARIES}  vocabulary ${Boost_LIBRARIES})
ADD_DEPENDENCIES(learnVocabularyKMeans flirt)
##SET_SOURCE_FILES_PROPERTIES(LearnVocabularyKMeans.cpp PROPERTIES COMPILE_FLAGS "-O0 -ggdb ")

ADD_EXECUTABLE(generateBoW GenerateBoW.cpp)
TARGET_LINK_LIBRARIES(generateBoW vocabulary ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})
ADD_DEPENDENCIES(generateBoW flirt)

ADD_EXECUTABLE(nnLoopClosingTest NNLoopClosingTest.cpp)
TARGET_LINK_LIBRARIES(nnLoopClosingTest ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})
ADD_DEPENDENCIES(nnLoopClosingTest flirt)
#SET_SOURCE_FILES_PROPERTIES(NNLoopClosingTest.cpp PROPERTIES COMPILE_FLAGS "-O0 -ggdb ")

ADD_EXECUTABLE(generateNN GenerateNN.cpp)
TARGET_LINK_LIBRARIES(generateNN vocabulary ${FLIRT_LIBRARIES} ${Boost_LIBRARIES} gflip )
ADD_DEPENDENCIES(generateNN flirt)
#SET_SOURCE_FILES_PROPERTIES(NNLoopClosingTest.cpp PROPERTIES COMPILE_FLAGS "-O0 -ggdb ")

ADD_EXECUTABLE(GFPLoopClosingTest GFPLoopClosingTest.cpp)
TARGET_LINK_LIBRARIES(GFPLoopClosingTest vocabulary ${FLIRT_LIBRARIES} ${Boost_LIBRARIES} gflip )
ADD_DEPENDENCIES(GFPLoopClosingTest flirt)
#SET_SOURCE_FILES_PROPERTIES(NNLoopClosingTest.cpp PROPERTIES COMPILE_FLAGS "-O0 -ggdb ")

ADD_EXECUTABLE(gflip_cl gflip_cl.cpp)
TARGET_LINK_LIBRARIES(gflip_cl gflip)

ADD_EXECUTABLE(gflip_cl_onequery gflip_cl_onequery.cpp)
TARGET_LINK_LIBRARIES(gflip_cl_onequery gflip)

install(TARGETS featureExtractor learnVocabularyKMeans generateBoW nnLoopClosingTest generateNN GFPLoopClosingTest gflip_cl gflip_cl_onequery
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
