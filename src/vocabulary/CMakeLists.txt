SET(vocabulary_SRCS 
  Vocabulary.cpp
) 

SET(vocabulary_HDRS 
  HierarchicalKMeansClustering.h
  HierarchicalKMeansClustering.hpp
  KMeansClustering.h
  KMeansClustering.hpp
  Vocabulary.h
) 

ADD_LIBRARY(gflip_vocabulary SHARED ${vocabulary_SRCS})
TARGET_LINK_LIBRARIES(gflip_vocabulary ${FLIRT_LIBRARIES})

install(TARGETS gflip_vocabulary
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(FILES ${vocabulary_HDRS} DESTINATION include/${PROJECT_NAME}/vocabulary)
