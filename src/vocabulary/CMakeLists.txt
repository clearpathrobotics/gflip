set(vocabulary_SRCS
  Vocabulary.cpp
)

set(vocabulary_HDRS
  HierarchicalKMeansClustering.h
  HierarchicalKMeansClustering.hpp
  KMeansClustering.h
  KMeansClustering.hpp
  Vocabulary.h
)

add_library(gflip_vocabulary SHARED ${vocabulary_SRCS})
target_link_libraries(gflip_vocabulary ${FLIRT_LIBRARIES} ${Boost_LIBRARIES})

install(TARGETS gflip_vocabulary
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(FILES ${vocabulary_HDRS} DESTINATION include/${PROJECT_NAME}/vocabulary)
