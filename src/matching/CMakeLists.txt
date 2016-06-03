set(matching_SRCS
  base_matcher.h
  base_matcher.cpp
  geometrical_phrases_matcher.h
  geometrical_phrases_matcher.cpp
  standard_bow_matcher.h
  standard_bow_matcher.cpp
  laser_scan_data.h
  laser_scan_data.cpp
)

set(matching_HDRS
  base_matcher.h
  geometrical_phrases_matcher.h
  standard_bow_matcher.h
  laser_scan_data.h
)

add_library(gflip_matching SHARED ${matching_SRCS})
target_link_libraries(gflip_matching ${Boost_LIBRARIES})


install(TARGETS gflip_matching
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(FILES ${matching_HDRS} DESTINATION include/${PROJECT_NAME}/matching)
