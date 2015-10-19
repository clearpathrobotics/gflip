set(gflip_SRCS
  gflip_engine.cpp
)

set(gflip_HDRS
  gflip_engine.hpp
)

add_library(gflip_gflip SHARED ${gflip_SRCS})
target_link_libraries(gflip_gflip ${Boost_LIBRARIES})

install(TARGETS gflip_gflip
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(FILES ${gflip_HDRS} DESTINATION include/${PROJECT_NAME}/gflip)
