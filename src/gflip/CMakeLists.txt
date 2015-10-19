SET(gflip_SRCS 
  gflip_engine.cpp
) 

SET(gflip_HDRS 
  gflip_engine.hpp
) 

ADD_LIBRARY(gflip_gflip SHARED ${gflip_SRCS})

install(TARGETS gflip_gflip
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(FILES ${gflip_HDRS} DESTINATION include/${PROJECT_NAME}/gflip)
