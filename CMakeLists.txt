cmake_minimum_required(VERSION 2.8.12)
project(Snaketest)

SET(SNAKE_SRC
  interpcont.cpp
  extForce.cpp
  imDev.cpp
  snake2D.cpp
  interF.cpp
  GVFimF.cpp
)

ADD_LIBRARY(mysnake
  ${SNAKE_SRC}
)

TARGET_LINK_LIBRARIES(mysnake ${LIBIGL_LIBRARIES} ${LIBIGL_VIEWER_EXTRA_LIBRARIES} ${LIBIGL_TRIANGLE_EXTRA_LIBRARIES} ${LIBIGL_PNG_EXTRA_LIBRARIES})

add_executable(${PROJECT_NAME}_bin main.cpp)
target_include_directories(${PROJECT_NAME}_bin PRIVATE ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(${PROJECT_NAME}_bin PRIVATE ${LIBIGL_DEFINITIONS})
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_VIEWER_EXTRA_LIBRARIES} ${LIBIGL_TRIANGLE_EXTRA_LIBRARIES} ${LIBIGL_PNG_EXTRA_LIBRARIES} mysnake)


