PROJECT(mciExample CXX)

FIND_PACKAGE(ITK REQUIRED COMPONENTS ITKMorphologicalContourInterpolation ITKSmoothing ITKIOImageBase)
INCLUDE(${ITK_USE_FILE})
  
ADD_EXECUTABLE(mciExample mciExample.cpp)

TARGET_LINK_LIBRARIES(mciExample ${ITK_LIBRARIES})
