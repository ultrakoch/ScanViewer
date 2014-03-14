cmake_minimum_required (VERSION 2.8)
include(FindPkgConfig)

############################
# General project settings #
############################

PROJECT(hokuyoEW)
SET(APPLICATIONS_VERSION_MAJOR 0)
SET(APPLICATIONS_VERSION_MINOR 1)

#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
PKG_CHECK_MODULES(XML libxml++-2.6)
find_package(PCL 1.2 REQUIRED)



INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/include 
                    /usr/local/include/
                    /usr/include/qt4/
                    /usr/include/lua5.1
         $ENV{OBVIOUSLY_ROOT}
                    ../obcore/math
                    ../obcore/filter
                    ../obcore
         ${PCL_INCLUDE_DIRS}
  )

LINK_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/lib
  $ENV{OBVIOUSLY_ROOT}/build/release/obcore  
  $ENV{OBVIOUSLY_ROOT}/build/release/obvision
  $ENV{OBVIOUSLY_ROOT}/build/release/obgraphic
  $ENV{OBVIOUSLY_ROOT}/build/release/obdevice
  ${PCL_LIBRARY_DIRS}
  )

ADD_LIBRARY(laserdev src/UTM-30LX_EW.cpp)


add_definitions(${PCL_DEFINITIONS})


SET(URGLIBS urg_c urg_cpp rt)
SET(CORELIBS obcore gsl gslcblas lua5.1)
SET(DEVICELIBS obdevice laserdev)
SET(GRAPHICLIBS obgraphic vtkHybrid glut GL jpeg)
SET(VISIONLIBS obvision ann flann)


ADD_EXECUTABLE(hokuyoSaveTest src/hokuyoSaveTest.cpp)
ADD_EXECUTABLE(hokuyoShow src/hokuyoShow.cpp)
ADD_EXECUTABLE(fileShow src/fileShow.cpp)
#ADD_EXECUTABLE(hokuyoClassTest src/hokuyoClassTest.cpp)

TARGET_LINK_LIBRARIES(hokuyoSaveTest ${URGLIBS} ${CORELIBS} ${VISIONLIBS} ${GRAPHICLIBS})
TARGET_LINK_LIBRARIES(hokuyoShow ${URGLIBS} ${GRAPHICLIBS} ${CORELIBS} ${VISIONLIBS} ${DEVICELIBS})
TARGET_LINK_LIBRARIES(fileShow ${URGLIBS} ${GRAPHICLIBS} ${CORELIBS} ${VISIONLIBS} ${DEVICELIBS} ${PCL_LIBRARIES})
#TARGET_LINK_LIBRARIES(hokuyoClassTest ${DEVICELIBS} ${URGLIBS} ${GRAPHICLIBS} ${CORELIBS} ${VISIONLIBS} )



