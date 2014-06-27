#David Coeurjolly <david.coeurjolly@liris.cnrs.fr>
# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(MyDGtalTools)

cmake_minimum_required (VERSION 2.6)


FIND_PACKAGE(DGtal 0.5 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

# -----------------------------------------------------------------------------
# Define variables and cmake parameters, and display some information
# -----------------------------------------------------------------------------
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
SET(DGtalTools_VERSION_MAJOR 0)
SET(DGtalTools_VERSION_MINOR 1)
SET(DGtalTools_VERSION_PATCH 0)
SET(DGTALTOOLS_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(PROJECT_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(VERSION ${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH})



SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
  SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies}
     ${Boost_LIBRAIRIES}
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )

#------Tools subdirectories ------

OPTION(VERBOSE "Verbose" ON)
IF(VERBOSE)
    ADD_DEFINITIONS(-DVERBOSE)
ENDIF(VERBOSE)

SET(DGTAL_TOOLS_SRC
  volFillExterior)


FOREACH(FILE ${DGTAL_TOOLS_SRC})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})
  install(TARGETS ${FILE}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)
