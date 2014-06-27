cmake_minimum_required(VERSION 2.6.2)

include_directories(${Simplevol_SOURCE_DIR}/libvol/src)
link_directories (${Simplevol_SOURCE_DIR}/libvol/src)

include_directories(${Simplevol_SOURCE_DIR}/liblongvol/src)
link_directories (${Simplevol_SOURCE_DIR}/liblongvol/src)

SET(optionalLibs ${optionalLibs} vol longvol )

ADD_EXECUTABLE(longvol2sphere longvol2sphere)
TARGET_LINK_LIBRARIES(longvol2sphere ${optionalLibs})
# INSTALL(TARGETS longvol2sphere DESTINATION bin)

#FOREACH(file ${SOURCES_TOOLS})
# ADD_EXECUTABLE(${file} ${file})
# TARGET_LINK_LIBRARIES(${file} ${optionalLibs})
# INSTALL(TARGETS ${file} DESTINATION bin)
#ENDFOREACH(file)


