# - Try to find Library grd
# Once done, this will define
#
#  grd_FOUND - system has grd module
#  grd_INCLUDE_DIRS - the grd include directories
#  grd_LIBRARY_DIRS - the grd library directories
#  grd_LIBRARIES - link these to use grd


# Uses  directory to search mrf_segmentation directory!
set(grd_PREFIX_DIR /usr/local)
message(STATUS "Searching grd in directory ${grd_PREFIX_DIR}." )

# Searches include directory /usr/local/include/grd
find_path(grd_INCLUDE_DIR grd ${grd_PREFIX_DIR}/include)
message(STATUS "    grd_INCLUDE_DIR ${grd_INCLUDE_DIR}." )
set(grd_INCLUDE_DIRS ${grd_INCLUDE_DIR})
  
# Searches library librimagraph.a in /usr/local/lib
find_path(grd_LIBRARY_DIR libgrd.a ${grd_PREFIX_DIR}/lib)
message(STATUS "    grd_LIBRARY_DIR ${grd_LIBRARY_DIR}." )
set(grd_LIBRARY_DIRS ${grd_PREFIX_DIR}/lib)

# Sets the names of library components (actually A name and A component)
find_library(grd_LIBRARY grd ${grd_LIBRARY_DIRS})
message(STATUS "    grd_LIBRARY ${grd_LIBRARY}." )
set(grd_LIBRARIES ${grd_LIBRARY})

if(("${grd_INCLUDE_DIR}" STREQUAL "grd_INCLUDE_DIR-NOTFOUND") OR
   ("${grd_LIBRARY_DIRS}" STREQUAL "grd_LIBRARY_DIRS-NOTFOUND") OR
   ("${grd_LIBRARY}" STREQUAL "grd_LIBRARY-NOTFOUND")
  )
  message(STATUS "Library grd NOT found")
  unset(grd_FOUND)
  unset(grd_INCLUDE_DIR)
  unset(grd_LIBRARY_DIR)
  unset(grd_LIBRARY)
  unset(grd_LIBRARIES)
endif()

mark_as_advanced(grd_INCLUDE_DIRS grd_LIBRARY_DIRS grd_LIBRARIES)

