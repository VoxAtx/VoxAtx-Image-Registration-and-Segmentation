PROJECT(AIRS)

# CMake version requirements
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# For CMake 3.0, prefer rpath over loader_path on OS X
IF(POLICY CMP0042)
  CMAKE_POLICY(SET CMP0042 NEW)
ENDIF(POLICY CMP0042)

#SET(CMAKE_SKIP_RPATH ON CACHE BOOL "rpath not included in .so files" FORCE)

# Look for a VTK build directory
FIND_PACKAGE(VTK REQUIRED)
IF (VTK_FOUND)
  INCLUDE (${VTK_USE_FILE})
ENDIF (VTK_FOUND)

# Look for DICOM add-ons to VTK
OPTION(A