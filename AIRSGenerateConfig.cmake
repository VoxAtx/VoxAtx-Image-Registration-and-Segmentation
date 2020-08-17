# Generate the AIRSConfig.cmake file in the build tree. This doesnot
# configure one for installation. The file tells external projects how to use
# AIRS.

# Help store a literal dollar in a string.  CMake 2.2 allows escaped
# dollars but we have to support CMake 2.0.
SET(DOLLAR "$")

#-----------------------------------------------------------------------------
# Settings for the build tree.

EXPORT(TARGETS ${AIRS_LIBRARIES}
  FILE ${AIRS_BINARY_DIR}/AIRSTargets.cmake)

# Set the source dir
SET(AIRS_SOURCE_DIR_CONFIG ${AIRS_SOURCE_DIR})

# The library dependencies file.
SET(AIRS_LIBRARY_DEPENDS_FILE
  ${AIRS_BINARY_DIR}/AIRSLibraryDepends.cmake)

INCLUDE(${CMAKE_ROOT}/Modules/CMakeExportBuildSettings.cmake)

CMAKE_EXPORT_BUILD_SETTINGS(
  ${AIRS_BINARY_DIR}/AIRSBuildSettings.cmake)

# The "use" file.
SET(AIRS_USE_FILE_CONFIG
  ${AIRS_BINARY_DIR}/UseAIRS.cmake)

# The build settings file.
SET(AIRS_BUILD_SETTINGS_FILE_CONFIG
  ${AIRS_BINARY_DIR}/AIRSBuildSettings.cmake)

# The target file
SET(AIRS_TARGET_FILE_CONFIG
  ${AIRS_BINARY_DIR}/AIRSTargets.cmake)

# The library directories.
SET(AIRS_LIBRARY_DIRS_CONFIG ${AIRS_LIBRARY_DIR})

# The kits.
SET(AIRS_KITS_CONFIG ${AIRS_KITS})

# The libraries.
SET(AIRS_LIBRARIES_CONFIG ${AIRS_LIBRARIES})

# The include directories.
SET(AIRS_INCLUDE_DIRS_CONFIG "")
FOREACH