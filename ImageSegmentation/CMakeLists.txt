#PROJECT(ImageSegmentation)

SET(KIT ImageSegmentation)
SET(UKIT IMAGESEGMENTATION)
IF(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(KIT_LIBS vtkHybrid)
ELSE(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(KIT_LIBS vtkImagingStencil vtkFiltersModeling vtkFiltersSources)
ENDIF(${VTK_MAJOR_VERSION} VERSION_LESS 6)

# Set your list of sources here.  Do not change the name of the
# Kit_SRCS variable.
SET(Kit_SR