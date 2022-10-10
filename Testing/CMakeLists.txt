include_directories(${AIRS_INCLUDE_DIRS})
set(AIRS_TESTING_TEMP_DIR ${AIRS_BINARY_DIR}/Testing/Temporary)
set(VTK_TESTING_DIRECTORY "${VTK_DIR}/ExternalData/Testing/")

if(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  set(VTK_LIBS vtkRendering vtkIO)
else(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  set(VTK_LIBS vtkRenderingCore vtkRenderingImage vtkInteractionStyle
      vtkRenderingFreeType vtkIOImage vtkIOMINC
      vtkImagingStencil vtkImagingMath vtkImagingMorphological
      vtkFiltersSources vtksys)
  set(VTK_EXTRA_LIBS vtkIOMPIImage
      vtkRenderingOpenGL vtkRenderingFreeTypeOpenGL
      vtkRenderingOpenGL2 vtkRenderingFreeTypeOpenGL2)
  foreach(TMP_LIB ${VTK_E