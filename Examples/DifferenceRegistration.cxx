/*=========================================================================

Program:   Atamai Image Registration and Segmentation
Module:    DifferenceRegistration.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

// This example registers two images (assumed to be of the same patient)
// and then subtracts first image from the second.

#include "AIRSConfig.h"

#include <vtkSmartPointer.h>

#include <vtkImageReslice.h>
#include <vtkImageResize.h>
#include <vtkImageSincInterpolator.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkMath.h>
#include <vtkGlobFileNames.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>

#include <vtkMINCImageReader.h>
#include <vtkDICOMImageReader.h>
#include <vtkMNITransformWriter.h>

// optional readers
#ifdef AIRS_USE_DICOM
#define AIRS_USE_NIFTI
#include <vtkNIFTIReader.h>
#include <vtkNIFTIWriter.h>
#include <vtkDICOMReader.h>
#include <vtkDICOMFileSorter.h>
#include <vtkDICOMMetaData.h>
#endif

#include <vtkRenderer.h>
#include <vtkCamera.h>
#inclu