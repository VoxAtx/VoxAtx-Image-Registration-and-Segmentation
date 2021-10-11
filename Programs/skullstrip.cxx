/*=========================================================================

Program:   Atamai Image Registration and Segmentation
Module:    skullstrip.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

// Apply Stephen Smith's BET algorithm to an MRI image.

#include <vtkSmartPointer.h>

#include <vtkImageReslice.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkROIStencilSource.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLookupTable.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#incl