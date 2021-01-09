
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGaussianInterpolator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkGaussianInterpolator.h"
#include "vtkImageInterpolatorInternals.h"
#include "vtkImageData.h"
#include "vtkDataArray.h"
#include "vtkObjectFactory.h"

#include "vtkTemplateAliasMacro.h"
// turn off 64-bit ints when templating over all types, because
// they cannot be faithfully represented by doubles
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

// masks for storing window and size in a single integer
#define VTK_INTERPOLATION_WINDOW_MASK        0x0000007f
#define VTK_INTERPOLATION_WINDOW_XBLUR_MASK  0x00008000
#define VTK_INTERPOLATION_WINDOW_XSIZE_MASK  0x00007f00
#define VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT 8
#define VTK_INTERPOLATION_WINDOW_YBLUR_MASK  0x00800000
#define VTK_INTERPOLATION_WINDOW_YSIZE_MASK  0x007f0000
#define VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT 16
#define VTK_INTERPOLATION_WINDOW_ZBLUR_MASK  0x80000000
#define VTK_INTERPOLATION_WINDOW_ZSIZE_MASK  0x7f000000
#define VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT 24

// kernel lookup table size must be 256*n where n is kernel half-width
// in order to provide sufficient precision for 16-bit images
#define VTK_GAUSS_KERNEL_TABLE_DIVISIONS 256

vtkStandardNewMacro(vtkGaussianInterpolator);
