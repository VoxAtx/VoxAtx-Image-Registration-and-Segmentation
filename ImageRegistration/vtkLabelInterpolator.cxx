
/*=========================================================================
  Program:   Atamai Image Registration and Segmentation
  Module:    vtkLabelInterpolator.h

  Copyright (c) 2014 David Gobbi
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  * Neither the name of the Calgary Image Processing and Analysis Centre
    (CIPAC), the University of Calgary, nor the names of any authors nor
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=========================================================================*/

#include "vtkLabelInterpolator.h"
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
#define VTK_LABEL_KERNEL_TABLE_DIVISIONS 256

vtkStandardNewMacro(vtkLabelInterpolator);

//----------------------------------------------------------------------------
vtkLabelInterpolator::vtkLabelInterpolator()
{
  this->RadiusFactors[0] = 3;
  this->RadiusFactors[1] = 3;
  this->RadiusFactors[2] = 3;
  this->KernelLookupTable[0] = NULL;
  this->KernelLookupTable[1] = NULL;
  this->KernelLookupTable[2] = NULL;
  this->KernelSize[0] = 6;
  this->KernelSize[1] = 6;
  this->KernelSize[2] = 6;
  this->Antialiasing = 0;
  this->BlurFactors[0] = 1.0;
  this->BlurFactors[1] = 1.0;
  this->BlurFactors[2] = 1.0;
  this->LastBlurFactors[0] = 1.0;
  this->LastBlurFactors[1] = 1.0;
  this->LastBlurFactors[2] = 1.0;
}

//----------------------------------------------------------------------------
vtkLabelInterpolator::~vtkLabelInterpolator()
{
  if (this->KernelLookupTable[0])
    {
    this->FreeKernelLookupTable();
    }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "RadiusFactors: " << this->RadiusFactors[0] << " "
     << this->RadiusFactors[1] << " " << this->RadiusFactors[2] << "\n";
  os << indent << "BlurFactors: " << this->BlurFactors[0] << " "
     << this->BlurFactors[1] << " " << this->BlurFactors[2] << "\n";
  os << indent << "Antialiasing: "
     << (this->Antialiasing ? "On\n" : "Off\n");
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::ComputeSupportSize(
  const double matrix[16], int size[3])
{
  // compute the default support size for when matrix is null
  if (this->Antialiasing)
    {
    size[0] = VTK_LABEL_KERNEL_SIZE_MAX;
    size[1] = VTK_LABEL_KERNEL_SIZE_MAX;
    size[2] = VTK_LABEL_KERNEL_SIZE_MAX;
    }
  else
    {
    for (int i = 0; i < 3; i++)
      {
      // use blur factors to compute support size
      size[i] = 2*static_cast<int>(
        this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
      double rowscale = this->BlurFactors[i];
      if (rowscale > (1.0 + VTK_INTERPOLATE_FLOOR_TOL))
        {
        size[i] = 2*static_cast<int>(
          rowscale*this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
        }
      }
    }

  if (matrix == NULL)
    {
    return;
    }

  if (this->Antialiasing)
    {
    // if antialiasing is on, initialize blur factors to 1
    for (int i = 0; i < 3; i++)
      {
      this->BlurFactors[i] = 1.0;
      this->KernelSize[i] = 2*static_cast<int>(
        this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
      }
    }
  else
    {
    // keep blur factors, use kernel size computed from blur factors
    this->KernelSize[0] = size[0];
    this->KernelSize[1] = size[1];
    this->KernelSize[2] = size[2];
    }

  // if matrix does perspective, use the defaults just computed
  if (matrix[12] != 0 || matrix[13] != 0 || matrix[14] != 0 ||
      matrix[15] != 1.0)
    {
    return;
    }

  // use matrix to compute blur factors and kernel size
  for (int i = 0; i < 3; i++)
    {
    double rowscale = 0.0;
    for (int j = 0; j < 3; j++)
      {
      // compute the scale from a row of the matrix
      double x = matrix[4*i + j];
      rowscale += x*x;

      // verify that the element is an integer:
      // check fraction that remains after floor operation
      double f;
      vtkInterpolationMath::Floor(x, f);
      }

    if (this->Antialiasing)
      {
      // rowscale is the subsampling factor in a particular direction
      rowscale = sqrt(rowscale);
      }
    else
      {
      // ignore computed value, use factor provided by SetBlurFactors()
      rowscale = this->BlurFactors[i];
      }

    // if scale is greater than one, expand kernel size
    if (rowscale > (1.0 + VTK_INTERPOLATE_FLOOR_TOL))
      {
      // need extra suport for antialiasing
      this->BlurFactors[i] = rowscale;
      int s = 2*static_cast<int>(
        rowscale*this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
      size[i] = s;
      this->KernelSize[i] = s;
      }
    }

  // rebuild the kernel lookup tables
  this->InternalUpdate();
}

//----------------------------------------------------------------------------
bool vtkLabelInterpolator::IsSeparable()
{
  return true;
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::SetRadiusFactors(
  double x, double y, double z)
{
  if (this->RadiusFactors[0] != x ||
      this->RadiusFactors[1] != y ||
      this->RadiusFactors[2] != z)
    {
    this->RadiusFactors[0] = x;
    this->RadiusFactors[1] = y;
    this->RadiusFactors[2] = z;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::SetBlurFactors(double x, double y, double z)
{
  if (this->BlurFactors[0] != x ||
      this->BlurFactors[1] != y ||
      this->BlurFactors[2] != z)
    {
    this->BlurFactors[0] = x;
    this->BlurFactors[1] = y;
    this->BlurFactors[2] = z;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::SetAntialiasing(int val)
{
  val = (val != 0);
  if (this->Antialiasing != val)
    {
    this->Antialiasing = val;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::InternalDeepCopy(
  vtkAbstractImageInterpolator *a)
{
  vtkLabelInterpolator *obj = vtkLabelInterpolator::SafeDownCast(a);
  if (obj)
    {
    this->SetRadiusFactors(obj->RadiusFactors);
    this->SetAntialiasing(obj->Antialiasing);
    if (this->Antialiasing)
      {
      // if blur factors were computed, then don't call "modified"
      obj->GetBlurFactors(this->BlurFactors);
      }
    else
      {
      this->SetBlurFactors(obj->BlurFactors);
      }
    }

  this->KernelSize[0] = 6;
  this->KernelSize[1] = 6;
  this->KernelSize[2] = 6;

  if (this->KernelLookupTable[0])
    {
    this->FreeKernelLookupTable();
    }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::InternalUpdate()
{
  bool blurchange = false;
  int mode = 0;
  int hsize[3];
  for (int i = 0; i < 3; i++)
    {
    static int minsize = 1;
    static int maxsize = VTK_LABEL_KERNEL_SIZE_MAX/2;
    int size = this->KernelSize[i]/2;
    size = ((size > minsize) ? size : minsize);
    size = ((size < maxsize) ? size : maxsize);
    hsize[i] = size;
    blurchange |= (fabs(this->BlurFactors[i] - this->LastBlurFactors[i]) >=
                   VTK_INTERPOLATE_FLOOR_TOL);
    }

  if (this->BlurFactors[0] > 1.0 + VTK_INTERPOLATE_FLOOR_TOL)
    {
    mode |= VTK_INTERPOLATION_WINDOW_XBLUR_MASK;
    }
  if (this->BlurFactors[1] > 1.0 + VTK_INTERPOLATE_FLOOR_TOL)
    {
    mode |= VTK_INTERPOLATION_WINDOW_YBLUR_MASK;
    }
  if (this->BlurFactors[2] > 1.0 + VTK_INTERPOLATE_FLOOR_TOL)
    {
    mode |= VTK_INTERPOLATION_WINDOW_ZBLUR_MASK;
    }

  mode |= (hsize[0] << VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT);
  mode |= (hsize[1] << VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT);
  mode |= (hsize[2] << VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT);

  if (this->InterpolationInfo->InterpolationMode != mode ||